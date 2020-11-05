#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 One-Particle Density Matrix Module.
 Bartosz Błasiak, Gundelfingen, May 2019
"""

import oepdev
import numpy
import psi4
import scipy.optimize
from abc import ABC, abstractmethod

__all__ = ["Density"]

class Density:
    """
 --------------------------------------------------------------------------------------------------------------------------
                                                Electron Density

 Handles the Electron Density Distribution.

 --------------------------------------------------------------------------------------------------------------------------

 Usage as container class: 
 
  1) Initialize container object:
     
     density = Density(D = None, jk = None)

     where:
       o D  - the density matrix in AO or MO basis
       o jk - the psi4::JK object for AO basis JK calculations

  2) Grab the density matrix

     D = density.matrix()

  3) Computations in AO basis:

     o compute 1-electron energy (does not require JK object to be set)

       e_1 = density.compute_1el_energy(D, V1)

     The below require jk to be set:

     o compute 2-electron energy (J-type expression)

       e_2j = density.compute_2el_energy(D_left, D_right, type='j')

     o compute 2-electron energy (K-type expression)

       e_2k = density.compute_2el_energy(D_left, D_right, type='k')

     o compute J matrix (or K matrix if type=='k'):

       J = density.generalized_JK(D, type='j')

 --------------------------------------------------------------------------------------------------------------------------

 Usage as method class. 

 Using 'Density' as a class of methods do not require object initialization.
 The list of class methods is given below:

   o Density.natural_orbitals     - compute natural orbitals
   o Density.generalized_density  - compute generalized OPDM
   o Density.orthogonalize_OPDM   - compute orthogonalized OPDM
   o Density.orthogonalizer       - compute orthogonalizer matrix
   o Density.deorthogonalizer     - compute deorthogonalizer rmatrix

 Usage: 
   result = Density.'class method name'

 See respective documentation for each of them for further details.

 --------------------------------------------------------------------------------------------------------------------------
                                                                                  Last Revision: Gundelfingen, 20 May 2019
"""
    def __init__(self, D=None, jk=None):
        self._D = None
        if D is not None:
           self._D         =  D.copy()
        self._global_jk = jk

    def matrix(self):
        return self._D

    def set_D(self, D): self._D = D.copy()

    def set_jk(self, jk):
        self._global_jk = jk


    # ---> Container Public Interface <--- #

    def compute_1el_energy(self, D, Hcore):
        "Compute generalized 1-electron energy"
        energy = numpy.dot(D, Hcore).trace()
        return energy

    def compute_2el_energy(self, D_left, D_right, type='j'):
        "Compute generalized 2-electron energy"
        JorK = self.generalized_JK(D_left, type)
        energy = numpy.dot(JorK, D_right).trace()
        return energy

    def generalized_JK(self, D, type='j'):
        "Compute J or K matrix in AO basis from OPDM in the same AO basis"
        self._global_jk.C_clear()                                           
        self._global_jk.C_left_add(psi4.core.Matrix.from_array(D, ""))
        I = numpy.identity(D.shape[0], numpy.float64)
        self._global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._global_jk.compute()
        if   type.lower() == 'j': JorK = self._global_jk.J()[0].to_array(dense=True)
        elif type.lower() == 'k': JorK = self._global_jk.K()[0].to_array(dense=True)
        else: raise ValueError("Incorrect type of JK matrix. Only J or K allowed.")
        return JorK


    # ---> Class Public Interface <--- #

    @classmethod
    def natural_orbitals(cls, D, 
                              S                      =  None         , 
                              C                      =  None         , 
                              orthogonalize_mo       =  True         , 
                              order                  = 'descending'  , 
                              return_ao_orthogonal   =  False        , 
                              renormalize            =  False        , 
                              no_cutoff              =  False        ,
                              ignore_large_n         =  False        ,
                              n_eps                  =  5.0E-5       ):
        """
 --------------------------------------------------------------------------------------------------------------------------
 Compute the Natural Orbitals from a given OPDM

 --------------------------------------------------------------------------------------------------------------------------

 Usage:

 n, c = Density.natural_orbitals(D, S = None, C = None, 
                                    orthogonalize_mo     = True,
                                    order                = 'descending',
                                    return_ao_orthogonal = False,
                                    renormalize          = False, 
                                    no_cutoff            = False, 
                                    ignore_large_n       = False, 
                                    n_eps                = 5.0E-5)

 where:
  o D - OPDM in AO or MO basis
  o S - overlap integrals in AO or MO basis
  o C - LCAO-MO transformation matrix
  o orthogonalize_mo     - whether to transform D from AO to certain MO basis and diagonalize
  o order                - order in which eigenvalues (occupancies) are sorted. Eigenvalues (NO's) are sorted accordingly.
  o return_ao_orthogonal - whether to return NO's in oAO basis set or not
  o renormalize          - renormalize to integer number of electrons
  o no_cutoff            - cut-off threshold for occupancies
  o ignore_large_n       - raise ValueError if (1.0 + n_eps) < n < (0.0 - n_eps)
  o n_eps                - tolerance for occupancy deviation

 --------------------------------------------------------------------------------------------------------------------------

 Examples:

  1) NO's in AO (non-orthogonal, original) basis from D in AO basis

     n, c = Density.natural_orbitals(D, S, C, orthogonalize_mo = True, n_eps = 0.001)

     D: ndarray of shape (AO x AO)
     S: ndarray of shape (AO x AO)
     C: ndarray of shape (AO x MO)

     --> transformation D (MO x MO) = C.T S D S C and its diagonalization
     --> transformation of transformation matrix from MO to AO basis

     n: ndarray of shape (NO)
     c: ndarray of shape (AO x NO)

  2) NO's in certain orthogonal MO basis from D in the same MO basis
    
     n, c = Density.natural_orbitals(D, None, None, orthogonalize_mo = False, n_eps = 0.001)

     D: ndarray of shape (MO x MO)

     --> diagonalization of D

     n: ndarray of shape (NO)
     c: ndarray of shape (MO x NO)

 --------------------------------------------------------------------------------------------------------------------------
                                                                                  Last Revision: Gundelfingen, 20 May 2019
 """
        # orthogonalize in MO basis
        if orthogonalize_mo is True:
            if S is not None:
                assert C is not None
                D_ = numpy.linalg.multi_dot([C.T, S, D, S, C])
            else:
                D_ = numpy.linalg.multi_dot([C.T,    D,    C])
        # orthogonalize in AO basis
        else:
            if S is not None:
               D_ = cls.orthogonalize_OPDM(D, S)
            else:
               D_ = D

        # Diagonalize density matrix in OAO or MO basis
        n, L = numpy.linalg.eigh(D_)
        #print(n.sum())

        # LCAO_NO matrix
        if orthogonalize_mo is True:
           U = numpy.dot(C, L)
        else:
           U = L
        if return_ao_orthogonal is False:
           if orthogonalize_mo is False:
              if S is not None:
               U = numpy.dot(cls.orthogonalizer(S), U)
              else: U = U
        else:
           raise NotImplementedError("Returning LCoAO_NO matrix is not supported yet if orthogonalization in MO basis was done!")

        # Warnings and sanity checks
        if not ignore_large_n:
           if n.max() > 1.0 or n.min() < 0.0:
              print(" Warning! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))
           if ((n.max() - 1.0) > n_eps or (n.min() < -n_eps)):
              raise ValueError("Unphysical NO populations detected! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))

        # Remove negative values
        n[numpy.where(n<0.0)] = 0.0

        # NO cutoff
        if no_cutoff is False: no_cutoff = self.no_cutoff
        if no_cutoff != 0.0:
           ids = numpy.where(n>=no_cutoff)
           n = n[ids]
           U =(U.T[ids]).T

        # Order according to occupations
        if order.lower() =='ascending': 
           pass
        elif order.lower() =='descending':
           n = n[  ::-1]
           U = U[:,::-1]
        else: raise ValueError("Incorrect order of NO orbitals. Possible only 'ascending' or 'descending'.")

        # Renormalize to match correct number of electrons
        if renormalize is True:
           if ( abs(n.sum() - numpy.round(n.sum())) > 1.e-7):
              print(" Warning: nsum=%14.4E delta=%14.4E" % (n.sum(), n.sum() - numpy.round(n.sum())))
           d = numpy.round(n.sum()) - n.sum()
           d/= numpy.float64(n.size)
           n+= d
           n[numpy.where(n<0.0)] = 0.0
           n[numpy.where(n>1.0)] = 1.0
        return n, U

    @classmethod
    def generalized_density(cls, n, c, g=1.0):
        "Generalized matrix G = C n**g C.T"
        ng = n.copy() if g==1.0 else n**g
        G = numpy.linalg.multi_dot([c, numpy.diag(ng), c.T])
        return G

    @classmethod
    def orthogonalize_OPDM(cls, D, S):
        "Transforms the one-particle density matrix to orthogonal space"
        Y = cls.deorthogonalizer(S)
        return numpy.dot(Y, numpy.dot(D, Y.T))

    @classmethod
    def deorthogonalizer(cls, S):
        "Compute the deorthogonalizer matrix from the overlap matrix"
        s, u = numpy.linalg.eig(S)
        s = numpy.sqrt(s)
        Y = numpy.dot(u, numpy.dot(numpy.diag(s), u.T))
        return Y

    @classmethod
    def orthogonalizer(cls, S):
        "Compute the orthogonalizer matrix from the overlap matrix"
        s, u = numpy.linalg.eig(S)
        sm = 1.0/numpy.sqrt(s)
        X = numpy.dot(u, numpy.dot(numpy.diag(sm), u.T))
        return X


class DensityProjection(ABC):
    """\
 Gradient Projection Algorithms.
 Ref.: Pernal, Cances, J. Chem. Phys. 2005

 Usage:
  proj = DensityProjection.create(np, dtype='p', S=None)
  n, c = proj.compute(n, c, S)
"""
    def __init__(self, np, S):
        self._np = np
        self._S  = S
        super(DensityProjection, self).__init__()


    # --> Public Interface <-- #

    @staticmethod
    def create(np, dtype='p', S=None):
        if   dtype.lower() == 'p': return Pset_DensityProjection(np, S)
        elif dtype.lower() == 'd': return Dset_DensityProjection(np, S)
        else: raise ValueError("Only projections onto D and P sets are possible.")

    def compute(self, n, c, perfect_pairing=False):
        return self._density_matrix_projection(n, c, perfect_pairing)


    # --> Protected Interface <-- #

    def _density_matrix_projection(self, n, c, perfect_pairing):
        "Find n_new and C_new such that new density matrix is N-representable"                              
        if self._S is None: pass # S = numpy.identity(len(n))
        else              : S = self._S.copy()
        # compute pre-density matrix                                                                       
        A = Density.generalized_density(n, c) # cannot be here self.D because it is pre-density matrix!
        #A = numpy.linalg.multi_dot([S, preA, S])
        if self._S is not None: A = Density.orthogonalize_OPDM(A, S)
        #a, b = scipy.linalg.eig(A, S)
        a, phi = numpy.linalg.eigh(A)
        #a = a.real; b = b.real
       #print(" Init sum = %14.6f" % (a**2).sum()) 
                                                                                                            
        muORnu = self._find_coef(a)
                                                                                                           
        # compute the projected density matrix
        n_new = self._eval_coef(a, muORnu)
       #print("OptSum= ", (n_new**2).sum())
        C_new = phi
                                                                                                            
        # sort (descending order)
        idx = numpy.argsort(n_new)[::-1]
        n_new = n_new [  idx]
        C_new = C_new [:,idx]
        if self._S is not None: C_new = numpy.dot(Density.orthogonalizer(S), C_new)

        # perfect-pairing projection (sorting will be changed to pp-descending)
        if perfect_pairing:
           n_new, C_new = self._perfect_pair_projection(n_new, C_new)
        return n_new.real, C_new.real

    @abstractmethod
    def _find_coef(self, n):
        pass

    @abstractmethod
    def _eval_coef(self, x, coeff):
        pass

    @abstractmethod
    def _perfect_pair_projection(self, n, c):
        pass

           
class Dset_DensityProjection(DensityProjection):
    """\
 Gradient Projection Algorithm on D-sets.
 Ref.: Pernal, Cances, J. Chem. Phys. 2005

 Notes:
  o Appropriate only for HF functional.
"""
    def __init__(self, np, S):
        super(Dset_DensityProjection, self).__init__(np, S)

    # --> Implementation <-- #

    def _find_coef(self, n):
        "Search for mu"                                                             
        options = {'disp': False, 'maxiter':1000}
        mu = 0.0
        def obj(mu, x):
            u = self._eval_coef(x, mu)
            Z = (u.sum() - self._np)**2
            return Z
        R = scipy.optimize.minimize(obj, mu, args=(n,), tol=1e-20, options=options)
        mu = R.x
        return mu

    def _eval_coef(self, a, mu):
        "Projected occupation numbers" 
        a_ = a.copy();
        for i in range(len(a)):
            u = a[i] + mu
            if   u <= 0.0: a_[i] = 0.0
            elif u >= 1.0: a_[i] = 1.0
            else: a_[i] = u
        return a_

    def _perfect_pair_projection(self, n, c):
        raise NotImplementedError("Perfect Pairing scheme is not implemented for D-set optimization!")


class Pset_DensityProjection(DensityProjection):
    """\
 Gradient Projection Algorithm on P-sets.
 Ref.: Pernal, Cances, J. Chem. Phys. 2005

 Notes:
  o Appropriate for any DMFT functional.
"""
    def __init__(self, np, S):
        super(Pset_DensityProjection, self).__init__(np, S)

    # --> Implementation <-- #

    def _find_coef(self, n):
        "Search for nu"                                                             
        options = {'disp': False, 'maxiter':2000, 'ftol': 1.0e-10}
        nu = 0.0
        def obj(nu, x):
            u = self._eval_coef(x, nu)
            Z = ((u*u).sum() - self._np)**2
            return Z

        min_v = 0.0
        for i in n:
            if i>= 0.0:
               min_v += min(i*i, 1.0)

        if min_v < self._np:
           nu = 0.0
        else:
           #R = scipy.optimize.minimize(obj, nu, args=(n,), tol=1e-20, options=options)
           bounds = [[0.0, None],]
           R = scipy.optimize.minimize(obj, nu, args=(n,), method='slsqp', tol=1.0e-50, bounds=bounds, options=options)
           nu = R.x
       #print("Here? Opt: Z= ", obj(nu, n), "nu= ", nu)

        return nu

    def _eval_coef(self, b, nu):
        "Projected occupation numbers" 
        b_ = b.copy();                 
        for i in range(len(b)):
            u = b[i]/(1.0 + nu)
            if   u <= 0.0: b_[i] = 0.0
            elif u >= 1.0: b_[i] = 1.0
            else: b_[i] = u
        return b_

    def __construct_p_from_g(self, g, dim):
        p_ = numpy.zeros(dim)
        sg = numpy.sin(g); sg = abs(sg)
        cg = numpy.cos(g); cg = abs(cg)
        for i in range(self._np):
            p_[           i] = sg[i]
            p_[self._np + i] = cg[i]
        return p_


    def _perfect_pair_projection(self, n, c): #TODO
        po = n
        Po = c @ numpy.diag(po) @ c.T

        def obj(g, P, dim):
            # construct perfect-pair occupancies
           # p_ = numpy.zeros(dim)
           # sg = numpy.sin(g); sg = abs(sg)
           # cg = numpy.cos(g); cg = abs(cg)
           # for i in range(self._np):
           #     p_[           i] = sg[i]
           #     p_[self._np + i] = cg[i]
            p_ = self.__construct_p_from_g(g, dim)
           
            # calculate orbitals 
            Q = psi4.core.Vector.from_array( numpy.einsum("i,jk->ijk", p_, P).ravel(), "")
            X = oepdev.calculate_unitary_uo_2(Q, dim).to_array(dense=True)

            # objective function value
            r = numpy.einsum("ai,bi,ab->i", X, X, P)
            Z = - r @ p_

            return Z

        # solve optimization problem for p_
        g_0 = numpy.zeros(self._np)
        bounds = None
        options = {'disp': False, 'maxiter':2000, 'ftol': 1.0e-10}
        R = scipy.optimize.minimize(obj, g_0, args=(Po, po.size), 
                                         method='slsqp', tol=1.0e-10, bounds=bounds, options=options)
        g = R.x
        p_new = self.__construct_p_from_g(g, po.size)

        # compute orbitals from p_
        Q = psi4.core.Vector.from_array( numpy.einsum("i,jk->ijk", p_new, Po).ravel(), "")
        c_new = oepdev.calculate_unitary_uo_2(Q, po.size).to_array(dense=True)

        # sort according to perfect-pair-dencending order
        idx = numpy.argsort(p_new[:self._np])[::-1]  # -> descending
        p_occ = p_new[  idx]
        p_vir = p_new[  idx[::-1] + self._np ]
        p_uno = p_new[  (2*self._np):]
        p_new = numpy.concatenate((p_occ, p_vir, p_uno))
        c_occ = c_new[:,idx]
        c_vir = c_new[:,idx[::-1] + self._np ]
        c_uno = c_new[:,(2*self._np):]
        c_new = numpy.concatenate((c_occ, c_vir, c_uno), axis=1)
        return p_new, c_new




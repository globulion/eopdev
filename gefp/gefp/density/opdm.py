#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 One-Particle Density Matrix Module.
 Bartosz Błasiak, Gundelfingen, May 2019
"""

import numpy
import psi4
import scipy.optimize
from abc import ABC, abstractmethod

__all__ = ["Density"]

class Density:
    """\
 Handles the Electron Density Distribution.

 Class methods:
  Density.natural_orbitals
  Density.generalized_density
  Density.orthogonalize_OPDM
  Density.orthogonalizer
  Density.deorthogonalizer

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

    @classmethod
    def natural_orbitals(cls, D, 
                              S                      =  None         , 
                              C                      =  None         , 
                              orthogonalize_mo       =  True         , 
                              order                  = 'descending'  , 
                              return_ao_orthogonal   =  False        , 
                              renormalize            =  False        , 
                              no_cutoff              =  False        ,
                              ignore_large_n         =  False        ):
        "Compute the Natural Orbitals from a given ODPM"
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
           if ((n.max() - 1.0) > 0.00001 or (n.min() < -0.00001)):
              raise ValueError("Unphysical NO populations detected! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))
        n[numpy.where(n<0.0)] = 0.0

        # NO cutoff
        if no_cutoff is False: no_cutoff = self.no_cutoff
        if no_cutoff != 0.0:
           ids = numpy.where(n>=self.no_cutoff)
           n = n[ids]
           U =(U.T[ids]).T

        # Order according to occupations
        if order=='ascending': 
           pass
        elif order=='descending':
           n = n[  ::-1]
           U = U[:,::-1]
        else: raise ValueError("Incorrect order of NO orbitals. Possible only ascending or descending.")

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

    #@classmethod
    def compute_1el_energy(self, D, Hcore):
        "Compute generalized 1-electron energy"
        energy = numpy.dot(D, Hcore).trace()
        return energy

    #@classmethod
    def compute_2el_energy(self, D_left, D_right, type='j'):
        "Compute generalized 2-electron energy"
        JorK = self.generalized_JK(D_left, type)
        energy = numpy.dot(JorK, D_right).trace()
        return energy

    @classmethod
    def generalized_density(cls, n, c, g=1.0):
        "Generalized matrix G = C n**g C.T"
        ng = n.copy() if g==1.0 else n**g
        G = numpy.linalg.multi_dot([c, numpy.diag(ng), c.T])
        return G

    #@classmethod
    def generalized_JK(self, D, type='j'):
        self._global_jk.C_clear()                                           
        self._global_jk.C_left_add(psi4.core.Matrix.from_array(D, ""))
        I = numpy.identity(D.shape[0], numpy.float64)
        self._global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._global_jk.compute()
        if   type.lower() == 'j': JorK = self._global_jk.J()[0].to_array(dense=True)
        elif type.lower() == 'k': JorK = self._global_jk.K()[0].to_array(dense=True)
        else: raise ValueError("Incorrect type of JK matrix. Only J or K allowed.")
        return JorK


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

    def compute(self, n, c):
        return self._density_matrix_projection(n, c)


    # --> Protected Interface <-- #

    def _density_matrix_projection(self, n, c):
        "Find n_new and C_new such that new density matrix is N-representable"                              
        if self._S is None: S = numpy.identity(len(n))
        else              : S = self._S.copy()
        # compute pre-density matrix                                                                       
        preD = Density.generalized_density(n, c) # cannot be here self.D because it is pre-density matrix!
        #A = numpy.linalg.multi_dot([S, preD, S])
        A = Density.orthogonalize_OPDM(preD, S)
        #a, b = scipy.linalg.eig(A, S)
        a, b = numpy.linalg.eigh(A)
        #a = a.real; b = b.real
        #print(" Init sum = %14.6f" % (a**2).sum()) 
                                                                                                            
        muORnu = self._find_coef(a)
                                                                                                           
        # compute the projected density matrix
        n_new = self._eval_coef(a, muORnu)
        #print((n_new**2).sum())
        C_new = b
                                                                                                            
        # sort (descending order)
        idx = numpy.argsort(n_new)[::-1]
        n_new = n_new [  idx]
        C_new = C_new [:,idx]
        C_new = numpy.dot(Density.orthogonalizer(S), C_new)
        return n_new.real, C_new.real

    @abstractmethod
    def _find_coef(self, n):
        pass

    @abstractmethod
    def _eval_coef(self, x, coeff):
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
        options = {'disp': False, 'maxiter':2000, 'ftol': 1.0e-20}
        nu = 0.0
        def obj(nu, x):
            u = self._eval_coef(x, nu)
            Z = ((u*u).sum() - self._np)**2
            return Z
        #R = scipy.optimize.minimize(obj, nu, args=(n,), tol=1e-20, options=options)
        R = scipy.optimize.minimize(obj, nu, args=(n,), method='slsqp', tol=1.0e-50, options=options)
        nu = R.x
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

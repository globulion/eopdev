#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT Solver module.
 Bartosz Błasiak, Gundelfingen, Mar 2019
"""

import os
import sys
import math
import numpy
import scipy.linalg
import scipy.optimize
import psi4
import oepdev
from abc import ABC, abstractmethod
from .functional import XCFunctional
from .partitioning import Density
from .parameters import Guess, rearrange_eigenpairs

__all__ = ["DMFT"]


def aaa(a, mu):
    "Projected occupation numbers"
    a_ = a.copy();
    for i in range(len(a)):
        u = a[i] + mu
        if   u <= 0.0: a_[i] = 0.0
        elif u >= 1.0: a_[i] = 1.0
        else: a_[i] = u
    return a_

def bbb(b, nu):
    "Projected occupation numbers"
    b_ = b.copy();
    for i in range(len(b)):
        u = b[i]/(1.0 + nu)
        if   u <= 0.0: b_[i] = 0.0
        elif u >= 1.0: b_[i] = 1.0
        else: b_[i] = u
    return b_

def find_mu(n, np):
    "Search for mu"
    options = {'disp': False, 'maxiter':250}
    mu = 0.0
    def obj(mu, x):
        u = aaa(x, mu)
        Z = (u.sum() - np)**2
        return Z
    R = scipy.optimize.minimize(obj, mu, args=(n,), tol=1e-20, options=options)
    mu = R.x
    return mu

def find_nu(n, np):
    "Search for mu"
    options = {'disp': False, 'maxiter':2050, 'ftol':1.0e-20}
    nu = 0.0
    def obj(nu, x):
        u = bbb(x, nu)
        Z = ((u*u).sum() - np)**2
        return Z
    R = scipy.optimize.minimize(obj, nu, args=(n,), method='slsqp', tol=1.0e-50, options=options)
    nu = R.x
    return nu

def density_matrix_projection(n, c, S, np, type='d'):
    "Find n_new and C_new such that new density matrix is N-representable"
    if type.lower() == 'd':
       func_find = find_mu
       func_coef = aaa
    else:
       func_find = find_nu
       func_coef = bbb
 
    # compute pre-density matrix                                                                       
    preD = Density.generalized_density(n, c) # cannot be here self.D because it is pre-density matrix!
    #A = numpy.linalg.multi_dot([S, preD, S])
    A = Density.orthogonalize_OPDM(preD, S)
    #a, b = scipy.linalg.eig(A, S)
    a, b = numpy.linalg.eigh(A)
    #a = a.real; b = b.real
    #print(" Init sum = %14.6f" % (a**2).sum()) 

    muORnu = func_find(a, np)
                                                                                                       
    # compute the projected density matrix
    n_new = func_coef(a, muORnu)
    #print((n_new**2).sum())
    C_new = b

    # sort (descending order)
    idx = numpy.argsort(n_new)[::-1]
    n_new = n_new [  idx]
    C_new = C_new [:,idx]
    C_new = numpy.dot(Density.orthogonalizer(S), C_new)
    return n_new.real, C_new.real



class ElectronCorrelation:
    """\
 The Electron Correlation: Dynamic and Non-dynamic Correlation.
"""
    @staticmethod
    def degree_of_correlation(dmft, scalar=True):
        "Compute scalar or matrix degree of correlation. Return I_d, I_n"
        n  = dmft.N
        ns = n.copy(); ns[ns<0.0] = 0.0
        #S   = numpy.sqrt(ns).sum()
        #N   = ns.sum()
        if scalar:
           I_n = (ns*(1.0 - ns)).sum()                            
           I_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - I_n
        else:
           c   = dmft.C
           I_n = (ns*(1.0 - ns))
           I_d = numpy.sqrt(abs(ns*(1.0 - ns))) / 2.0 - I_n
           I_n = numpy.linalg.multi_dot([c, numpy.diag(I_n), c.T])
           I_d = numpy.linalg.multi_dot([c, numpy.diag(I_d), c.T])
        return I_d, I_n



class DMFT(ABC,ElectronCorrelation):
    """\
 The Density Matrix Functional Theory.
"""
    # Defaults
    default_xc_functional = XCFunctional.create('hf') # Exchange-Correlation Functional
    default_algorithm     = 'proj-d'                  # Projected Gradient Algorithm Abbreviation            
    default_convergence   =  0.00001                  # Convergence in total energy [A.U.]
    default_maxiter       =  100                      # Maximum number of iterations
    default_verbose_run   =  True                     # Show additional information?
    default_g0            =  0.0001                   # Steepest-descents step scale in the first iteration
    default_v_ext         =  None                     # External 1-electron potential
    default_guess         = 'hcore'                   # Guess for ODPM
    default_step_mode     = 'search'                  # Steepest-Descents step (simple line search)

    def __init__(self, wfn, xc_functional, v_ext=default_v_ext, guess=default_guess, step_mode=default_step_mode):
        "Initialize"
        # Protected variable namespace
        self._current_energy           = None         # Total Energy                                        
        self._current_orbitals         = None         # Natural Orbitals
        self._current_occupancies      = None         # Natural Orbital Occupation Numbers
        self._current_density          = None         # 1-Particle Density Matrix in AO Basis
                                                                                                            
        self._xc_functional            = None         # XC Functional Object
        self._mol                      = None         # Molecule Object
        self._wfn                      = None         # Wavefunction Object
        self._bfs                      = None         # Basis Set Object
        self._bfs_name                 = None         # Basis Set Abbreviation
        self._mints                    = None         # Integral Calculator Object
        self._ints                     = None         # Integral Transform Object (for MO-SCF basis)
        self._jk                       = None         # JK Object (for AO basis)
        self._e_nuc                    = None         # Nuclear Repulsion Energy
        self._np                       = None         # Number of Electron Pairs
        self._Ca                       = None         # AO-MO Matrix at Hartree-Fock Level (LCAO-MO)
        self._S                        = None         # AO Overlap Matrix
        self._X                        = None         # AO Orthogonalizer 
        self._Y                        = None         # AO Deorthogonalizer
        self._H                        = None         # AO Core Hamiltonian Matrix
        self._V_ext                    = None         # AO External Potential Matrix
        self._H_mo                     = None         # MO-SCF Core Hamiltonian + External Potential Matrix
        self._step_mode                = None         # Mode of Steepest-Descents Minimization Steps
        self._mode_xc_gradient         = None         # Mode of E_XC gradient

        # Initialize all variables
        self.__common_init(wfn, xc_functional, v_ext, guess, step_mode)

        # Sanity checks
        if self._xc_functional.abbr.lower() != 'hf' and self.abbr.lower() == 'dmft-projd':
           print("""\
 Warning! The D-set is not Lipshitz with %s functional. 
 This will probably result in lack of convergence. Use P-set instead.""" % xc_functional.abbr.upper())

        # Abstract base meta class
        super(DMFT, self).__init__()


    # ---- Public Interface ----- #

    @classmethod
    def create(cls, wfn, xc_functional = default_xc_functional, 
                         v_ext         = default_v_ext        ,
                         guess         = default_guess        ,
                         algorithm     = default_algorithm    ,
                         step_mode     = default_step_mode    ,        **kwargs):
        """\
 Create DMFT solver. 
"""
        args = [wfn, xc_functional, v_ext, guess, step_mode]

        if   algorithm.lower() == 'proj-d': solver = DMFT_ProjD(*args, **kwargs)
        elif algorithm.lower() == 'proj-p': solver = DMFT_ProjP(*args, **kwargs)
        elif algorithm.lower() == 'nc'    : solver = DMFT_NC   (*args, **kwargs)
        elif algorithm.lower() == 'pc'    : solver = DMFT_PC   (*args, **kwargs)
        else: raise ValueError("Chosen algorithm is not available! Mistyped?")
        return solver

    def run(self, conv    = default_convergence, 
                  maxit   = default_maxiter    , 
                  verbose = default_verbose_run, 
                  g_0     = default_g0         , **kwargs):
        """\
 Run the DMFT calculations.
"""

        # Run!
        success = self._run_dmft(conv, maxit, verbose, g_0, **kwargs)

        # Crash?
        if not success: 
            print(" DMFT Iterations failed!")
            sys.exit(1)

        return success


    # ---- Public Interface: Properties ----- #

    def set_gradient_mode(self, exact=False, approx=False, num=False):
        "Set the mode of computation of derivatives of XC energy wrt P matrix (not D!)"
        assert any((exact, approx, num)) is True, "You must set one of the following to True: exact, approx, num"
        assert all((exact, approx)) is False
        assert all((exact, num)) is False
        assert all((num, approx)) is False
        if exact  is True: self._mode_xc_gradient = 'exact'
        if approx is True: self._mode_xc_gradient = 'approximate'
        if num    is True: self._mode_xc_gradient = 'numerical'


    @property
    def E(self): 
        "Total Energy"
        return self._current_energy
    @property
    def D(self): 
        "1-Particle Density Matrix"
        return self._current_density.matrix()
    @property
    def C(self): 
        "NO Orbitals"
        return self._current_orbitals
    @property
    def N(self): 
        "NO Occupations"
        return self._current_occupancies


    @property
    def scalar_correlation(self):
        "Scalar Degree of Electron Correlation"
        return self.degree_of_correlation(self, scalar=True)
    @property
    def matrix_correlation(self):
        "Matrix Degree of Electron Correlation"
        return self.degree_of_correlation(self, scalar=False)


    @property
    def abbr(self): 
        "Abbreviation symbol for the DMFT optimization algorithm"
        pass

    @staticmethod
    @abstractmethod
    def name(): 
        "Full name for the method implemented"
        pass


    # --- Protected Interface --- #

    def _run_dmft(self, conv, maxit, verbose, g_0, **kwargs):
        "DMFT Iterations"

        restart = True if 'restart' in kwargs.keys() else False

        if verbose and not restart: 
           print(" Running %s:%s/%s" % (self.abbr, self._xc_functional.abbr, self._bfs_name))
        if verbose and restart:
           print(" Restarting %s:%s/%s" % (self.abbr, self._xc_functional.abbr, self._bfs_name))

        if not restart:

           # [0] Initialize                                                                  
           self._iteration = 0                                                             
           success = False
                                                                                             
           # [1] Compute Guess
           x_0 = self._guess()
                                                                                             
           # [2] Starting energy
           self._current_energy = self._minimizer(x_0)
           self._E_old = self._current_energy
           if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (self._iteration, self._E_old))
                                                                                             
           # [3] First iteration
           self._iteration += 1
           self._x_old_2    = x_0
           self._x_old_1    = self._step_0(self._x_old_2, g_0)
           self._x_old_1    = self._density(self._x_old_1)
                                                                                             
           self._current_energy = self._minimizer(self._x_old_1)
           self._E_new = self._current_energy
           if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (self._iteration, self._E_new))

        # [4] Further iterations
        self._iteration += 1
        stop       = False
        while stop is False:

            # [4.1] New guess
            x_new = self._step(self._x_old_1, self._x_old_2)
            x_new = self._density(x_new)

            # [4.2] Current energy
            self._current_energy = self._minimizer(x_new)
            self._E_new = self._current_energy
            if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (self._iteration, self._E_new))
                                                                                  
            # [4.3] Converged?
            if abs(self._E_new-self._E_old) < conv: 
               stop    = True
               success = True
                                                                                  
            # [4.4] Maximum iterations exceeded?
            if self._iteration >= maxit: 
               stop    = True
               success = False
                                                                                  
            # [4.5] Prepare for next iteration
            self._iteration += 1
            self._E_old      = self._E_new
            self._x_old_2    = self._x_old_1.copy()
            self._x_old_1    = x_new  .copy()
        
        # [5] Finish
        if verbose and success:
           print(" DMFT iterations converged.")
           print(" Final Energy = %14.8f" % self._current_energy)

        if verbose and not success:
           print(" DMFT iterations did not converge.")

        return success

    def _compute_no_exchange_energy(self):
        "Compute Total Energy without Exchange-Correlation Energy"
        E_N = self._e_nuc
        E_1 = self._compute_hcore_energy()
        E_H = self._compute_hartree_energy()
        E   = E_N + E_1 + E_H
        return E

    def _compute_no_exchange_gradient_D(self, x):
        "1-electron and Hartree part of gradient wrt D. Returned in SCF-MO basis."
        D = self._current_density.matrix() # Must be in MO-SCF basis!
        #D = x.matrix()
        gradient_H  = self._H_mo
        gradient_J  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        gradient  = Guess.create(matrix = gradient_H + 2.0 * gradient_J)
        return gradient

    def _compute_no_exchange_gradient_P(self, x):
        P  = x.matrix()
        P2 = numpy.dot(P, P)
        J  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P2, ""))[0].to_array(dense=True)
        v  = self._H_mo + 2.0 * J
        gradient = numpy.dot(v, P)
        gradient+= gradient.T
        gradient  = Guess.create(matrix = 2.0 * gradient)
        return gradient

    def _compute_no_exchange_gradient_nc(self, x):
        "Compute Gradient excluding the Exchange-Correlation Part"
        grad_n = self.__grad_n_no_exchange(*x.unpack())
        grad_c = self.__grad_c_no_exchange(*x.unpack())
        grad = Guess.create(grad_n, grad_c, None, 'nc')
        return grad

    def _compute_no_exchange_gradient_pc(self, x):#OK
        "Compute Gradient excluding the Exchange-Correlation Part"
        grad_p, grad_c = self.__grad_pc_no_exchange(*x.unpack())#OK
        grad = Guess.create(grad_p, grad_c, None, 'nc')
        return grad


    def _estimate_step_size(self, dx, dg):
        "Estimate step length of steepest descents"
        DX = dx.matrix().ravel()
        DG = dg.matrix().ravel()
        norm = numpy.dot(DG, DG)
        g    = numpy.dot(DX, DG)/ norm
        g    = abs(g)
        return g

    def _correct_negative_occupancies(self, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns

                  
    @abstractmethod
    def _compute_hcore_energy(self):
        "1-Electron Energy: HCore + Vext"
        pass
           
    @abstractmethod
    def _compute_hartree_energy(self):
        "2-Electron Energy: Hartree"
        pass

    @abstractmethod
    def _compute_initial_NOs(self):
        "Initial natural orbital analysis"
        pass


    @abstractmethod
    def _minimizer(self, x):
        "Minimizer function"
        pass

    @abstractmethod
    def _guess(self):
        "Initial guess"
        pass

    @abstractmethod
    def _density(self, x):
        "1-particle density matrix in AO basis"
        pass

    @abstractmethod
    def _gradient(self, x):
        "Gradient"
        pass

    @abstractmethod
    def _step_0(self, x0, g0):
        "Steepest-descents initial step"
        pass

    @abstractmethod
    def _step(self, x1, x2):
        "Steepest-descents further step"
        pass



    # --- Private Interface --- #

    def __common_init(self, wfn, xc_functional, V_ext, guess, step_mode):
        "Initialize"

        # ---->  Basic Objects <---- #

        # Wavefunction
        self._wfn = wfn
        # Basis set name
        self._bfs_name = wfn.basisset().name()
        # XC functional
        self._xc_functional = xc_functional
        self._mode_xc_gradient = 'exact'

        # Molecule
        self._mol = self._wfn.molecule()
        # Basis set
        self._bfs = self._wfn.basisset()
        # Integral calculator
        self._mints = psi4.core.MintsHelper(self._bfs)
        #print("Computing AO ERI")
        #self._ao_eri= self._mints.ao_eri()
        #print("Done!")
        # JK object
        self._jk = psi4.core.JK.build(self._bfs, jk_type="Direct")
        self._jk.set_memory(int(5e8))
        self._jk.initialize()
        self._xc_functional.set_jk(self._jk)
        self._xc_functional.set_wfn(self._wfn)
        #self._xc_functional._ao_eri = self._ao_eri
        # Nuclear repulsion energy
        self._e_nuc = self._mol.nuclear_repulsion_energy()
        # Number of electron pairs
        self._np = self._wfn.nalpha()
        # Step mode
        self._step_mode = step_mode
                                                                                               
        # ---->  Constant AO matrices <---- #

        # SCF LCAO-MO coefficients
        self._Ca= self._wfn.Ca_subset("AO","ALL").to_array(dense=True)
        # Overlap integrals and orthogonalizer
        self._S = self._wfn.S().to_array(dense=True)
        self._X = Density.  orthogonalizer(self._S)
        self._Y = Density.deorthogonalizer(self._S)
        # Hcore matrix
        self._H = self._wfn.H().to_array(dense=True)
        # External potential
        self._V_ext = V_ext
        if V_ext is not None and guess == 'current':
           raise ValueError(" External potential can only be set for 'HCore' guess!")

        # ---->  MO integrals in MO-SCF basis <---- #

        # Integral Transform
        psi4.check_iwl_file_from_scf_type(psi4.core.get_global_option('SCF_TYPE'), self._wfn)
        mo_all = psi4.core.MOSpace.all()
        spaces = [mo_all]
        trans_type = psi4.core.IntegralTransform.TransformationType.Restricted
        self._ints = psi4.core.IntegralTransform(self._wfn, spaces, trans_type)
        self._ints.transform_tei(mo_all, mo_all, mo_all, mo_all)
        psi4.core.print_out('Integral transformation complete!\n')
        self._xc_functional.set_ints(self._ints)

        # ---->  Current OPDM and total energy <---- #

        # Guess based on current density matrix from the input wavefunction
        if guess == 'current':
           assert self._V_ext is None, "External potential cannot be set with guess=current!"
           E = self._wfn.energy()
           D = self._wfn.Da().to_array(dense=True)
        # Guess based on the one-electron Hamiltonian
        elif guess == 'hcore':
           if self._V_ext is not None: self._H += self._V_ext
           e, c = numpy.linalg.eigh(numpy.linalg.multi_dot([self._X, self._H, self._X]))
           c = c[:,::-1]
           c = numpy.dot(self._X, c)
           E = self._e_nuc + 2.0 * e.sum()
           D = numpy.zeros((self._bfs.nbf(), self._bfs.nbf()), numpy.float64)
           for i in range(self._np):
               D += numpy.outer(c[:,i], c[:,i])
        else:
           raise ValueError("Only 'Current' or 'HCore' ODPM's are supported as starting points.")
        self._current_energy  = E
        self._current_density = Density(D, self._jk)
        # Initial Natural Orbitals
        self._current_occupancies, self._current_orbitals = self._compute_initial_NOs() # it also changes _current_density
        # H_core + V_ext in MO-SCF basis
        self._H_mo = numpy.linalg.multi_dot([self._Ca.T, self._H, self._Ca])
        return

    def __grad_n_no_exchange(self, n, c):
        "Energy gradient wrt NO occupation numbers"
        nn = len(n)
                                                                     
        # 1-electron contributon
        Hmm = numpy.linalg.multi_dot([c.T, self._H, c]).diagonal()
        grad = 2.0 * Hmm
                                                                     
        # 2-electron contribution
        I = numpy.identity(nn, numpy.float64)
                                                                     
        # J-type
        self._jk.C_clear()
        self._jk.C_left_add(psi4.core.Matrix.from_array(self._current_density.matrix(), ""))
        self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._jk.compute()
        J = self._jk.J()[0].to_array(dense=True)
        grad += 4.0 * numpy.linalg.multi_dot([c.T, J, c]).diagonal()
                                                                     
        # K-type
        # ---> moved to XC_Functional
        return grad

    def __grad_pc_no_exchange(self, p, c):#OK
        "Energy gradient wrt NO occupation numbers"
             
        # ===> Gradient wrt p <=== #

        # 1-electron contributon
        Ham = numpy.dot(self._H_mo, c)
        Hmm = numpy.dot(c.T, Ham).diagonal()
        grad_p = 4.0 * Hmm * p

                                                                     
        # 2-electron contribution: J-type
        P2 = Density.generalized_density(p, c, 2.0)
        P2_psi = psi4.core.Matrix.from_array(P2, "")
        J = oepdev.calculate_JK_r(self._wfn, self._ints, P2_psi)[0].to_array(dense=True)
        Jam = numpy.dot(J, c)
        Jmm = numpy.dot(c.T, Jam).diagonal()
        grad_p += 8.0 * Jmm * p

        # K-type
        # ---> in XC_Functional

        # ===> Gradient wrt c <=== #

        pp = p*p
        # 1-electron contribution
        grad_c = 4.0 * Ham * pp
        # 2-electron contribution
        grad_c+= 8.0 * Jam * pp
        #grad_c = grad_c * pp

        return grad_p, grad_c


    def __grad_c_no_exchange(self, n, c):
        "Energy gradient wrt LCAO-NO wavefunction coefficients"
        nn = len(n)
                                                                                             
        # 1-electron contributon
        Ham = numpy.dot(self._H, c)
        grad = Ham * n

        # 2-electron contribution
        I = numpy.identity(nn, numpy.float64)

        # J-type
        self._jk.C_clear()                                           
                                                                                             
        self._jk.C_left_add(psi4.core.Matrix.from_array(self._current_density.matrix(), ""))
        self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._jk.compute()
        J = self._jk.J()[0].to_array(dense=True)
        J = 2.0 * numpy.dot(J, c)
        for m in range(nn): grad[:,m] += n[m] * J[:,m]

        # K-type
        # ---> moved to XC_Functional
                                                                                             
        # finalize
        grad *= 4.0
                                                                                             
        # transform gradient to oAO basis
        grad = numpy.dot(self._X, grad)
        return grad


class DMFT_AO(DMFT):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_AO, self).__init__(wfn, xc_functional, v_ext, guess, step)

    # --- Implementation (Protected Interface) --- #

    def _compute_hcore_energy(self):
        "1-Electron Energy"
        H = self._H.copy()
        D = self._current_density.matrix()
        E = 2.0 * self._current_density.compute_1el_energy(D, H)
        return E

    def _compute_hartree_energy(self):
        "2-Electron Energy: Hartree ---> AO basis"
        D = self._current_density.matrix()
        E = 2.0 * self._current_density.compute_2el_energy(D, D, type='j')
        return E

    def _compute_initial_NOs(self):
        "C: AO to NO matrix"
        n, c = Density.natural_orbitals(self._current_density.matrix(), self._S, self._Ca, orthogonalize_mo=True,
                                       order='descending', no_cutoff=0.0, renormalize=False, return_ao_orthogonal=False,
                                       ignore_large_n=False)
        return n, c

    def _step(self, x1, x2):
        "Steepest-descents step. Back-transforms to OAO basis for simplicity"
        n1, c1 = x1.unpack() 
        n2, c2 = x2.unpack()

        # transform C to OAO basis
        c1_ = numpy.dot(self._Y, c1)
        c2_ = numpy.dot(self._Y, c2)

        x_old_1 = Guess.create(n1, c1_, None, 'nc')
        x_old_2 = Guess.create(n2, c2_, None, 'nc')
                                                                                       
        gradient_1 = self._gradient(Guess.create(n1, c1, None, 'nc'))
        gradient_2 = self._gradient(Guess.create(n2, c2, None, 'nc'))
                                                                                       
        g = 0.1
        if self._step_mode.lower() != 'constant': 
           g = self._estimate_step_size(x_old_1 - x_old_2, gradient_1 - gradient_2)
        x_new = x_old_1 - g * gradient_1

        # transform back to AO basis
        n, c_ = x_new.unpack()
        x_new = Guess.create(n, numpy.dot(self._X, c_), None, 'nc')
        return x_new

    def _step_0(self, x0, g0):
        "Steepest-descents step. Back-transforms to OAO basis for simplicity"
        n2, c2 = x0.unpack()

        # transform C to OAO basis
        c2_ = numpy.dot(self._Y, c2)
        x_old_2 = Guess.create(n2, c2_, None, 'nc')

        # compute new guess
        gradient_2 = self._gradient(Guess.create(n2, c2, None, 'nc'))
        x_new = x_old_2 - g0 * gradient_2
        
        # transform back to AO basis
        n, c_ = x_new.unpack()
        x_new = Guess.create(n, numpy.dot(self._X, c_), None, 'nc')
        return x_new





class DMFT_MO(DMFT):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_MO, self).__init__(wfn, xc_functional, v_ext, guess, step)

    # --- Implementation (Protected Interface) --- #

    def _compute_hcore_energy(self):
        "1-Electron Energy"
        H = self._H_mo
        D = self._current_density.matrix()
        E = 2.0 * numpy.dot(H, D).trace()
        return E

    def _compute_hartree_energy(self):
        "2-Electron Energy: Hartree ---> MO-SCF basis"
        D = self._current_density.matrix()
        J = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        E = 2.0 * numpy.dot(J, D).trace()
        return E

    def _compute_initial_NOs(self):
        "C: MO-SCF to NO matrix"
        D = self._current_density.matrix()
        D_mo = numpy.linalg.multi_dot([self._Ca.T, self._S, D, self._S, self._Ca])
        self._current_density.set_D(D_mo)
        n, c = Density.natural_orbitals(D_mo, None, None, orthogonalize_mo=False,
                                       order='descending', no_cutoff=0.0, renormalize=False, return_ao_orthogonal=False,
                                       ignore_large_n=True) # change to False!
        return n, c

    def _step(self, x1, x2):
        "Steepest-descents step."
        gradient_1 = self._gradient(x1)
        gradient_2 = self._gradient(x2)

        g = 0.0000000001
        if self._step_mode.lower() != 'constant': 
           g = self._estimate_step_size(x1 - x2, gradient_1 - gradient_2)
        x_new = x1 - g * gradient_1
        x_new.update()
        return x_new

    def _step_0(self, x0, g0):
        "Steepest-descents step."
        gradient_2 = self._gradient(x0)
        x_new = x0 - g0 * gradient_2
        x_new.update()
        return x_new




class DMFT_NC(DMFT_AO):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_NC, self).__init__(wfn, xc_functional, v_ext, guess, step)

    # --- Implementation (Public) --- #

    @staticmethod
    def name(): return "DMFT with Gradient Projection on D Set: Optimization in N and C set."

    @property
    def abbr(self): return "DMFT-NC"


    # --- Implementation (Protected Interface) --- #

    def _minimizer(self, x):
        "Minimizer function: Total Energy"
        E_H  = self._compute_no_exchange_energy()
        E_XC = self._xc_functional.energy_D(x, mode='ao')
        E    = E_H + E_XC
        return E

    def _guess(self):
        "Initial guess"
        x = Guess.create(n=self._current_occupancies, c=self._current_orbitals, t='nc') 
        return x

    def _density(self, x):
        "1-particle density matrix in AO basis"
        n, c = x.unpack()
        n, c = density_matrix_projection(n, c, self._S, self._np, type='d')
        self._current_occupancies = n
        self._current_orbitals    = c
        D = self._current_density.generalized_density(n, c)
        self._current_density.set_D(D)
        new_guess = Guess.create(n, c, None, 'nc')
        return new_guess

    def _gradient(self, x):
        "Gradient"
        dE_n   , dE_c    = self._compute_no_exchange_gradient_nc(x).unpack()
        dE_n_xc, dE_c_xc = self._xc_functional      .gradient_nc(x).unpack()
        dE_c_xc = numpy.dot(self._X, dE_c_xc) # oAO basis
        dE_n += dE_n_xc
        dE_c += dE_c_xc
        gradient = Guess.create(n=dE_n, c=dE_c, t='nc')
        return gradient


class DMFT_PC(DMFT_MO):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_PC, self).__init__(wfn, xc_functional, v_ext, guess, step)
        self.g = 0.1


    @staticmethod
    def name(): return "DMFT with Gradient Projection on PC Set: MO-SCF basis"

    @property
    def abbr(self): return "DMFT-PC"


    # --- Implementation (Protected Interface) --- #

    def _minimizer(self, x):#OK
        "Minimizer function: Total Energy"
        E_H  = self._compute_no_exchange_energy()
        E_XC = self._xc_functional.energy_pc(x)
        E    = E_H + E_XC
        return E

    def _guess(self):#OK
        "Initial guess"
        p = numpy.sqrt(self._correct_negative_occupancies(self._current_occupancies))
        x = Guess.create(n=p, c=self._current_orbitals, t='nc') 
        return x

    def _density(self, x):#OK
        "1-particle density matrix in MO basis: P-projection of PC set"
        p, c = x.unpack()
        p, c = density_matrix_projection(p, c, numpy.identity(len(p)), self._np, type='p')
        self._current_occupancies = p*p
        self._current_orbitals    = c
        D = Density.generalized_density(p, c, 2.0)
        self._current_density.set_D(D)
        new_guess = Guess.create(n=p, c=self._current_orbitals, t='nc')
        return new_guess

    def _gradient(self, x):#OK
        "Gradient"
        gradient = self._compute_no_exchange_gradient_pc(x) #OK
        gradient+= self._xc_functional.gradient_pc(x)
        return gradient

    def _step(self, x1, x2):#OK
        "Steepest-descents step."

        #x2 = self.__rearrange_eigenpairs(x2, x1)
        #x1 = self.__rearrange_eigenpairs(x1, x2)

        gradient_1 = self._gradient(x1)
        gradient_2 = self._gradient(x2)

        self.g *= 0.99
        g = self.g
        if self._step_mode.lower() != 'constant': 
           g = self._estimate_step_size(x1 - x2, gradient_1 - gradient_2)
        x_new = x1 - g * gradient_1
        x_new.update()
        return x_new

    def _step_0(self, x0, g0):#OK
        "Steepest-descents step."
        gradient_2 = self._gradient(x0)
        x_new = x0 - g0 * gradient_2
        x_new.update()
        return x_new

    
    # ---> Private interface <---- #

    def __rearrange_eigenpairs(self, x2, x1):
        p1, c1 = x1.unpack()
        p2, c2 = x2.unpack()
        p2, c2 = rearrange_eigenpairs(p2, c2, c1)
        x2_new = Guess.create(p2, c2, None, 'nc')
        return x2_new






class DMFT_ProjD(DMFT_MO):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_ProjD, self).__init__(wfn, xc_functional, v_ext, guess, step)


    @staticmethod
    def name(): return "DMFT with Gradient Projection on D Set"

    @property
    def abbr(self): return "DMFT-ProjD"


    # --- Implementation (Protected Interface) --- #

    def _minimizer(self, x):
        "Minimizer function: Total Energy"
        E_H  = self._compute_no_exchange_energy()
        E_XC = self._xc_functional.energy_D(x, mode='scf-mo')
        E    = E_H + E_XC
        return E

    def _guess(self):
        "Initial guess"
        x = Guess.create(n=self._current_occupancies, c=self._current_orbitals, t='matrix') 
        return x

    def _density(self, x):
        "1-particle density matrix in MO basis: D-projection"
        n, c = x.unpack()
        n, c = density_matrix_projection(n, c, numpy.identity(len(n)), self._np, type='d')
        self._current_occupancies = n
        self._current_orbitals    = c
        D = self._current_density.generalized_density(n, c)
        self._current_density.set_D(D)
        new_guess = Guess.create(n=self._current_occupancies, c=self._current_orbitals, t='matrix')
        return new_guess

    def _gradient(self, x):
        "Gradient"
        gradient = self._compute_no_exchange_gradient_D(x)
        gradient+= self._xc_functional.gradient_D(x)
        return gradient

    def _step(self, x1, x2):
        "Steepest-descents step."
        gradient_1 = self._gradient(x1)
        gradient_2 = self._gradient(x2)

        g = 0.5
        if self._step_mode.lower() != 'constant': 
           g = self._estimate_step_size(x1 - x2, gradient_1 - gradient_2)
        x_new = x1 - g * gradient_1
        x_new.update()
        return x_new

    def _step_0(self, x0, g0):
        "Steepest-descents step."
        gradient_2 = self._gradient(x0)
        x_new = x0 - g0 * gradient_2
        x_new.update()
        return x_new




class DMFT_ProjP(DMFT_MO):
    def __init__(self, wfn, xc_functional, v_ext, guess, step):
        super(DMFT_ProjP, self).__init__(wfn, xc_functional, v_ext, guess, step)


    @staticmethod
    def name(): return "DMFT with Gradient Projection on P Set"

    @property
    def abbr(self): return "DMFT-ProjP"


    # --- Implementation (Protected Interface) --- #

    def _minimizer(self, x):
        "Minimizer function: Total Energy"
        E_H  = self._compute_no_exchange_energy()
        E_XC = self._xc_functional.energy_P(x)
        E    = E_H + E_XC
        return E

    def _guess(self):
        "Initial guess"
        p = numpy.sqrt(self._correct_negative_occupancies(self._current_occupancies))
        x = Guess.create(n=p, c=self._current_orbitals, t='matrix') 
        return x

    def _density(self, x):
        "1-particle density matrix in MO basis: P-projection"
        p, c = x.unpack()
        p, c = density_matrix_projection(p, c, numpy.identity(len(p)), self._np, type='p')
        self._current_occupancies = p**2
        self._current_orbitals    = c
        D = self._current_density.generalized_density(p, c, 2.0)
        self._current_density.set_D(D)
        new_guess = Guess.create(n=p, c=self._current_orbitals, t='matrix')
        return new_guess

    def _gradient(self, x):
        "Gradient"
        gradient = self._compute_no_exchange_gradient_P(x)
        if   self._mode_xc_gradient == 'numerical'  : gradient+= self._xc_functional.gradient_P_numerical(x)
        elif self._mode_xc_gradient == 'approximate': gradient+= self._xc_functional.gradient_P_approximate(x)
        else:                                         gradient+= self._xc_functional.gradient_P(x)
        return gradient

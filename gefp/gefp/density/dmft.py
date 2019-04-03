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
from abc import ABC, abstractmethod
from .functional import XCFunctional
from .partitioning import Density

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
    mu = 0.0
    def obj(mu, x):
        u = aaa(x, mu)
        Z = (u.sum() - np)**2
        return Z
    R = scipy.optimize.minimize(obj, mu, args=(n,))
    mu = R.x
    return mu

def find_nu(n, np):
    "Search for mu"
    nu = 0.0
    def obj(nu, x):
        u = bbb(x, nu)
        Z = ((u*u).sum() - np)**2
        return Z
    R = scipy.optimize.minimize(obj, nu, args=(n,))
    nu = R.x
    return nu

def density_matrix_projection(n, c, S, np, type='d'):#OK
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
    a = a.real; b = b.real
    #print(" Init sum = %14.6f" % a.sum()) 
                                                                                                     
    muORnu = func_find(a, np)
                                                                                                       
    # compute the projected density matrix
    n_new = func_coef(a, muORnu)
    C_new = b

    # sort (descending order)
    idx = numpy.argsort(n_new)[::-1]
    n_new = n_new [  idx]
    C_new = C_new [:,idx]
    C_new = numpy.dot(Density.orthogonalizer(S), C_new)
    return n_new.real, C_new.real



class DMFT(ABC):
    """\
 The Density Matrix Functional Theory.
"""
    # Defaults
    default_xc_functional = XCFunctional.create('hf') # XC functional abbreviation
    default_algorithm     = 'proj-d'           # Projected gradient algorithm abbreviation           
    default_convergence   =  0.00001           # Convergence in total energy [A.U.]
    default_maxiter       =  100               # Maximum number of iterations
    default_verbose_run   =  True              # Show additional information?
    default_g0            =  0.0001            # Steepest-descents step scale in the first iteration
    default_v_ext         =  None              # External 1-electron potential
    default_guess         = 'hcore'            # Guess for ODPM

    def __init__(self, wfn, xc_functional, v_ext=default_v_ext, guess=default_guess):
        "Initialize"
        # Protected variable namespace
        self._current_energy           = None       # Total Energy
        self._current_orbitals         = None       # Natural Orbitals
        self._current_occupancies      = None       # Natural Orbital Occupation Numbers
        self._current_density          = None       # 1-Particle Density Matrix in AO Basis

        self._xc_functional            = None       # XC Functional Object
        self._mol                      = None       # Molecule Object
        self._wfn                      = None       # Wavefunction Object
        self._bfs                      = None       # Basis Set Object
        self._bfs_name                 = None       # Basis Set Abbreviation
        self._mints                    = None       # Integral Calculator Objects
        self._jk                       = None       # JK Object
        self._e_nuc                    = None       # Nuclear Repulsion Energy
        self._np                       = None       # Number of Electron Pairs
        self._Ca                       = None       # AO-MO Matrix at Hartree-Fock Level (LCAO-MO)
        self._S                        = None       # AO Overlap Matrix
        self._X                        = None       # AO Orthogonalizer 
        self._Y                        = None       # AO Deorthogonalizer
        self._H                        = None       # AO Core Hamiltonian Matrix
        self._V_ext                    = None       # AO External Potential Matrix

        # Initialize all variables
        self.__common_init(wfn, xc_functional, v_ext, guess)

        # Sanity checks
        if self._xc_functional.abbr.lower() != 'hf' and self.abbr.lower() == 'dmft-projd':
           print("""\
 Warning! The D-set is not Lipshitz with %s functional. 
 This will probably result in lack of convergence. Use P-set instead.""" % xc_functional.abbr.upper())

        # Call constructor for meta class
        super(DMFT, self).__init__()


    # ---- Public Interface ----- #

    @classmethod
    def create(cls, wfn, xc_functional = default_xc_functional, 
                         v_ext         = default_v_ext        ,
                         guess         = default_guess        ,
                         algorithm     = default_algorithm    ,        **kwargs):
        """\
 Create DMFT solver. 
"""
        args = [wfn, xc_functional, v_ext, guess]
        if   algorithm.lower() == 'proj-d': solver = DMFT_ProjD(*args, **kwargs)
        elif algorithm.lower() == 'proj-p': solver = DMFT_ProjP(*args, **kwargs)
        elif algorithm.lower() == 'nc'    : solver = DMFT_NC   (*args, **kwargs)
        else: raise ValueError("Chosen algorithm is not available! Mistyped?")
        return solver

    def run(self, conv    = default_convergence, 
                  maxit   = default_maxiter    , 
                  verbose = default_verbose_run, 
                  g_0     = default_g0         , **kwargs):
        """\
 Run the DMFT calculations.
"""

        if verbose: print(" Running %s:%s/%s" % (self.abbr, self._xc_functional.abbr, self._bfs_name))

        # Run!
        success = self._run_dmft(conv, maxit, verbose, g_0, **kwargs)
        return success


    # ---- Public Interface: Properties ----- #

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
        return self._current_occupations

    @property
    @abstractmethod
    def abbr(self): 
        "Abbreviation symbol for the DMFT optimization algorithm"
        pass

    @staticmethod
    @abstractmethod
    def name(): pass


    # --- Protected Interface --- #

    def _run_dmft(self, conv, maxit, verbose, g_0, **kwargs):#OK
        "DMFT Iterations"

        iteration = 0                                                             
        success = False

        # [1] Compute Guess
        x_0 = self._guess()
                                                                                  
        # [2] Starting energy
        nn=self._bfs.nbf()
        self._current_energy = self._minimizer(x_0)
        E_old = self._current_energy
        if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (iteration, E_old))
                                                                                  
        # [3] First iteration
        iteration += 1
        x_old_2    = x_0
        x_old_1    = self._step_0(x_old_2, g_0)
        x_old_1    = self._density(x_old_1)

        self._current_energy = self._minimizer(x_old_1)
        E_new = self._current_energy
        if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (iteration, E_new))
                                                                                  
        # [4] Further iterations
        iteration += 1
        stop       = False
        while stop is False:

            # [4.1] New guess
            x_new = self._step(x_old_1, x_old_2)
            x_new = self._density(x_new)

            # [4.2] Current energy
            self._current_energy = self._minimizer(x_new)
            E_new = self._current_energy
            if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (iteration, E_new))
                                                                                  
            # [4.3] Converged?
            if abs(E_new-E_old) < conv: 
               stop    = True
               success = True
                                                                                  
            # [4.4] Maximum iterations exceeded?
            if iteration >= maxit: 
               stop    = True
               success = False
                                                                                  
            # [4.5] Prepare for next iteration
            iteration += 1
            E_old      = E_new
            x_old_2    = x_old_1.copy()
            x_old_1    = x_new  .copy()
        
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

    def _compute_no_exchange_gradient_D(self):#TODO
        return NotImplementedError
    def _compute_no_exchange_gradient_P(self):#TODO
        return NotImplementedError
    def _compute_no_exchange_gradient_nc(self, n, c):
        "Compute Gradient excluding the Exchange-Correlation Part"
        grad_n = self.__grad_n_no_exchange(n, c)
        grad_c = self.__grad_c_no_exchange(n, c)
        grad = numpy.hstack([grad_n, grad_c.ravel()])
        return grad
                                                                           
    def _compute_hcore_energy(self):
        "1-Electron Energy"
        H = self._H.copy()
        D = self._current_density.matrix()
        E = 2.0 * self._current_density.compute_1el_energy(D, H)
        return E
                                                                           
    def _compute_hartree_energy(self):
        "2-Electron Energy: Hartree"
        D = self._current_density.matrix()
        E = 2.0 * self._current_density.compute_2el_energy(D, D, type='j')
        return E


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


    def _pack(self, n, c):
        "Pack the n and c parameters into one hypervector"
        x   = numpy.hstack([n, c.ravel()])
        return x


    def _unpack(self, x):
        "Unpack n and c parameters from the hypervector"
        dim = self._bfs.nbf()
        n   = x[:dim]
        c   = x[dim:].reshape(dim,dim)
        return n, c


    # --- Private Interface --- #

    def __common_init(self, wfn, xc_functional, V_ext, guess):#OK--->JK!!!
        "Initialize"
        # Wavefunction
        self._wfn = wfn
        # Basis set name
        self._bfs_name = wfn.basisset().name
        # XC functional
        self._xc_functional = xc_functional

        # Molecule                                                                             
        self._mol = self._wfn.molecule()
        # Basis set
        self._bfs = self._wfn.basisset()
        # Integral calculator
        self._mints = psi4.core.MintsHelper(self._bfs)
        # JK object
        self._jk = psi4.core.JK.build(self._bfs, jk_type="Direct")
        self._jk.set_memory(int(5e8))
        self._jk.initialize()
        self._xc_functional.set_jk(self._jk)
        self._xc_functional.set_wfn(self._wfn)
        # Nuclear repulsion energy
        self._e_nuc = self._mol.nuclear_repulsion_energy()
        # Number of electron pairs
        self._np = self._wfn.nalpha()
                                                                                               
        ### Constant AO matrices
        # SCF LCAO-MO coefficients
        self._Ca= self._wfn.Ca_subset("AO","ALL").to_array(dense=True)
        # Overlap integrals and orthogonalizer
        self._S = self._wfn.S().to_array(dense=True)
        self._X = Density.  orthogonalizer(self._S)
        self._Y = Density.deorthogonalizer(self._S)
        #self._Y = numpy.linalg.inv(self._X)
        # Hcore matrix
        self._H = self._wfn.H().to_array(dense=True)
        # External potential
        self._V_ext = V_ext
        if V_ext is not None and guess == 'current':
           raise ValueError(" External potential can only be set for 'HCore' guess!")

        ### Current OPDM and total energy                                                                
        if guess == 'current':  # Guess based on current density matrix from the input wavefunction
           assert self._V_ext is None
           E = self._wfn.energy()
           D = self._wfn.Da().to_array(dense=True)
        elif guess == 'hcore':  # Guess based on one-electron Hamiltonian
           if self._V_ext is not None: self._H += self._V_ext
           e, c = numpy.linalg.eigh(numpy.linalg.multi_dot([self._X, self._H, self._X]))
           c = c[:,::-1]
           c = numpy.dot(self._X, c)
           E = 2.0*e.sum() + self._e_nuc
           D = numpy.zeros((self._bfs.nbf(), self._bfs.nbf()), numpy.float64)
           for i in range(self._np):
               D += numpy.outer(c[:,i], c[:,i])
        else:
           raise ValueError("Only 'Current' or 'HCore' ODPM's are supported as starting points.")
        self._current_energy  = E
        self._current_density = Density(D, self._jk)
        # Natural Orbitals
        self._current_occupancies, self._current_orbitals = \
                 Density.natural_orbitals(self._current_density.matrix(), self._S, self._Ca, orthogonalize_mo=True,
                                       order='descending', no_cutoff=0.0, renormalize=False, return_ao_orthogonal=False)
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




class DMFT_NC(DMFT):
    def __init__(self, wfn, xc_functional, v_ext, guess):
        super(DMFT_NC, self).__init__(wfn, xc_functional, v_ext, guess)

    # --- Implementation (Public) --- #

    @staticmethod
    def name(): return "DMFT with Gradient Projection on D Set: Optimization in N and C set."

    @property
    def abbr(self): return "DMFT-NC"


    # --- Implementation (Protected Interface) --- #


    def _minimizer(self, x):
        "Minimizer function: Total Energy"
        n, c = self._unpack(x)
        E_H  = self._compute_no_exchange_energy()
        E_XC = self._xc_functional.energy(n, c)
        E    = E_H + E_XC
        return E

    def _guess(self):
        "Initial guess"
        x = self._pack(self._current_occupancies, self._current_orbitals) 
        return x

    def _density(self, x):
        "1-particle density matrix in AO basis"
        n, c = self._unpack(x)
        n, c = density_matrix_projection(n, c, self._S, self._np, type='d')
        self._current_occupancies = n
        self._current_orbitals    = c
        D = self._current_density.generalized_density(n, c)
        self._current_density.set_D(D)
        x_new = self._pack(n, c)
        return x_new

    def _gradient(self, x):
        "Gradient"
        n, c = self._unpack(x)
        dE_n   , dE_c    = self._unpack(self._compute_no_exchange_gradient_nc(n, c))
        dE_n_xc, dE_c_xc = self._unpack(self._xc_functional      .gradient_nc(n, c))
        dE_c_xc = numpy.dot(self._X, dE_c_xc) # OAO basis
        dE_n += dE_n_xc
        dE_c += dE_c_xc
        gradient = self._pack(dE_n, dE_c)
        return gradient

    def _step(self, x1, x2):
        "Steepest-descents step. Back-transforms to OAO basis for simplicity"
        n1, c1 = self._unpack(x1)
        n2, c2 = self._unpack(x2)

        # transform C to OAO basis
        c1_ = numpy.dot(self._Y, c1)
        c2_ = numpy.dot(self._Y, c2)

        x_old_1 = self._pack(n1, c1_)
        x_old_2 = self._pack(n2, c2_)
                                                                                       
        gradient_1 = self._gradient(self._pack(n1, c1))
        gradient_2 = self._gradient(self._pack(n2, c2))
                                                                                       
        norm = numpy.linalg.norm(gradient_1 - gradient_2)
        g = numpy.dot(x_old_1 - x_old_2, gradient_1 - gradient_2) / norm**2
        #print(" G = ", g)
        g = abs(g)
        x_new = x_old_1 - g * gradient_1

        # transform back to AO basis
        n, c_ = self._unpack(x_new)
        x_new = self._pack(n, numpy.dot(self._X, c_))
        return x_new

    def _step_0(self, x0, g0):
        "Steepest-descents step. Back-transforms to OAO basis for simplicity"
        n2, c2 = self._unpack(x0)

        # transform C to OAO basis
        c2_ = numpy.dot(self._Y, c2)
        x_old_2 = self._pack(n2, c2_)

        # compute new guess
        gradient_2 = self._gradient(self._pack(n2, c2))
        x_new = x_old_2 - g0 * gradient_2
        
        # transform back to AO basis
        n, c_ = self._unpack(x_new)
        x_new = self._pack(n, numpy.dot(self._X, c_))
        return x_new




class DMFT_ProjD(DMFT):
    def __init__(self, wfn, xc_functional, v_ext, guess):
        super(DMFT_ProjD, self).__init__(wfn, xc_functional, v_ext, guess)
        raise NotImplementedError


    @staticmethod
    def name(): return "DMFT with Gradient Projection on D Set"

    @property
    def abbr(self): return "DMFT-ProjD"

class DMFT_ProjP(DMFT):
    def __init__(self, wfn, xc_functional, v_ext, guess):
        super(DMFT_ProjP, self).__init__(wfn, xc_functional, v_ext, guess)
        raise NotImplementedError


    @staticmethod
    def name(): return "DMFT with Gradient Projection on P Set"

    @property
    def abbr(self): return "DMFT-ProjP"

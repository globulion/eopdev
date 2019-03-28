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
import psi4
from abc import ABC, abstractmethod
from .functional import XCFunctional

__all__ = ["DMFT"]


class DMFT(ABC):
    """\
 The Density Matrix Functional Theory.
"""
    # Defaults
    default_xc_functional = 'hf'         # XC functional abbreviation
    default_algorithm     = 'proj-d'     # projected gradient algorithm abbreviation
    default_convergence   =  0.00001     # convergence in total energy [A.U.]
    default_maxiter       =  100         # maximum number of iterations
    default_verbose_run   =  True        # show additional information?
    default_g0            =  0.0001      # steepest-descents step scale in the first iteration

    def __init__(self, wfn, xc_functional):
        "Initialize"
        self.wfn = wfn
        self.basis_name = wfn.bassiset().name
        self.xc_functional = xc_functional
        self._current_energy = None

        # check if XC functional and algorithm are properly chosen
        if self.xc_functional.abbr.lower() != 'hf' and self.abbr.lower() == 'dmft-projd':
           print(""" Warning! The D-set is not Lipshitz with %s functional. 
 This will probably result in lack of convergence. Use P-set instead.""" % xc_functional.abbr.upper())

        super(DMFT, self).__init__()

    @classmethod
    def create(cls, wfn, xc_functional = default_xc_functional, 
                             algorithm = default_algorithm, **kwargs):
        """\
 Create DMFT solver. 
"""
        if   algorithm.lower() == 'proj-d': solver = DMFT_ProjD(wfn, xc_functional, **kwargs)
        elif algorithm.lower() == 'proj-p': solver = DMFT_ProjP(wfn, xc_functional, **kwargs)
        else: raise ValueError("Chosen algorithm is not available! Mistyped?")
        return solver

    def run(self, conv    = default_convergence, 
                  maxit   = default_maxiter    , 
                  verbose = default_verbose_run, 
                  g_0     = default_g0         , **kwargs):
        "Run the DMFT calculations"

        if verbose:
           print(" Running %s:%s/%s" % (self.abbr, xc_functional.abbr, self.basis_name))

        # Run!
        success = self._run_dmft(conv, maxit, verbose, g_0, **kwargs)
        return success

    @property
    def E(self): 
        "Total Energy"
        return self._current_energy

    @property
    @abstractmethod
    def abbr(self): 
        "Abbreviation symbol for the DMFT optimization algorithm"
        pass



    @staticmethod
    @abstractmethod
    def name(): pass


    # --- Protected Interface --- #

    def _run_dmft(self, conv, maxit, verbose, g_0, **kwargs):
        "DMFT Iterations"

        iteration = 0                                                             
        success = False

        # [1] Compute Guess
        x_0 = self._guess()
                                                                                  
        # [2] Starting energy
        self._current_energy = self._minimizer(x_0)
        E_old = self._current_energy
        if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (iteration, E_old))
                                                                                  
        # [3] First iteration
        iteration += 1
        x_old_2    = x_0
        grad_2     = self._gradient(x_0)
        x_old_1    = x_old_2 - g_0 * grad_2
        D          = self._density(x_old_1)
                                                                                  
        E_new = self._minimizer(x_old_1)
        if verbose: print(" @DMFT Iter %2d. E = %14.8f" % (iteration, E_new))
                                                                                  
        # [4] Further iterations
        iteration += 1
        stop       = False
        while stop is False:

            # [4.1] New guess
            x_new = self._step(x_old_1, x_old_2)
            D     = self._density(x_new)
                                                                                  
            # [4.2] Current energy
            E_new = self._minimizer(x_new)
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
            self._current_energy = E_old
                                                                                  
            x_old_2    = x_old_1.copy()
            x_old_1    = x_new  .copy()
        
        # [5] Finish
        if verbose and success:
           print(" DMFT iterations converged.")
           print(" Final Energy = %14.8f" % self._current_energy)

        if verbose and not success:
           print(" DMFT iterations did not converge.")

        return success


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
    def _step(self, x1, x2):
        "Steepest-descents step"
        pass








class DMFT_ProjD(DMFT):
    def __init__(self, wfn):
        super(DMFT_ProjD, self).__init__(wfn)


    @staticmethod
    def name(): return "DMFT with Gradient Projection on D Set"

    @property
    def abbr(self): return "DMFT-ProjD"

class DMFT_ProjP(DMFT):
    def __init__(self, wfn):
        super(DMFT_ProjP, self).__init__(wfn)
        raise NotImplementedError


    @staticmethod
    def name(): return "DMFT with Gradient Projection on P Set"

    @property
    def abbr(self): return "DMFT-ProjP"



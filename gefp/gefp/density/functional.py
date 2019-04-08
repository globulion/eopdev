#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT and DFT Functional module.
 Bartosz Błasiak, Gundelfingen, Mar 2019
"""

import os
import sys
import math
import numpy
from abc import ABC, abstractmethod
from .partitioning import Density
from .parameters import Guess
import oepdev
import psi4

__all__ = ["XCFunctional"]


class XCFunctional(ABC, Density):
    """\
 The Exchange-Correlation DMFT functional.
"""
    default = 'hf'

    def __init__(self, jk=None, wfn=None, ints=None):
        ABC.__init__(self)
        Density.__init__(self, None, jk)
        self._jk = self._global_jk
        self._wfn = wfn
        self._ints = ints
        #
        self._Ca   = None
        self._S    = None

    def set_jk(self, jk): 
        self._jk = jk
        self._global_jk = jk
    def set_wfn(self, wfn):
        self._wfn = wfn
        self._Ca  = wfn.Ca_subset("AO","ALL").to_array(dense=True)
        self._S   = wfn.S().to_array(dense=True)
    def set_ints(self, ints):
        self._ints = ints

    @classmethod
    def create(cls, name=default, **kwargs):
        """\
 Create a functional. Available functionals:
"""
        if   name.lower() == 'hf' :   xc_functional = HF_XCFunctional()
        elif name.lower() == 'mbb':   xc_functional = MBB_XCFunctional()
        else: raise ValueError("Chosen XC functional is not available! Mistyped?")
        return xc_functional

    @staticmethod
    @abstractmethod
    def fij(n): pass

    @staticmethod
    @abstractmethod
    def fij_1(n, m): pass

    @staticmethod
    @abstractmethod
    def name(): pass

    @abstractmethod
    def energy(self, n, c): 
        "Exchange-correlation energy"
        pass

    @abstractmethod
    def gradient_D(self, n, c): 
        "Gradient with respect to density matrix"
        pass

    @abstractmethod
    def gradient_P(self, n, c): 
        "Gradient with respect to P matrix"
        pass

    @abstractmethod
    def gradient_nc(self, n, c): 
        "Gradient with respect to N and C"
        pass


    @property
    def abbr(self): return default.upper()

    def _correct_negative_occupancies(self, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns



class HF_XCFunctional(XCFunctional):
    """
 Hartree-Fock Exchange-Correlation Functional
"""
    def __init__(self):
        super(HF_XCFunctional, self).__init__()
    @staticmethod
    def fij(n): return numpy.outer(n, n)
    @staticmethod
    def fij_1(n, m): 
        nn = len(n)                                  
        fij_1 = numpy.zeros((nn,nn), numpy.float64)
        for i in range(nn):
            for j in range(nn):
                v = 0.0
                if   (i==m) and (j!=m): v = n[j]
                elif (j==m) and (i!=m): v = n[i]
                elif (i==j==m)        : v = n[m]*2.0
                fij_1[i, j] = v
        return fij_1

    @staticmethod
    def name(): return "Hartree-Fock Functional XC for closed-shell systems"
    @property
    def abbr(self): return "HF"

    def energy(self, n, c, mode='scf-mo'): 
        "Exchange-correlation energy"
        D = self.generalized_density(n, c, 1.0)
        if mode.lower() == 'scf-mo':
           K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
           xc_energy = -numpy.dot(K, D).trace()
        elif mode.lower() == 'ao':
           xc_energy = -self.compute_2el_energy(D, D, type='k')
        else: raise ValueError("Only mode=ao or scf-mo is supported as for now. Mistyped?")
        return xc_energy

    def gradient_D(self, n, c): 
        "Gradient with respect to density matrix: MO-SCF basis"
        #D = self.generalized_density(n, c, 1.0)
        #K = self.generalized_JK(D, type='k')
        #gradient = Guess.create(matrix=-K) # ---> must be in MO basis!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #c_psi4 = self._wfn.Ca_subset("AO","ALL")
        #gradient = -oepdev.calculate_JK(self._wfn, c_psi4)[1].to_array(dense=True)

        D = Density.generalized_density(n, c)
        gradient_K  = -oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        gradient = Guess.create(matrix=gradient_K)
        return gradient

    def gradient_P(self, n, c):
        "Gradient with respect to P matrix"
        raise NotImplementedError("Gradient of HF energy are not implemented for P sets.")

    def gradient_nc(self, n, c):
        "Gradient with respect to N and C"
        nn = len(n)

        # gradient wrt n
        grad_n = numpy.zeros(nn, numpy.float64)
        D = Density.generalized_density(n, c)
        D = numpy.linalg.multi_dot([self._Ca.T, self._S, D, self._S, self._Ca])
        Kij = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)

        for m in range(nn):
            fij_1 = self.fij_1(n, m)
            grad_n[m] -=(numpy.dot(Kij, fij_1)).trace()

        # gradient wrt c
        grad_c = numpy.zeros((nn, nn), numpy.float64)

        self._jk.C_clear()
        fij = self.fij(n); I = numpy.identity(nn)
        for m in range(nn):
            Bm = numpy.linalg.multi_dot([c, numpy.diag(fij[:,m]), c.T])
            self._jk.C_left_add(psi4.core.Matrix.from_array(Bm, ""))
            self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._jk.compute()
                                                                                            
        for m in range(nn):
            Km = self._jk.K()[m].to_array(dense=True)
            grad_c[:,m] -= numpy.dot(Km, c[:,m])
        grad_c *= 4.0

        # pack
        gradient = Guess.create(grad_n, grad_c, None, 'nc')
        return gradient


 

class MBB_XCFunctional(XCFunctional):
    """
 Muller-Buijse-Baerends Exchange-Correlation Functional
"""
    def __init__(self):
        super(MBB_XCFunctional, self).__init__()
    @staticmethod
    def fij(n): 
        ns = numpy.sqrt(self.correct_negative_occupancies(n))
        return numpy.outer(ns, ns)
    @staticmethod
    def fij_1(n, m): 
        raise NotImplementedError("Derivatives of fij for MBB functional were not implemented.")
    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems"
    @property
    def abbr(self): return "MBB"

    def energy(self, n, c): 
        "Exchange-correlation energy"
        D = self.generalized_density(n, c, 0.5)
        xc_energy = -self.compute_2el_energy(D, D, type='k')
        return xc_energy

    def gradient_D(self, n, c): 
        "Gradient with respect to density matrix"
        raise NotImplementedError("Gradient of MBB XC energy are not implemented for D sets.")

    def gradient_P(self, n, c): 
        "Gradient with respect to P matrix"
        ns = self.correct_negative_occupancies(n)
        D = self.generalized_density(ns, c, 0.5)
        K = self.generalized_JK(D, type='k')
        gradient = -K
        return gradient

    def gradient_nc(self, n, c): 
        "Gradient with respect to N and C"
        raise NotImplementedError("Gradient of MBB XC energy are not implemented for nc sets.")

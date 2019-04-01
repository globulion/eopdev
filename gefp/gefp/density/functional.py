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

__all__ = ["XCFunctional"]


class XCFunctional(ABC, Density):
    """\
 The Exchange-Correlation DMFT functional.
"""
    default = 'hf'

    def __init__(self, jk=None):
        ABC.__init__(self)
        Density.__init__(self, None, jk)
        self._jk = self._global_jk

    def set_jk(self, jk): 
        self._jk = jk
        self._global_jk = jk

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
    def fij_m(n, m): pass

    @staticmethod
    @abstractmethod
    def name(): pass

    @classmethod
    @abstractmethod
    def energy(cls, n, c): 
        "Exchange-correlation energy"
        pass

    @classmethod
    @abstractmethod
    def gradient_D(cls, n, c): 
        "Gradient with respect to density matrix"
        pass

    @classmethod
    @abstractmethod
    def gradient_P(cls, n, c): 
        "Gradient with respect to P matrix"
        pass

    @classmethod
    @abstractmethod
    def gradient_nc(cls, n, c): 
        "Gradient with respect to N and C"
        pass


    @property
    def abbr(self): return default.upper()

    @classmethod
    def _correct_negative_occupancies(cls, n):
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
    def fij_m(n, m): 
        nn = len(n)                                  
        fij_m = numpy.zeros((nn,nn), numpy.float64)
        for i in range(nn):
            for j in range(nn):
                v = 0.0
                if   (i==m) and (j!=m): v = n[j]
                elif (j==m) and (i!=m): v = n[i]
                elif (i==j==m)        : v = n[m]*2.0
                fij_m[i, j] = v
        return fij_m

    @staticmethod
    def name(): return "Hartree-Fock Functional XC for closed-shell systems"
    @property
    def abbr(self): return "HF"

    @classmethod
    def energy(cls, n, c): 
        "Exchange-correlation energy"
        D = self.generalized_density(n, c, 1.0)
        xc_energy = -self.compute_2el_energy(D, D, type='k')
        return xc_energy

    @classmethod
    def gradient_D(cls, n, c): 
        "Gradient with respect to density matrix"
        D = self.generalized_density(n, c, 1.0)
        K = self.generalized_JK(D, type='k')
        gradient = -2.0 * K
        return gradient

    @classmethod
    def gradient_P(cls, n, c):
        "Gradient with respect to P matrix"
        raise NotImplementedError("Gradient of HF energy are not implemented for P sets.")

    @classmethod
    def gradient_nc(cls, n, c):#OK --> JK!!!
        "Gradient with respect to N and C"
        nn = len(n)

        # gradient wrt n
        grad_n = numpy.zeros(nn, numpy.float64)
        Kij = oepdev.calculate_Kij(self._wfn, c).to_array(dense=True)
        psi4.core.clean()

        for m in range(nn):
            fij_m = self.fij_1(n, m)
            grad_n[m] -=(numpy.dot(Kij, fij_m)).trace()

        # gradient wrt c
        grad_c = numpy.zeros((nn, nn), numpy.float64)

        self._jk.C_clear()                                           
        for m in range(nn):
            Bm = numpy.linalg.multi_dot([c, numpy,diag(fij[:,m]), c.T])
            self._jk.C_left_add(psi4.core.Matrix.from_array(Bm, ""))
            self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self._jk.compute()
                                                                                            
        for m in range(nn):
            Km = self._jk.K()[m].to_array(dense=True)
            grad_c[:,m] -= numpy.dot(Km, c[:,m])
        grad_c *= 4.0

        # pack
        gradient = numpy.hstack([grad_n, grad_c.ravel()])
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
    def fij_m(n, m): 
        raise NotImplementedError("Derivatives of fij for MBB functional were not implemented.")
    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems"
    @property
    def abbr(self): return "MBB"

    @classmethod
    def energy(cls, n, c): 
        "Exchange-correlation energy"
        D = self.generalized_density(n, c, 0.5)
        xc_energy = -self.compute_2el_energy(D, D, type='k')
        return xc_energy

    @classmethod
    def gradient_D(cls, n, c): 
        "Gradient with respect to density matrix"
        raise NotImplementedError("Gradient of MBB XC energy are not implemented for D sets.")

    @classmethod
    def gradient_P(cls, n, c): 
        "Gradient with respect to P matrix"
        ns = self.correct_negative_occupancies(n)
        D = self.generalized_density(ns, c, 0.5)
        K = self.generalized_JK(D, type='k')
        gradient = -K
        return gradient

    @classmethod
    def gradient_nc(cls, n, c): 
        "Gradient with respect to N and C"
        raise NotImplementedError("Gradient of MBB XC energy are not implemented for nc sets.")

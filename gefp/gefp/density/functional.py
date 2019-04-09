#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT Functional module.
 Bartosz Błasiak, Gundelfingen, Mar 2019
"""

import os
import sys
import math
import numpy
import oepdev
import psi4
from abc import ABC, abstractmethod
from .partitioning import Density
from .parameters import Guess

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

    def set_wfn(self, wfn):
        "Set wavefunction"
        self._wfn = wfn
        self._Ca  = wfn.Ca_subset("AO","ALL").to_array(dense=True)
        self._S   = wfn.S().to_array(dense=True)

    def set_jk(self, jk):
        "Set JK object for AO tensor contractions"
        self._jk = jk
        self._global_jk = jk

    def set_ints(self, ints):
        "Set integral transform object for MO tensor contractions"
        self._ints = ints

    @classmethod
    def create(cls, name=default, **kwargs):
        """\
 Create a density matrix exchange-correlation functional. 

 Available functionals:                              Sets:
  o 'HF'  - the Hartree-Fock functional (default)    D, NC
  o 'CHF' - the corrected Hartree-Fock functional    P
  o 'MBB' - the Muller-Buijse-Baerends functional    P
  o 'GU'  - the Goedecker-Urmigar functional         P

"""
        if   name.lower() == 'hf' :   xc_functional =  HF_XCFunctional()
        elif name.lower() == 'mbb':   xc_functional = MBB_XCFunctional()
        elif name.lower() == 'gu' :   xc_functional =  GU_XCFunctional()
        else: raise ValueError("Chosen XC functional is not available! Mistyped?")
        return xc_functional

    @staticmethod
    @abstractmethod
    def name(): pass

    @staticmethod
    @abstractmethod
    def fij(n): pass

    @staticmethod
    @abstractmethod
    def fij_1(n, m): pass


    def energy_D(self, x, mode): 
        "Exchange-correlation energy"
        raise NotImplementedError("%s energy is not implemented for D sets." % self.abbr.upper())

    def energy_P(self, x): 
        "Exchange-correlation energy"
        raise NotImplementedError("%s energy is not implemented for P sets." % self.abbr.upper())

    def gradient_D(self, x): 
        "Gradient with respect to density matrix"
        raise NotImplementedError("Gradient of %s energy is not implemented for D sets." % self.abbr.upper())

    def gradient_P(self, x): 
        "Gradient with respect to P matrix"
        raise NotImplementedError("Gradient of %s energy is not implemented for P sets." % self.abbr.upper())

    def gradient_nc(self, x): 
        "Gradient with respect to N and C"
        raise NotImplementedError("Gradient of %s energy is not implemented for NC sets." % self.abbr.upper())


    @property
    def abbr(self): 
        "Functional name abbreviation"
        return default.upper()


    # ----> Protected Interface (utilities) <---- #

    def _correct_negative_occupancies(self, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns



class HF_XCFunctional(XCFunctional):
    """
 The Hartree-Fock Exchange-Correlation Functional
"""
    def __init__(self):
        super(HF_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Hartree-Fock Functional XC for closed-shell systems"

    @property
    def abbr(self): return "HF"

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

    def energy_D(self, x, mode='scf-mo'): 
        "Exchange-correlation energy"
        n, c = x.unpack()
        D = self.generalized_density(n, c, 1.0)
        if mode.lower() == 'scf-mo':
           K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
           xc_energy = -numpy.dot(K, D).trace()
        elif mode.lower() == 'ao':
           xc_energy = -self.compute_2el_energy(D, D, type='k')
        else: raise ValueError("Only mode=ao or scf-mo is supported as for now. Mistyped?")
        return xc_energy

    def gradient_D(self, x): 
        "Gradient with respect to density matrix: MO-SCF basis"
        D = Density.generalized_density(*x.unpack())
        gradient_K  = -oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        gradient = Guess.create(matrix=gradient_K)
        return gradient

    def gradient_nc(self, x):
        "Gradient with respect to N and C"
        n, c = x.unpack()
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
 The Muller-Buijse-Baerends Exchange-Correlation Functional.
"""
    def __init__(self):
        super(MBB_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems"

    @property
    def abbr(self): return "MBB"

    @staticmethod
    def fij(n): 
        ns = numpy.sqrt(self.correct_negative_occupancies(n))
        return numpy.outer(ns, ns)

    @staticmethod
    def fij_1(n, m): 
        raise NotImplementedError("Derivatives of fij for MBB functional were not implemented since are unnecessary.")

    def energy_P(self, x): 
        "Exchange-correlation energy: Practical expression is for P-sets."
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        xc_energy = -numpy.dot(K, P).trace()
        return xc_energy

    def gradient_P(self, x):
        "Gradient with respect to P matrix"
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        return Guess.create(matrix=-2.0*K)


    # ----> Additional Interface (illustrative, not practical) <---- #

    def energy_D(self, x, mode='scf-mo'): 
        "Exchange-correlation energy: Not useful, only for illustrative purpose."
        D = self.generalized_density(*x.unpack(), 0.5)
        if mode.lower() == 'scf-mo':
           K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
           xc_energy = -numpy.dot(K, D).trace()
        elif mode.lower() == 'ao':
           xc_energy = -self.compute_2el_energy(D, D, type='k')
        else: raise ValueError("Only mode=ao or scf-mo is supported as for now. Mistyped?")
        return xc_energy


class GU_XCFunctional(XCFunctional):
    """
 The Goedecker-Urmigar Exchange-Correlation Functional.
"""
    def __init__(self):
        super(GU_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Goedecker-Urmigar XC Functional for closed-shell systems"

    @property
    def abbr(self): return "GU"

    @staticmethod
    def fij(n): 
        ns = n.copy()
        ns[ns<0.0] = 0.0
        #ns = numpy.sqrt(self._correct_negative_occupancies(n))
        f = numpy.outer(ns, ns) - numpy.diag(n - n*n)
        return f
        
    @staticmethod
    def fij_1(n, m):
        "P-version"
        nn = len(n)
        d = numpy.zeros((nn,nn))
        d[m,m] -= 2.0 * n[m] - 4.0 * n[m]**3
        for i in range(nn):
            for j in range(nn):
                if i==m: d[i,j] += n[j]
                if j==m: d[i,j] += n[i]
        return d

    def energy_P(self, x):#TODO! -> do the integrals look OK?
        "Exchange-correlation energy: Practical expression is for P-sets."
        p, c = x.unpack()
        n = p**2
        f = self.fij(n)
        C = numpy.dot(self._Ca, c)
        K  = oepdev.calculate_JK(self._wfn, psi4.core.Matrix.from_array(C, ""))[1].to_array(dense=True)
        xc_energy = -numpy.dot(K, f).trace()
        return xc_energy

    def gradient_P(self, x):#TODO!
        "Gradient with respect to P matrix"
        p, c = x.unpack() # C: MO(SCF)-MO(new)
        nn=len(p)

        A_nj = numpy.zeros((nn,nn))
        for i in range(nn):
            for j in range(nn):
                A_nj[i,j] = self.fij_1(p, j)[i,j]

        An_bd = numpy.einsum("nj,bj,dj->nbd",A_nj,c,c)
        An_bd_psi = []
        for i in range(len(An_bd)):
            An_bd_psi.append(psi4.core.Matrix.from_array(An_bd[i].copy(), ""))
        del An_bd
       
        C_psi = psi4.core.Matrix.from_array(c.T, "") # MO(new)-MO(SCF)
        gradient = oepdev.calculate_der_D(self._wfn, self._ints, C_psi, An_bd_psi)[1].to_array(dense=True)
        return Guess.create(matrix=gradient)

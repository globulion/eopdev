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
import scipy.special
from abc import ABC, abstractmethod
from .partitioning import Density
from .parameters import Guess
from ..math.matrix import matrix_power, matrix_power_derivative

__all__ = ["XCFunctional"]


class XCFunctional(ABC, Density):
    """\
 The Exchange-Correlation DMFT functional.
"""
    # Default Functional
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

 Available functionals:                              Sets:     Analytic Derivatives:
  o 'HF'  - the Hartree-Fock functional (default)    D, NC     Yes
  o 'MBB' - the Muller-Buijse-Baerends functional    P         Yes
  o 'GU'  - the Goedecker-Urmigar functional         P         Approximate
  o 'BBC1'- the BBC1 functional                      P         No
  o 'BBC2'- the BBC2 functional                      P         No
  o 'MEDI'- the monotonous exponential decay         P         No
            of interpolates between MBB and
            MBB with zero exchange.                  
  o 'OEDI'- the oscillatory exponential decay        P         No
            of interpolates between MBB and
            MBB with zero exchange.                  
"""
        # Stored Parameterizations                                                                                    
        parameters_p2_medi = [-0.637709, -81.9733, -3.75826, -16.0301, -2.29926, 68.4848, -2.22155, 46.483, 0.918441]
        # A_V2_MEDI -- FCI/STO-3G
        # a0              = -0.637709        +/- 0.01133      (1.776%)
        # a1              = -81.9733         +/- 39.4         (48.06%)
        # a2              = -3.75826         +/- 1.416        (37.67%)
        # a3              = -16.0301         +/- 10.78        (67.25%)
        # a4              = -2.29926         +/- 0.9788       (42.57%)
        # b1              = 68.4848          +/- 36.61        (53.45%)
        # b2              = -2.22155         +/- 1.542        (69.39%)
        # b3              = 46.483           +/- 24.82        (53.4%)
        # b4              = 0.918441         +/- 0.3866       (42.09%)

        if   name.lower() == 'hf'   : xc_functional =        HF_XCFunctional()
        elif name.lower() == 'mbb'  : xc_functional =       MBB_XCFunctional()
        elif name.lower() == 'gu'   : xc_functional =        GU_XCFunctional()
        elif name.lower() == 'bbc1' : xc_functional =      BBC1_XCFunctional()
        elif name.lower() == 'bbc2' : xc_functional =      BBC2_XCFunctional()
        elif name.lower() == 'a1medi': xc_functional =    A_V1_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        elif name.lower() == 'a2medi': xc_functional =    A_V2_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        elif name.lower() == 'p2medi': xc_functional =    P_V2_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        elif name.lower() == 'par_p2medi': xc_functional =P_V2_MEDI_XCFunctional(kwargs['kmax'], coeff={'pade_coefficients':parameters_p2_medi})
        #elif name.lower() ==  'oedi': xc_functional =   AB_OEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        #elif name.lower() == 'pmedi': xc_functional = Pade_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
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

    def energy_pc(self, x): 
        "Exchange-correlation energy"
        raise NotImplementedError("%s energy is not implemented for PC sets." % self.abbr.upper())

    def gradient_D(self, x): 
        "Gradient with respect to density matrix"
        raise NotImplementedError("Gradient of %s energy is not implemented for D sets." % self.abbr.upper())

    def gradient_P(self, x): 
        "Exact analytical gradient with respect to P matrix"
        raise NotImplementedError("Exact analytical gradient of %s energy is not implemented for P sets." % self.abbr.upper())

    def gradient_P_approximate(self, x): 
        "Approximate analytical gradient with respect to P matrix"
        raise NotImplementedError("Approximate analytical gradient of %s energy is not implemented for P sets." % self.abbr.upper())

    def gradient_P_numerical(self, x): 
        "Numerical gradient with respect to P matrix"
        raise NotImplementedError("Numerical gradient of %s energy is not implemented for P sets." % self.abbr.upper())

    def gradient_nc(self, x): 
        "Gradient with respect to N and C"
        raise NotImplementedError("Gradient of %s energy is not implemented for NC sets." % self.abbr.upper())

    def gradient_pc(self, x): 
        "Gradient with respect to P and C"
        raise NotImplementedError("Gradient of %s energy is not implemented for PC sets." % self.abbr.upper())


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

    def _compute_deriv_P(self, n, m, P, E0):
        step = 0.000008; h = 2.0*step
        P1= P.copy()
        P1[m,n] += step; P1[n,m] += step

        guess = Guess.create(matrix=P1)
        guess.update()
        E1 = self.energy_P(guess)
        der = (E1 - E0) / h
        return der

    def gradient_P_numerical(self, x):
        E0 = self.energy_P(x)
        P  = x.matrix()
        N  = P.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_P(i, j, P, E0)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)




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
        ns = n.copy()
        ns[ns<0.0] = 0.0
        ns = numpy.sqrt(ns)
        #ns = numpy.sqrt(self.correct_negative_occupancies(n))
        return numpy.outer(ns, ns)

    @staticmethod
    def fij_1(n, m): 
        raise NotImplementedError("Derivatives of fij for MBB functional were not implemented since are unnecessary.")

    def energy_P(self, x): 
        "Exchange-correlation energy: Practical expression is for P-sets."
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        xc_energy = -numpy.dot(K, P).trace()
        #p, c = x.unpack()
        #f = self.fij(p**2)
        #psi_f = psi4.core.Matrix.from_array(f, "")
        #psi_c = psi4.core.Matrix.from_array(c, "")
        #xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

    def energy_pc(self, x): 
        "Exchange-correlation energy: Practical expression is for PC-sets."
        xc_energy = self.energy_P(x)
        return xc_energy

    def gradient_P(self, x):
        "Gradient with respect to P matrix"
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        return Guess.create(matrix=-2.0*K)

    def gradient_pc(self, x):
        "Gradient with respect to PC matrix"
        P = x.matrix()
        p, c = x.unpack()
        K    = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        Kmm  = numpy.linalg.multi_dot([c.T, K, c]).diagonal()
        grad_p =-2.0 * Kmm
        grad_c =-4.0 * numpy.dot(K, c) * p
        return Guess.create(grad_p, grad_c, None, 'nc')



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
        ns = numpy.sqrt(ns)
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

    def energy_P(self, x):
        "Exchange-correlation energy: Practical expression is for P-sets."
        p, c = x.unpack()
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

    def __gradient_P_old(self, x):#deprecate!
        "Gradient with respect to P matrix"
        p, c = x.unpack() # C: MO(SCF)-MO(new)
        nn=len(p)

        A_nj = numpy.zeros((nn,nn))
        for i in range(nn):
            for j in range(nn):
                A_nj[i,j] = self.fij_1(p, j)[i,j]

        An_bd = numpy.einsum("nj,bj,dj->nbd",A_nj.T,c,c)
        An_bd_psi = []
        for i in range(len(An_bd)):
            An_bd_psi.append(psi4.core.Matrix.from_array(An_bd[i].copy(), ""))
        del An_bd

        C_psi = psi4.core.Matrix.from_array(c.T, "") # MO(new)-MO(SCF)
        gradient = oepdev.calculate_der_D(self._wfn, self._ints, C_psi, An_bd_psi).to_array(dense=True)
        return Guess.create(matrix=gradient)

    def gradient_P(self, x):
        "Approximate gradient with respect to P matrix"
        X = psi4.core.Matrix.from_array(x.matrix(), "")
        gradient = -2.0*oepdev.calculate_JK_r(self._wfn, self._ints, X)[1].to_array(dense=True)

        p, c = x.unpack() # C: MO(SCF)-MO(new)
        nn=len(p)
        s = 2.0 * p * (2.0 * p*p - 1.0)
        C = numpy.dot(self._Ca, c)
        #aaai = numpy.einsum("ijkl,ia,ja,ka,lb->ab", self._ao_eri, C, C, C, C)
        #aaa = aaai.sum(axis=1) * s
        #s = aaa
        s = numpy.einsum("ijkl,ia,ja,ka,la->a", self._ao_eri, C, C, C, C) * s
        gradient -= Density.generalized_density(s, c)
        return Guess.create(matrix=gradient)

class BBC1_XCFunctional(XCFunctional):
    """
 The BBC1 Exchange-Correlation Functional.
"""
    def __init__(self):
        super(BBC1_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "BBC1 XC Functional for closed-shell systems"

    @property
    def abbr(self): return "BBC1"

    @staticmethod
    def fij(n): 
        f = MBB_XCFunctional.fij(n) * BBC1_XCFunctional.__pc(n)
        return f

    @staticmethod
    def __pc(n):                                   
        "Phase correction according to BBC1 functional"
        m = n.copy(); m.fill(0.0)
        m[numpy.where(n>=0.5)] =  1.0  # strong
        m[numpy.where(n< 0.5)] = -1.0  # weak
        pc = numpy.ones((len(n), len(n)))
        for i in range(len(n)):
            for j in range(len(n)):
                if i!=j:
                  if m[i] < 0 and m[j] < 0:
                    pc[i,j] = -1.0
        return pc

        
    @staticmethod
    def fij_1(n, m):
        "P-version"
        raise NotImplementedError

    def energy_P(self, x):
        "Exchange-correlation energy: Practical expression is for P-sets."
        p, c = x.unpack()
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

class BBC2_XCFunctional(BBC1_XCFunctional):
    """
 The BBC2 Exchange-Correlation Functional.
"""
    def __init__(self):
        super(BBC2_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "BBC2 XC Functional for closed-shell systems"

    @property
    def abbr(self): return "BBC2"

    @staticmethod
    def fij(n): 
        f = BBC1_XCFunctional.fij(n)
        f_hf = HF_XCFunctional.fij(n)                  
        m = n.copy(); m.fill(0.0)
        m[numpy.where(n>=0.5)] =  1.0  # strong
        m[numpy.where(n< 0.5)] = -1.0  # weak
        for i in range(len(n)):
            for j in range(len(n)):
                if i!=j:
                   if m[i] > 0 and m[j] > 0:
                      f[i,j] = f_hf[i,j]
        return f


class Interpolation_XCFunctional(XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: Interpolation Functionals

 They differ in the model for the interpolation decay.
"""
    def __init__(self, coeff, kmax=10):
        super(Interpolation_XCFunctional, self).__init__()
        self._coeff  = coeff
        self._kmax   = kmax

    @abstractmethod
    def compute_a0(self, n):
        "First coefficient in the interpolates"
        pass

    @abstractmethod
    def compute_ak(self, k, t):
        "Further coefficients in the interpolates"
        pass


    @abstractmethod
    def fij(self, n):
        pass

    @staticmethod
    def _fij_bbbk_1(n, k, eps=1.0e-20):
        "The BBB-(k) Functional"
        W = 1.0/float(k+1)
        K = float(k*(k+1))
        U = float(k+1)
        f = MBB_XCFunctional.fij(n)
        f = 2.0**W * ( f**U )
        de= ( n[:,numpy.newaxis]**K + n[numpy.newaxis,:]**K )**W
        for i in range(len(n)):
            for j in range(len(n)):
                de_ij = de[i,j]
                if de_ij < eps: f[i,j] = 0.0
                else: f[i,j] = f[i,j] / de_ij
        return f

    @staticmethod
    def _fij_bbbk_2(n, k):
        "The BBB-(k) Functional - 2-version"
        K = float(k)
        #W = k/(k+1.0)**2
        #X = 0.5 * (1.0 - k/(k+1))
        #ni= n**W
        #f = 2.0**(-K-1.0) * numpy.outer(n,n)**X * (ni[:,numpy.newaxis] + ni[numpy.newaxis,:])**(K+1.0)
        #p = K/(K+1)**1

        #p = 1.0 # K/(K+1.0)
        #ni = n**(p/(K+1.0))
        #f = 2.0**(-K-1.0) * numpy.outer(n,n)**(0.5*(1.0-p)) * (ni[:,numpy.newaxis] + ni[numpy.newaxis,:])**(K+1.0)

        #nn = len(n)
        #f = numpy.zeros((nn,nn))
        #for t in range(0, k+2):
        #    c = scipy.special.binom(k+1, t)
        #    ni = n**(float(k+1-t)/float(k+1))
        #    nj = n**(float(t)/float(k+1))
        #    f += c * numpy.outer(ni, nj)
        #f*= 2.0**(-K-1.0)

        #K = K - 1.0
        ni = n**(1.0/(K+1.0))
        f = 2.0**(-K-1.0) * (ni[:,numpy.newaxis] + ni[numpy.newaxis,:])**(K+1.0)
        return f

        
    @staticmethod
    def fij_1(n, m):
        raise NotImplementedError

    def energy_P(self, x):
        "Exchange-correlation energy: Practical expression is for P-sets."
        p, c = x.unpack()
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

    def gradient_P(self, x):
        "Approximate gradient with respect to P matrix"
        raise NotImplementedError


class MEDI_XCFunctional(Interpolation_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay.

 The decay in the interpolates is modelled by the monotonous decay

 a_k = a_0 exp(a_0 log k)
"""
    def __init__(self, coeff, kmax):
        super(MEDI_XCFunctional, self).__init__(coeff, kmax)

    @abstractmethod
    def _interpolate_function(self, n, k, **kwargs): pass

    @abstractmethod
    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        pass

    def fij(self, n): 
        "The MBB-MBB0 Interpolation Functional with Monotonous Exponential Decay (MBB/MEDI)"
        a0 = self.compute_a0(n)
        # First term
        f = a0 * MBB_XCFunctional.fij(n)
        #a_sum = a0
        # Other terms
        t = self.compute_t(n)
        for k in range(1, self._kmax+1):
            #ak = a0 * math.exp(k*math.log(1.0 - a0))
            #ak = self.compute_ak(k, self._coeff['a0'])
            ak = self.compute_ak(k, t)
            f += ak * self._interpolate_function()(n, k)
            #print(ak, a0)
            #Interpolation_XCFunctional._fij_bbbk(n, k, eps=1.0e-20)
            #a_sum += ak
        #print( " Sum of a: %13.4f" % a_sum)
        return f

class V1_MEDI_XCFunctional(MEDI_XCFunctional):
    def __init__(self, coeff, kmax):
        super(V1_MEDI_XCFunctional, self).__init__(coeff, kmax)

    def _interpolate_function(self): 
       return Interpolation_XCFunctional._fij_bbbk_1

    def compute_ak(self, k, t):
        return t * math.exp(k*math.log(1.0 - t))

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        return self._coeff['a0']


class V2_MEDI_XCFunctional(MEDI_XCFunctional):
    def __init__(self, coeff, kmax):
        super(V2_MEDI_XCFunctional, self).__init__(coeff, kmax)

    def _interpolate_function(self): 
       return Interpolation_XCFunctional._fij_bbbk_2

    def compute_ak(self, k, t):
        A = 100.0
        L =-0.4
        return A * scipy.special.gegenbauer(k, L)(t/2.0) * t**k

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        return self._coeff['t']

    def compute_ak_derivative_t(self, k, t):
        A = 100.0
        L =-0.4
        d = k*scipy.special.gegenbauer(k, L)(t/2.0) + L*t*scipy.special.gegenbauer(k-1, L+1.0)(t/2.0)
        d*= A*t**(k-1)
        return d





#class OEDI_XCFunctional(Interpolation_XCFunctional):
#    """
# The New Class of Exchange-Correlation Functionals: 
# Interpolation Functionals with Oscillatory Exponential Decay.
#
# The decay in the interpolates is modelled by the oscillatory decay
#
# a_k = a_0 exp(-k B ) cos(k C)
#"""
#    def __init__(self, coeff, kmax):
#        super(OEDI_XCFunctional, self).__init__(coeff, kmax)
#
#    @abstractmethod
#    def compute_c(self, n):
#        pass
#
#    def fij(self, n): 
#        "The MBB-MBB0 Interpolation Functional with Oscillatory Exponential Decay (MBB/OEDI)"
#        a0 = self.compute_a0(n)
#        C  = self.compute_c(n)
#        B  = self.compute_b(a0, C)
#        # First term
#        f = a0 * MBB_XCFunctional.fij(n)
#        # Other terms
#        for k in range(1,self._kmax+1):
#            ak = a0 * math.exp(-k*B) * math.cos(k*C)
#            f += ak * self._fij_bbbk(n, k, eps=1.0e-20)
#        return f
#
#    def compute_b(self, a0, C):
#        "Decay rate"
#        r  = 2./a0 - 1.0
#        b  = 2.0*r*math.cos(C) +  math.sqrt(2.0 * (2.0 - r*r + r*r*math.cos(2.0*C)))
#        b /= 2.0*(r-1.0)
#        #b = r * math.cos(C) - math.sqrt(1.0 + r*r * math.sin(C))
#        #b = abs(b/ (1.0 + r))
#        b  = math.log(b)
#        return b


class A_V1_MEDI_XCFunctional(V1_MEDI_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay.
"""
    def __init__(self, coeff, kmax):
        super(A_V1_MEDI_XCFunctional, self).__init__(coeff, kmax)

    @staticmethod
    def name(): return "Monotonous Exponential Decay of Interpolates XC Functional for closed-shell systems"

    @property
    def abbr(self): return "MEDI"

    def compute_a0(self, n):
        "First coefficient in the interpolates from Pade approximant of universal function"
        a0 = self._coeff['a0']
        return a0

class A_V2_MEDI_XCFunctional(V2_MEDI_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay.
"""
    def __init__(self, coeff, kmax):
        super(A_V2_MEDI_XCFunctional, self).__init__(coeff, kmax)

    @staticmethod
    def name(): return "Monotonous Exponential Decay of Interpolates XC Functional for closed-shell systems"

    @property
    def abbr(self): return "MEDI"

    def compute_a0(self, n):
        "First coefficient in the interpolates from Pade approximant of universal function"
        return 1.0

    def energy_P_costly(self, x): 
        "Exchange-correlation energy: This expression is more costly"
        # contribution from MBB term
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        xc_energy = -numpy.dot(K, P).trace()
        del K
        p, c = x.unpack()
        t = self.compute_t(p*p, c)
        # contributions from interpolates
        for k in range(1, self._kmax+1):
            K = float(k)
            a_k = self.compute_ak(k, t)
            a_k/= 2.0**(K+1.0)
            for n in range(0, k+2):
                N = float(n)
                c_nk = a_k * scipy.special.binom(k+1, n)
                p_nk = 2.0 * N / (K+1.0)
                A_nk = matrix_power(P,       p_nk)
                B_nk = matrix_power(P, 2.0 - p_nk)
                A    = psi4.core.Matrix.from_array(A_nk, "")
                KA   = oepdev.calculate_JK_r(self._wfn, self._ints, A)[1].to_array(dense=True)
                xc_energy -= c_nk * numpy.dot(KA, B_nk).trace()

        return xc_energy

    def gradient_P(self, x):
        "Exacr Analytical Gradient with respect to P matrix"
        # contribution from the MBB term
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        grad = -2.0*K
        del K
        p, c = x.unpack()
        t = self.compute_t(p*p, c)
        # contributions from the interpolates                                                  
        for k in range(1, self._kmax+1):
            K = float(k)
            a_k = self.compute_ak(k, t)
            a_k/= 2.0**(K+1.0)
            for n in range(0, k+2):
                N = float(n)
                c_nk = a_k * scipy.special.binom(k+1, n)
                p_nk = 2.0 * N / (K + 1.0)
                A_nk = matrix_power(P,       p_nk)
                B_nk = matrix_power(P, 2.0 - p_nk)
                dA   = matrix_power_derivative(P,       p_nk, step=0.000001, approx=False)
                dB   = matrix_power_derivative(P, 2.0 - p_nk, step=0.000001, approx=False)
                A    = psi4.core.Matrix.from_array(A_nk, "")
                B    = psi4.core.Matrix.from_array(B_nk, "")
                KA   = oepdev.calculate_JK_r(self._wfn, self._ints, A)[1].to_array(dense=True)
                KB   = oepdev.calculate_JK_r(self._wfn, self._ints, B)[1].to_array(dense=True)
                a = numpy.einsum("ijkl,kl->ij", dA, KB)
                b = numpy.einsum("ijkl,kl->ij", dB, KA)
                grad-= c_nk * (a+b)
        return Guess.create(matrix=grad)


    def gradient_P_approximate(self, x):
        "Approximate analytical gradient with respect to P matrix"
        # contribution from the MBB term
        P = x.matrix()
        K  = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)
        grad = -2.0*K
        del K
        p, c = x.unpack()
        t = self.compute_t(p*p, c)
        # contributions from the interpolates                                                  
        for k in range(1, self._kmax+1):
            K = float(k)
            a_k = self.compute_ak(k, t)
            a_k/= 2.0**(K+1.0)
            for n in range(0, k+2):
                N = float(n)
                c_nk = a_k * scipy.special.binom(k+1, n)
                p_nk = 2.0 * N / (K + 1.0)
                A_nk = matrix_power(P,       p_nk)
                B_nk = matrix_power(P, 2.0 - p_nk)
                dA_nk= matrix_power_derivative(P,       p_nk, step=0.000001, approx=True)
                dB_nk= matrix_power_derivative(P, 2.0 - p_nk, step=0.000001, approx=True)
                A    = psi4.core.Matrix.from_array(A_nk, "")
                B    = psi4.core.Matrix.from_array(B_nk, "")
                KA   = oepdev.calculate_JK_r(self._wfn, self._ints, A)[1].to_array(dense=True)
                KB   = oepdev.calculate_JK_r(self._wfn, self._ints, B)[1].to_array(dense=True)
                grad-= c_nk * (numpy.dot(dA_nk, KB) + numpy.dot(KA, dB_nk))
        return Guess.create(matrix=grad)



#class AB_OEDI_XCFunctional(OEDI_XCFunctional):
#    """
# The New Class of Exchange-Correlation Functionals: 
# Interpolation Functionals with Oscillatory Exponential Decay.
#"""
#    def __init__(self, coeff, kmax):
#        super(AB_OEDI_XCFunctional, self).__init__(coeff, kmax)
#
#    @staticmethod
#    def name(): return "Oscillatory Exponential Decay of Interpolates XC Functional for closed-shell systems"
#
#    @property
#    def abbr(self): return "OEDI"
#
#    def compute_a0(self, n):
#        "First coefficient in the interpolates from Pade approximant of universal function"
#        a0 = self._coeff['a0']
#        return a0
#
#    def compute_c(self, n):
#        "Oscillation period"
#        c = self._coeff['c']
#        return c

class P_V2_MEDI_XCFunctional(A_V2_MEDI_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay: Pade approximant for universal function.
"""
    def __init__(self, coeff, kmax):
        super(P_V2_MEDI_XCFunctional, self).__init__(coeff, kmax)

    @staticmethod
    def name(): return "Monotonous Exponential Decay of Interpolates XC Functional for closed-shell systems"

    @property
    def abbr(self): return "MEDI"

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        a0, a1, a2, a3, a4, b1, b2, b3, b4 = self._coeff['pade_coefficients']
        ns = n.copy(); ns[ns<0.0] = 0.0
        i_n = (ns*(1.0 - ns)).sum()
        i_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - i_n
        i_n = abs(i_n)
        i_d = abs(i_d)
        N   = self._wfn.nalpha()
        x = math.log(i_d/N + 1.0)
        y = math.log(2.0 * i_n/N) if i_n > 0.0 else 1.e10
        t = 0.5*math.erf( (a0 + a1*x + a2*y + a3*x*y + a4*y*y)/(1.0+b1*x + b2*y + b3*x*y + b4*y*y) ) + 0.5
        return t





#class Pade_MEDI_XCFunctional(MEDI_XCFunctional):
#    """
# The New Class of Exchange-Correlation Functionals: 
# Interpolation Functionals with Monotonous Exponential Decay: Pade approximant for universal function.
#"""
#    def __init__(self, coeff, kmax):
#        super(Pade_MEDI_XCFunctional, self).__init__(coeff, kmax)
#
#    @staticmethod
#    def name(): return "Monotonous Exponential Decay of Interpolates XC Functional for closed-shell systems"
#
#    @property
#    def abbr(self): return "MEDI"
#
#    #f(x,y) = 0.5+0.5*erf( (a0 + a*x + b*y + c*x*y + d*y*y)/(1.0 + e*x + f*y + g*x*y + h*y*y) )
#    #a               = -163.922         +/- 18.85        (11.5%)
#    #b               = -10.348          +/- 1.888        (18.24%)
#    #c               = 160.828          +/- 20.28        (12.61%)
#    #d               = 1.02675          +/- 0.2923       (28.47%)
#    #e               = -240.342         +/- 9.075        (3.776%)
#    #f               = 17.9345          +/- 1.195        (6.664%)
#    #g               = 55.4836          +/- 23.3         (42%)
#    #h               = -3.54743         +/- 0.134        (3.776%)
#    #a0              = 10.5995          +/- 1.064        (10.04%)
#
#    def compute_a0(self, n):
#        "First coefficient in the interpolates from Pade approximant of universal function"
#        # Pade parameters
#        A = self._coeff['A']
#        B = self._coeff['B']
#        A0, A1, A2, A3, A4 = A
#        B1, B2, B3, B4     = B
#
#        # Dynamic and non-dynamic correlation
#        ns = n.copy(); ns[ns<0.0] = 0.0
#        I_n = (ns*(1.0 - ns)).sum()
#        I_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - I_n
#        S   = numpy.sqrt(ns).sum()
#        N   = ns.sum()
#        #print(ns)
#        #print(N, S, I_n, I_d)
#
#        # Universal phase space
#        x = math.log(I_d/S + 1.0)
#        c = abs(2.0 * I_n/S)
#        if c < 1.0e-20: y = 1e100
#        else: y =-math.log(c)
#
#        # Coefficient
#        a_0 = (A0 + A1*x + A2*y + A3*x*y + A4*y*y)/\
#              (1.0+ B1*x + B2*y + B3*x*y + B4*y*y)
#        a_0 = 0.500*(math.erf(a_0) + 1.0)
#        return a_0
#
#    def __compute_a0_old(self, n):
#        "First coefficient in the interpolates from Pade approximant of universal function"
#        # Pade parameters
#        A = self._coeff['A']
#        B = self._coeff['B']
#        A0, A1, A2, A3, A5, A6 = A
#        B1, B2, B3, B5, B6     = B
#
#        # Dynamic and non-dynamic correlation
#        ns = n.copy(); ns[ns<0.0] = 0.0
#        I_n = (ns*(1.0 - ns)).sum()
#        I_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - I_n
#        S   = numpy.sqrt(ns).sum()
#        N   = ns.sum()
#        #print(ns)
#        #print(N, S, I_n, I_d)
#
#        # Universal phase space
#        x = math.log(I_d/S + 1.0)
#        y =-math.log(abs(2.0 * I_n/S))
#
#        # Coefficient
#        a_0 = (A0 + A1*x + A2*y + A3*x*y + A5*y*y + A6*x*y*y)/\
#              (1.0+ B1*x + B2*y + B3*x*y + B5*y*y + B6*x*y*y)
#        a_0 = 0.500*(math.erf(a_0) + 1.0)
#        #print(a_0)
#        return a_0

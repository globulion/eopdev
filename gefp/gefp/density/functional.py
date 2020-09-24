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
from .opdm import Density
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
  o 'IDF1'- interpolating density matrix functional  P         No
            ???
"""
        if   name.lower() == 'hf'   : xc_functional =        HF_XCFunctional()
        elif name.lower() == 'mbb'  : xc_functional =       MBB_XCFunctional()
        elif name.lower() == 'mbb-1': xc_functional =     MBB_1_XCFunctional()
        elif name.lower() == 'mbb-2': xc_functional =     MBB_2_XCFunctional()
        elif name.lower() == 'gu'   : xc_functional =        GU_XCFunctional()
        elif name.lower() == 'chf'  : xc_functional =       CHF_XCFunctional()
        elif name.lower() == 'ohf'  : xc_functional =       OHF_XCFunctional()
        elif name.lower() == 'bbc1' : xc_functional =      BBC1_XCFunctional()
        elif name.lower() == 'bbc2' : xc_functional =      BBC2_XCFunctional()
        elif name.lower() == 'idf1' : xc_functional =      IDF1_XCFunctional(kwargs['parameters']) #TODO
        elif name.lower() == 'a1medi': xc_functional =    A_V1_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        elif name.lower() == 'a2medi': xc_functional =    A_V2_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        elif name.lower() == 'p2medi': xc_functional =    P_V2_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        #elif name.lower() == 'par_p2medi': xc_functional =P_V2_MEDI_XCFunctional(kwargs['kmax'], coeff={'pade_coefficients':parameters_p2_medi})
        #elif name.lower() ==  'oedi': xc_functional =   AB_OEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        #elif name.lower() == 'pmedi': xc_functional = Pade_MEDI_XCFunctional(kwargs['coeff'], kwargs['kmax'])
        else: raise ValueError("Chosen XC functional is not available! Mistyped?")
        return xc_functional

    @staticmethod
    @abstractmethod
    def name(): pass

    @property
    def abbr(self): 
        "Functional name abbreviation"
        return default.upper()

    @staticmethod
    @abstractmethod
    def fij(n): pass

    @staticmethod
    def fij_1(n, m): 
        raise NotImplementedError("Derivatives of fij for %s functional were not implemented since they are unnecessary." %\
                                  self.abbr.upper())

    def energy_D(self, x, mode):
        "Exchange-correlation energy"
        raise NotImplementedError("%s energy is not implemented for D sets." % self.abbr.upper())

    def energy_D_add(self, x): #ADD
        "Exchange-correlation energy"
        n, c = x.unpack()
        f = self.fij(n)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy


    def energy_P(self, x):
        "Exchange-correlation energy"
        p, c = x.unpack()
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

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

    def gradient_D_numerical(self, x): #ADD
        "Numerical gradient with respect to P matrix"
        E0 = self.energy_D_add(x)
        D  = x.matrix()
        N  = D.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_D(i, j, D, E0)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)


    def gradient_nc(self, x): 
        "Gradient with respect to N and C"
        raise NotImplementedError("Gradient of %s energy is not implemented for NC sets." % self.abbr.upper())

    def gradient_pc(self, x): 
        "Gradient with respect to P and C"
        raise NotImplementedError("Gradient of %s energy is not implemented for PC sets." % self.abbr.upper())



    # ----> Protected Interface (utilities) <---- #

    def _correct_negative_occupancies(self, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns

    def _compute_deriv_P(self, n, m, P, E0):
        step = 0.0000002; h = 2.0*step
        P1= P.copy()
        P1[m,n] += step; P1[n,m] += step

        guess = Guess.create(matrix=P1)
        guess.update()
        E1 = self.energy_P(guess)
        der = (E1 - E0) / h
        return der 

    def _compute_deriv_D(self, n, m, D, E0): #ADD
        step = 0.0000002; h = 2.0*step
        D1= D.copy()
        D1[m,n] += step; D1[n,m] += step

        guess = Guess.create(matrix=D1)
        guess.update()
        E1 = self.energy_D_add(guess)
        der = (E1 - E0) / h
        return der 




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
        gradient_K  = -2.0 * oepdev.calculate_JK_r(self._wfn, self._ints, 
                                                   psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        gradient = Guess.create(matrix=gradient_K)
       #Dp = numpy.linalg.pinv(D)
       #print(D.diagonal())
       #print(Dp.diagonal())
       #print("Łojojoj!! ", self.energy_D(x), (0.5*Dp @ gradient_K).trace())
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

    def energy_P(self, x): 
        "Exchange-correlation energy"
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

class MBB_1_XCFunctional(MBB_XCFunctional):
    """
 The Muller-Buijse-Baerends Exchange-Correlation Functional + first term from unfolding
"""
    def __init__(self):
        super(MBB_1_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems + 1st unfolding term"

    @property
    def abbr(self): return "MBB-1"

    @staticmethod
    def fij(n, eps=1e-20): 
        a = 0.2 # MBB part in the functional: adjust manually
        b = 1.0 - a
        s2 = math.sqrt(2.0)
        fij_mbb = MBB_XCFunctional.fij(n)
        f2 = n*n 
        fij_hf = numpy.outer(n,n)
        d = numpy.sqrt( f2[:,numpy.newaxis] + f2[numpy.newaxis,:] )
        f = fij_hf.copy()
        for i in range(len(n)):
            for j in range(len(n)):
                d_ij = d[i,j]
                if abs(d_ij) < eps: f[i,j] = a*fij_mbb[i,j]
                else: f[i,j] = a*fij_mbb[i,j] + b*s2 * f[i,j] / d_ij
        return f

    def energy_P(self, x):
        "Exchange-correlation energy"
        p, c = x.unpack()
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        return xc_energy

    def gradient_P(self, x):
        "Gradient with respect to P matrix"
        raise NotImplementedError

    def gradient_pc(self, x):
        "Gradient with respect to PC matrix"
        raise NotImplementedError

class MBB_2_XCFunctional(MBB_1_XCFunctional):
    """
 The Muller-Buijse-Baerends Exchange-Correlation Functional + first and second terms from unfolding
"""
    def __init__(self):
        super(MBB_2_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems + 1st unfolding term"

    @property
    def abbr(self): return "MBB-2"

    @staticmethod
    def fij(n, eps=1e-20): 
        a = 0.2 # MBB_1 part in the functional: adjust manually
        b = 1.0 - a
        s2 = 2.0**(1./3.)
        fij_mbb_1 = MBB_1_XCFunctional.fij(n)
        f6 = n**6
        fij_q = numpy.outer(n,n)**(3./2.)
        d = ( f6[:,numpy.newaxis] + f6[numpy.newaxis,:] )**(1./3.)
        f = fij_q.copy()
        for i in range(len(n)):
            for j in range(len(n)):
                d_ij = d[i,j]
                if abs(d_ij) < eps: f[i,j] = a*fij_mbb_1[i,j]
                else: f[i,j] = a*fij_mbb_1[i,j] + b*s2 * f[i,j] / d_ij
        return f

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

    def gradient_P_approximate_old(self, x):#deprecate!
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

    def gradient_P_approximate(self, x):
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

class CHF_XCFunctional(XCFunctional):
    """
 The Corrected Hartree-Fock Exchange-Correlation Functional.
"""
    def __init__(self):
        super(CHF_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Corrected Hartree-Fock XC Functional for closed-shell systems"

    @property
    def abbr(self): return "CHF"

    @staticmethod
    def fij(n): 
        x = n*(1.0 - n)
        f = numpy.outer(n,n)
        f+= numpy.sqrt(abs(numpy.outer(x,x)))
        return f

class OHF_XCFunctional(XCFunctional):
    """
 The Overrepulsive Hartree-Fock Exchange-Correlation Functional.
"""
    def __init__(self):
        super(OHF_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Overrepulsive Hartree-Fock XC Functional for closed-shell systems"

    @property
    def abbr(self): return "OHF"

    @staticmethod
    def fij(n): 
        raise NotImplementedError
 
    def energy_D_add(self, x): #ADD
        "Exchange-correlation energy"
        n, c = x.unpack()
        p    = numpy.sqrt(abs(n))
        P    = c @ numpy.diag(p) @ c.T
      # P = self.generalized_density(n, c, 0.5)
        y    = Guess.create(p, c, P, "matrix")
        xc_energy = self.energy_P(y)
        return xc_energy
    
    def energy_P(self, x):
        "Exchange-correlation energy of OHF functional"
        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T
        P = x.matrix()
        Kd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        Jp= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[0].to_array(dense=True)
        
       #f = numpy.outer(n,n)
       #psi_f = psi4.core.Matrix.from_array(f, "")
       #psi_c = psi4.core.Matrix.from_array(c, "")
       #xc_energy_hf_test = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);

        xc_energy =      -numpy.dot(Kd, D).trace()
       #print("Test= ", xc_energy, xc_energy_hf_test)
        xc_energy-= 2.0 * numpy.dot(Jd, D).trace()
        xc_energy+= 2.0 * numpy.dot(Jp, P).trace()
        return xc_energy
       

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

class Interpolating_XCFunctional(XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: Interpolation Functionals: Abstract base.
"""
    def __init__(self, parameters):
        super(Interpolating_XCFunctional, self).__init__()
        self.parameters = parameters

    @staticmethod 
    def phi(n):
        return (n*n).sum()

    @staticmethod
    def fij(n): 
        raise NotImplementedError

    @staticmethod
    def fij_deprecate(n, z, eps=1.0e-20): 
        "z - unfolding parameter"
        raise NotImplementedError
        k = z
        W = 1.0/float(k+1.0)
        K = float(k*(k+1.0))
        U = float(k+1.0)
        f = MBB_XCFunctional.fij(n)
        f = 2.0**W * ( f**U )
        de= ( n[:,numpy.newaxis]**K + n[numpy.newaxis,:]**K )**W
        for i in range(len(n)):
            for j in range(len(n)):
                if i==j: f[i,j] = n[i]
                else:
                  de_ij = de[i,j]
                  if de_ij < eps: f[i,j] = 0.0
                  else: f[i,j] = f[i,j] / de_ij
        return f




class IDF1_XCFunctional(Interpolating_XCFunctional):
    """
 4-Component Interpolation Functional with MBB lower and HF upper bound.
"""
    def __init__(self, parameters):
        super(IDF1_XCFunctional, self).__init__(parameters)

    @staticmethod
    def name(): return "IDF(MBB/HF) XC functional for closed-shell systems"

    @property
    def abbr(self): return "IDF1"

    @staticmethod
    def aN(phi, N):
        x = 2.0 * phi / N
        a0 = (0.0417644 * N - 0.00223544 * N*N) / (1.0 -1.04842 * N + 0.292082 * N*N)
        a1 = (1.71556   * N + 0.00409398 * N*N) / (1.0 -1.42459 * N + 1.02131  * N*N)
        a2 = (0.00024548* N - 0.00000168 * N*N) / (1.0 -0.923738* N + 0.215108 * N*N)
        a = a0 + a1 * x + a2 * x*x
       #print("a= ",a)
        return a
    @staticmethod
    def xiN(an, N):
        kn = abs(an * (N-1.0)/N)
        kns= math.sqrt(kn)
        xi = (2.0 + 2.0*math.exp(kn)*(2.0*kns * scipy.special.dawsn(kns) - 1.0) ) / an
       #print(kn, kns, scipy.special.dawsn(kns), an, math.exp(kn))
       #print("xi= ",xi)
       #exit()
        return xi

    def energy_D_add(self, x): #ADD
        "Exchange-correlation energy"
        n, c = x.unpack()
        p    = numpy.sqrt(n)
        P    = c @ numpy.diag(p) @ c.T
        # P = self.generalized_density(n, c, 0.5)
        y    = Guess.create(p, c, P, "matrix")
        xc_energy = self.energy_P(y)
        return xc_energy

    
    def energy_P(self, x, n_quad=12, eps=1.0e-19):
        "Exchange-correlation energy of IDF functional"
        p, c = x.unpack()
       #print("p= ",p*p)
        n    = p*p ; n[numpy.where(n < 0.0)] = 0.0
        N    = n.sum().round()*2.0
        Phi  = Interpolating_XCFunctional.phi(n)
        a_n  = IDF1_XCFunctional.aN(Phi, N)
        xi_n = IDF1_XCFunctional.xiN(a_n, N)

        kn   = abs(a_n * (N-1.0)/N)
        si, ci = scipy.special.shichi(kn)
        A = (si + ci -numpy.log(kn) -numpy.euler_gamma - xi_n*a_n) / a_n

       #print("N= ", N)
        x, w = scipy.special.roots_legendre(n_quad); wnorm = 1./w.sum()
        o    = 0.5*(1.0 + x)
        o    = psi4.core.Vector.from_array(o)
        w    = psi4.core.Vector.from_array(w)

        # f and g matrices
        f = []
        g = []
        zeta = 0.00001 # adjust!!!

        L = numpy.dot(self._Ca.T, self._S)                        
        R = numpy.dot(numpy.linalg.inv(numpy.dot(L.T, L)), L.T).T
       #R = L.copy()

        for i in range(n_quad):
            omega = o.get(i)
            z = -zeta * numpy.log(omega)
           #print("Quad: ", i+1, omega, z)
            fx = abs(n**(z*(z+1.0)))
            vij = abs( (fx[numpy.newaxis,:] + fx[:,numpy.newaxis])**(1.0/(z+1.0)) )
            wij = abs( (2.0)**(1.0/(z+1.0)) * numpy.outer(n,n)**((z+1.0)/2.0) )
            Fmo = vij.copy()
           #print(n)
           #print(vij, wij)
            for I in range(len(wij)):
                for J in range(len(wij)):
                    if I==J: Fmo[I,J] = n[I]
                    else:
                      dv_ij = vij[I,J]
                      if abs(dv_ij) < eps: Fmo[I,J] = 0.0
                      else: Fmo[I,J] = wij[I,J] / dv_ij

            Gmo = scipy.linalg.fractional_matrix_power(Fmo.real, 0.5)

            Fao = psi4.core.Matrix.from_array(R.T @ Fmo.real @ R)
            Gao = psi4.core.Matrix.from_array(R.T @ Gmo.real @ R)
            f.append(Fao)
            g.append(Gao)

        # current density matrix (AO)
        Dmo = c @ numpy.diag(n) @ c.T

        if 1:
           Dao = numpy.linalg.multi_dot([R.T, Dmo, R])
        else:
           Dao = self._Ca.T @ self._S @ Dmo @ self._S @ self._Ca
        Dao = psi4.core.Matrix.from_array(Dao)
        
        xc_energy = oepdev.calculate_idf_xc_energy(self._wfn, Dao, f, g, w, o, N, a_n, xi_n, A, wnorm)
        return xc_energy


class JKOnly_Interpolating_XCFunctional(Interpolating_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: Interpolation Functionals for JK-Only ansatz

 They differ in the model for the interpolation decay.
 Each functional has to provide its own coefficients, 
 handled by 'coeff' dictionary in the constructor.

 Eg.: coeff = {'coefficient_name': coefficient_object}
"""
    def __init__(self, coeff, kmax=10):
        super(JKOnly_Interpolating_XCFunctional, self).__init__(parameters=None) # change it because parameters are still used but in the old interface
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
        ni = n**(1.0/(K+1.0))
        f = 2.0**(-K-1.0) * (ni[:,numpy.newaxis] + ni[numpy.newaxis,:])**(K+1.0)

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
        return f

    def gradient_P(self, x):
        "Approximate gradient with respect to P matrix"
        raise NotImplementedError


class MEDI_XCFunctional(JKOnly_Interpolating_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay.

 The decay in the interpolates is modelled by the monotonous decay

 a_k = a_0 exp(k ln{1 - a0})

 Functional coefficients:
  o 'a0'  - first term in the interpolates coefficient 
            (the one that multiplies MBB functional term - for k=0).
            Parameter 'a0' has to be between 0.0 and 1.0. The smaller 'a0',
            the more terms need to be taken (larger kmax).
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
            ak = self.compute_ak(k, t)
            f += ak * self._interpolate_function()(n, k)
            #a_sum += ak
        #print( " Sum of a: %13.4f" % a_sum)
        return f

class V1_MEDI_XCFunctional(MEDI_XCFunctional):
    """\
 Version 1 of MEDI functional. It is a proper functional of density matrix. 
 However, it cannot be represented explicitely in terms of density matrix.

 Functional coefficients:
  o 't'  - t parameter in Gegenbauer Series. Parameter t to be between 0.0 and 1.0.
           The larger t, the more terms should be included (higher kmax).
"""
    def __init__(self, coeff, kmax):
        super(V1_MEDI_XCFunctional, self).__init__(coeff, kmax)

    def _interpolate_function(self): 
       return JKOnly_Interpolating_XCFunctional._fij_bbbk_1

    def compute_ak(self, k, t):
        return t * math.exp(k*math.log(1.0 - t))

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        return self._coeff['a0']


class V2_MEDI_XCFunctional(MEDI_XCFunctional):
    """\
 Version 2 of MEDI functional. It is a proper functional of density matrix. 
 Also, it can be represented explicitely in terms of density matrix.

 Functional coefficients:
  o 't'  - t parameter in Gegenbauer Series. Parameter t to be between 0.0 and 1.0.
           The larger t, the more terms should be included (higher kmax).

"""

    # Gegenbauer Series Parameters
    parameter_A = 100.0
    parameter_L =-0.4 

    def __init__(self, coeff, kmax):
        super(V2_MEDI_XCFunctional, self).__init__(coeff, kmax)

    def _interpolate_function(self): 
       return JKOnly_Interpolating_XCFunctional._fij_bbbk_2

    def compute_ak(self, k, t):
        "Computes ak coefficients in Gegenbauer series"
        A = V2_MEDI_XCFunctional.parameter_A
        L = V2_MEDI_XCFunctional.parameter_L
        return A * scipy.special.gegenbauer(k, L)(t/2.0) * t**k

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        return self._coeff['t']

    def compute_ak_derivative_t(self, k, t):
        "Derivative of ak coefficient wrt t"
        A = V2_MEDI_XCFunctional.parameter_A
        L = V2_MEDI_XCFunctional.parameter_L
        d = k*scipy.special.gegenbauer(k, L)(t/2.0) + L*t*scipy.special.gegenbauer(k-1, L+1.0)(t/2.0)
        d*= A*t**(k-1)
        return d


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
        "Exact Analytical Gradient with respect to P matrix"
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


class P_V2_MEDI_XCFunctional(A_V2_MEDI_XCFunctional):
    """
 The New Class of Exchange-Correlation Functionals: 
 Interpolation Functionals with Monotonous Exponential Decay: 2D Pade Approximant for Universal Function.

 Functional coefficients:
  o 'pade_coefficients'  - 2D Pade Coefficients that fit the 't' parameter.
"""
    def __init__(self, coeff, kmax):
        super(P_V2_MEDI_XCFunctional, self).__init__(coeff, kmax)

    @staticmethod
    def name(): return "Pade-Guided Monotonous Exponential Decay of Interpolates XC Functional for closed-shell systems"

    @property
    def abbr(self): return "P-MEDI"

    def compute_t(self, n=None, c=None):
        "t parameter in Gegenbauer numbers generator"
        pade = self._coeff['pade_coefficients']
        ns = n.copy(); ns[ns<0.0] = 0.0
        i_n = (ns*(1.0 - ns)).sum()
        i_d = numpy.sqrt(abs(ns*(1.0 - ns))).sum() / 2.0 - i_n
        i_n = abs(i_n)
        i_d = abs(i_d)
        N   = self._wfn.nalpha()
        x = math.log(i_d/N + 1.0)
        y = math.log(2.0 * i_n/N) if i_n > 0.0 else 1.e10
        #t = 0.5*math.erf( (a0 + a1*x + a2*y + a3*x*y + a4*y*y)/(1.0+b1*x + b2*y + b3*x*y + b4*y*y) ) + 0.5
        t = pade.value(x, y)
        return t





#class OEDI_XCFunctional(JKOnly_Interpolating_XCFunctional):
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

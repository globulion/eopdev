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
        #
        self._ffstep = 0.00000002

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

    def set_ffstep(self, step):
        "Set the finite difference step for numerical differentiation wrt D or P matrices"
        self._ffstep = step

    def load(self):
        "Load miscellanea"
        return

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
        elif name.lower() == 'mbb-anti': xc_functional =     MBB_ANTI_XCFunctional()
        elif name.lower() == 'gu'   : xc_functional =        GU_XCFunctional()
        elif name.lower() == 'chf'  : xc_functional =       CHF_XCFunctional()
        elif name.lower() == 'ohf'  : xc_functional =       OHF_XCFunctional()
        elif name.lower() == 'hf-mbb' : xc_functional =       HF_MBB_XCFunctional()
        elif name.lower() == 'ohf-mbb' : xc_functional =     OHF_MBB_XCFunctional()
        elif name.lower() == 'smbb-up' : xc_functional =     SMBB_UP_XCFunctional()
        elif name.lower() == 'bbc1' : xc_functional =      BBC1_XCFunctional()
        elif name.lower() == 'bbc2' : xc_functional =      BBC2_XCFunctional()
        elif name.lower() == 'apsg' : xc_functional =      APSG_XCFunctional()
        elif name.lower() == 'idf1' : xc_functional =      IDF1_XCFunctional(kwargs['parameters']) #TODO
        elif name.lower() == 'idfw' : xc_functional =      IDFW_XCFunctional(kwargs['omega']) #TODO
        elif name.lower() == 'idf-jk' : xc_functional =      IDF_JK_XCFunctional(kwargs['omega']) #TODO
        elif name.lower() == 'idf-jk+' : xc_functional =      IDF_JKPlus_XCFunctional(kwargs['omega']) #TODO
        elif name.lower() == 'idf-jk-di' : xc_functional =      IDF_JK_Di_XCFunctional(kwargs['ngrid'], kwargs['constant']) #TODO
        elif name.lower() == 'idf-alpha' : xc_functional =      IDF_Alpha_XCFunctional(kwargs['alpha'], kwargs['ngrid']) #TODO
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


    def potential_D(self, x):
        "Exchange-Correlation Potential in MO basis: D-sets"
        raise NotImplementedError("XC Potential in MO basis for D-sets is not implemented.")

    def potential_P(self, x):
        "Exchange-Correlation Potential in MO basis: P-sets"
        raise NotImplementedError("XC Potential in MO basis for P-sets is not implemented.")


    # ----> Protected Interface (utilities) <---- #

    def _correct_negative_occupancies(self, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns

    def _compute_deriv_P(self, n, m, P, E0):
        step = self._ffstep; h = 2.0*step
        P1= P.copy()
        P1[m,n] += step; P1[n,m] += step

        guess = Guess.create(matrix=P1)
        guess.update()
        E1 = self.energy_P(guess)
        der = (E1 - E0) / h
        return der 

    def _compute_deriv_D(self, n, m, D, E0): #ADD
        step = self._ffstep; h = 2.0*step
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


    def potential_D(self, x):
        "Exchange-Correlation Potential in MO basis"
        D = Density.generalized_density(*x.unpack())
        K = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        V = -K
        return Guess.create(matrix=V)


 

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


    def potential_P(self, x):
        "Exchange-Correlation Potential in MO basis"
        p, c = x.unpack()
        n = p*p
        y = Guess.create(n=n, c=c, t="matrix")
        V = -self.gradient_D_numerical(y).matrix()
        return Guess.create(matrix=V)  # does return -self.gradient_D_numerical(y) work too?

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
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems + 1st and 2nd unfolding terms"

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

    def load(self):
        print(" Computing AO ERIs")
        self._ao_eri = numpy.array(psi4.core.MintsHelper(self._wfn.basisset()).ao_eri())

    def gradient_P  (self, x):
        ""
        p, c = x.unpack()
       # for i in range(p.size):
       #     if p[i]   > 0.99999999999:
       #        p[i] = p[i] - (i+1)*0.0000000000000000000006
       #     elif p[i] < 0.00000000001:
       #        p[i] = p[i] + (i+1)*0.0000000000000000000006

        P_psi = psi4.core.Matrix.from_array(x.matrix(), "")

        KP = oepdev.calculate_JK_r(self._wfn, self._ints, P_psi)[1].to_array(dense=True)
        KP = c.T @ KP @ c # NO basis

        s = 2.*p*(1.0 - 2.*p*p)

        AJ = numpy.zeros((p.size,p.size))
        aJ = numpy.zeros((p.size,p.size))
        AKL=-self.fij(p*p)
       #aKL=-numpy.diag(2.*KP.diagonal() - s)
        aKL= numpy.diag(s) - 2.0* p[:,numpy.newaxis] 

        psi_p  = psi4.core.Vector.from_array(p  , "")
        psi_c  = psi4.core.Matrix.from_array(c  , "")
        psi_AJ = psi4.core.Matrix.from_array(AJ , "")
        psi_AKL= psi4.core.Matrix.from_array(AKL, "")
        psi_aJ = psi4.core.Matrix.from_array(aJ , "")
        psi_aKL= psi4.core.Matrix.from_array(aKL, "")

        gradient   = oepdev.calculate_de_apsg(self._wfn, self._ints, psi_p, psi_AJ, psi_AKL, psi_aJ, psi_aKL, psi_c).to_array(dense=True)

        # transform to MO basis
        gradient = c @ gradient @ c.T

        return Guess.create(matrix=gradient)


    def gradient_P_i  (self, x):
        ""
        p, c = x.unpack()
       # for i in range(p.size):
       #     if p[i] > 0.999999:
       #        p[i] = p[i] - (i+1)*0.000000006
       #     elif p[i] < 0.00001:
       #        p[i] = p[i] + (i+1)*0.000000006

       #P_psi = psi4.core.Matrix.from_array(x.matrix(), "")

       #KP = oepdev.calculate_JK_r(self._wfn, self._ints, P_psi)[1].to_array(dense=True)
       #KP = c.T @ KP @ c # NO basis

        s = 2.*p*(1.0 - 2.*p*p)

        C = self._Ca @ c

        eri_no = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", self._ao_eri, C, C, C, C)
        t = numpy.einsum("mmmm->m", eri_no) * s

        f = self.fij(p*p)

        A = numpy.einsum("mjjn,mj->mn", eri_no, f); A -= A.T
        A*=-2.0

        kp = numpy.einsum("mjjm,j->m", eri_no, p)

        # gradient in NO basis
        gradient = numpy.zeros((p.size,p.size))
        for i in range(p.size):
           #gradient[i,i] = -2.0*KP[i,i] + t[i]
            gradient[i,i] = -2.0*kp[i  ] + t[i]
            for j in range(p.size): 
                if i!=j:
                   dp = p[i] - p[j]
                   if abs(dp) > 1.e-80:
                      gradient[i,j] = A[i,j] / dp 
                      gradient[j,i] = gradient[i,j]

        # transform to MO basis
        gradient = c @ gradient @ c.T

        return Guess.create(matrix=gradient)


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

class HF_MBB_XCFunctional(XCFunctional):
    """
 The Corrected Hartree-Fock Exchange-Correlation Functional.
"""
    def __init__(self):
        super(HF_MBB_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "50-50 Mix of Hartree-Fock and MBB XC Functional for closed-shell systems"

    @property
    def abbr(self): return "HF-MBB"

    @staticmethod
    def fij(n): 
        f = 0.5 * (HF_XCFunctional.fij(n) + MBB_XCFunctional.fij(n))
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

class OHF_MBB_XCFunctional(XCFunctional):
    """
 The 50-50 Mix of Overrepulsive Hartree-Fock and MBB Exchange-Correlation Functional.
"""
    def __init__(self):
        super(OHF_MBB_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "50-50 Overrepulsive Hartree-Fock MBB XC Functional for closed-shell systems"

    @property
    def abbr(self): return "OHF-MBB"

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
        n = p*p; # print("N= %14.5f"%n.sum())
        D = c @ numpy.diag(n) @ c.T
        P = x.matrix()
        Kd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        Jp= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[0].to_array(dense=True)
        
        xc_energy =      -numpy.dot(Kd, D).trace()
       #print("Test= ", xc_energy, xc_energy_hf_test)
        xc_energy-= 2.0 * numpy.dot(Jd, D).trace()
        xc_energy+= 2.0 * numpy.dot(Jp, P).trace()

        # MBB part
        f = MBB_XCFunctional.fij(n)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy += oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);
        xc_energy*= 0.5
        return xc_energy

class MBB_ANTI_XCFunctional(XCFunctional):
    """
Antisymmetrized MBB XC Functional
"""
    def __init__(self):
        super(MBB_ANTI_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Antisymmetrized MBB XC Functional for closed-shell systems"

    @property
    def abbr(self): return "MBB-ANTI"

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
        "Exchange-correlation energy of MBB-ANTI functional"
        p, c = x.unpack()
        n = p*p; # print("N= %14.5f"%n.sum())
        D = c @ numpy.diag(n) @ c.T
        P = x.matrix()

        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        Kd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[1].to_array(dense=True)
        Jp= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[0].to_array(dense=True)
        Kp= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P, ""))[1].to_array(dense=True)

        xc_energy = (Jp @ P).trace() - (Jd @ D).trace() - (Kd @ D).trace() - (Kp @ P).trace()
        xc_energy*= 0.5

        return xc_energy


class SMBB_UP_XCFunctional(XCFunctional):
    """
"""
    def __init__(self):
        super(SMBB_UP_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "Upper Bound for SMBB XC Functional for closed-shell systems"

    @property
    def abbr(self): return "SMBB-UP"

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
        n = p*p; # print("N= %14.5f"%n.sum())
        D = c @ numpy.diag(n) @ c.T
        P = x.matrix()

        # 2xMBB part
        f = MBB_XCFunctional.fij(n)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = 1.0*oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);

        # 4-rank tensor part
       #f = numpy.sqrt(f)        
        F = c @ f @ c.T
        Kf= oepdev.calculate_JK_rb(self._wfn, self._ints, psi4.core.Matrix.from_array(F, ""))[1].to_array(dense=True)
        
        de = numpy.dot(Kf, F).trace()
        xc_energy += de

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
        super(BBC2Log_XCFunctional, self).__init__()

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

class APSG_XCFunctional(XCFunctional):
    """
 The APSG Exchange-Correlation Functional.

 Perfect pairing (PP) is assumed.
"""
    def __init__(self):
        super(APSG_XCFunctional, self).__init__()

    @staticmethod
    def name(): return "APSG XC Functional for closed-shell systems"

    @property
    def abbr(self): return "APSG"

    @staticmethod
    def fij(n): raise NotImplementedError

    def load(self):
        self._NG= self._wfn.nalpha()
        self._N = self._NG * 2
        r1 = numpy.arange(self._NG)
        r2 = self._N - r1  - 1
        self._geminal_idx = numpy.array(list(zip(r1,r2)), int)
        #
        mints = psi4.core.MintsHelper(self._wfn.basisset())
        self._eri_ao = numpy.array(mints.ao_eri(), numpy.float64)


    @staticmethod
    def geminals(n):
        ns = abs(n); nmo = n.size
        n_strong = ns[numpy.where(ns>=0.5)]

        NG = n_strong.size
        geminal_idx = []
        for i in range(NG):
            J = 2*NG - i - 1
            indices = [i,J]
            geminal_idx.append(numpy.array(indices,int))

        print(geminal_idx[0])
        return geminal_idx
            

    def geminal_index(self, i):
        G = -1
        for g in range(self._NG):
            idx = self._geminal_idx[g]
            if i in idx:
               G = g 
               break
        return G

    def fij_JKL(self, n): 
        ns = abs(n); p = numpy.sqrt(ns)

        geminal_idx = self._geminal_idx

        f_J = numpy.zeros((n.size,n.size))
        f_K = numpy.zeros((n.size,n.size))
        f_L = numpy.zeros((n.size,n.size))

        for i in range(n.size):
          G_i = self.geminal_index(i)
          if G_i >= 0:
            I_G_i = self._geminal_idx[G_i]

            phase_i = 1.0 if i == I_G_i[0] else -1.0

            for j in range(n.size):
              G_j = self.geminal_index(j)
              if G_j >= 0:
                 I_G_j = self._geminal_idx[G_j]

                 phase_j = 1.0 if j == I_G_j[0] else -1.0

                 if G_i != G_j:
                    f_J[i,j]+= ns[i] * ns[j]
                    f_K[i,j]-= ns[i] * ns[j]
                 if G_i == G_j:
                    f_L[i,j]+= p[i] * p[j] * phase_i * phase_j

        f_J *= 2.0
        f_K *= 1.0
        return f_J, f_K, f_L

    def fj_JKL(self, n): 
        ns = abs(n); p = numpy.sqrt(ns)

        geminal_idx = self._geminal_idx

        f_J = numpy.zeros((n.size,n.size))
        f_K = numpy.zeros((n.size,n.size))
        f_L = numpy.zeros((n.size,n.size))

        for i in range(n.size):
          G_i = self.geminal_index(i)
          if G_i >= 0:
            I_G_i = self._geminal_idx[G_i]

            phase_i = 1.0 if i == I_G_i[0] else -1.0

            for j in range(n.size):
              G_j = self.geminal_index(j)
              if G_j >= 0:
                 I_G_j = self._geminal_idx[G_j]

                 phase_j = 1.0 if j == I_G_j[0] else -1.0

                 if G_i != G_j:
                    f_J[i,j]+= p[i] * ns[j]
                    f_K[i,j]-= p[i] * ns[j]
                 if G_i == G_j:
                    f_L[i,j]+= p[j] * phase_i * phase_j

        f_L *= 2.0
        f_J *= 4.0
        f_K *= 2.0
        return f_J.T, f_K.T, f_L.T



    def energy_P(self, x):
        "Exchange-correlation energy of APSG functional"
        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T

        f_J, f_K, f_L = self.fij_JKL(n)

        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D,""))[0].to_array(dense=True)

        psi_fJ  = psi4.core.Matrix.from_array(f_J      , "")
        psi_fKL = psi4.core.Matrix.from_array(f_K + f_L, "")
        psi_c   = psi4.core.Matrix.from_array(c        , "")

        if 1:

           xc_energy = oepdev.calculate_e_apsg(self._wfn, self._ints, psi_fJ, psi_fKL, psi_c)
           xc_energy-= 2.0 * numpy.dot(Jd, D).trace()

        else:
           C = self._Ca @ c
           eri_no = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", self._eri_ao, C, C, C, C)

           xc_energy = numpy.einsum("ijij,ij", eri_no, f_K+f_L) + numpy.einsum("iijj,ij", eri_no, f_J)
           xc_energy-= 2.0 * numpy.dot(Jd, D).trace()

        return xc_energy

    def gradient_P(self, x):
        "Analytical gradient"
        p, c = x.unpack()
        P = x.matrix()
        P2= P @ P

        A_J, A_K, A_L = self.fij_JKL(p*p)
        a_J, a_K, a_L = self.fj_JKL(p*p)

       # A_J /= 2.0
       # a_J /= 2.0
       # A_K /= 2.0
       # a_K /= 2.0
       # A_L /= 1.0
       # a_L /= 1.0

        A_KL = A_K + A_L
        a_KL = a_K + a_L

       # A_KL /= 4.0
       # a_KL /= 4.0
       # A_J  /= 4.0
       # a_J  /= 4.0

        psi_AJ  = psi4.core.Matrix.from_array(A_J      , "")                                                                              
        psi_AKL = psi4.core.Matrix.from_array(A_KL     , "")
        psi_aJ  = psi4.core.Matrix.from_array( a_J      , "")
        psi_aKL = psi4.core.Matrix.from_array(a_KL     , "")
        psi_c   = psi4.core.Matrix.from_array(c        , "")
        psi_p   = psi4.core.Vector.from_array(p        , "")
                                                                                                                                          
        JP= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P2, ""))[0].to_array(dense=True)

        if 1:
           gradient   = oepdev.calculate_de_apsg(self._wfn, self._ints, psi_p, psi_AJ, psi_AKL, psi_aJ, psi_aKL, psi_c).to_array(dense=True)
           gradient   = c @ gradient @ c.T

        else:
           C = self._Ca @ c
           eri_no = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", self._eri_ao, C, C, C, C)

           g_mm = numpy.einsum("mjmj,jm->m", eri_no, a_KL)
           g_mm+= numpy.einsum("mmjj,jm->m", eri_no, a_J)

           t_mn = numpy.einsum("mjnj,mj->mn", eri_no, A_KL)
           t_mn+= numpy.einsum("mnjj,mj->mn", eri_no, A_J )
           t_mn-= t_mn.T
           g_mn = numpy.zeros((p.size, p.size))
           for i in range(p.size):
               g_mn[i,i] = g_mm[i]
               for j in range(i):
                 d_ij = p[i] - p[j]
                 if abs(d_ij) > 1.e-20:
                   g_mn[i,j] = 2./d_ij * t_mn[i,j]
                   g_mn[j,i] = g_mn[i,j]

           gradient = c @ g_mn @ c.T
          

        gradient_H = JP @ P
        gradient_H+= gradient_H.T
        gradient_H*= 4.0

        gradient-= gradient_H

        return Guess.create(matrix=gradient)
 

class IDFW_XCFunctional(XCFunctional):
    """
 The IDF test version of functional for particular omega value.
"""
    def __init__(self, omega):
        super(IDFW_XCFunctional, self).__init__()
        self._omega = omega
        self._p     = (3.+numpy.cos(omega/2.))/4.
        self._q     = (3.-numpy.cos(omega/2.))/4.
        print(" P = %14.5f  Q = %14.5f" % (self._p, self._q))


    @staticmethod
    def name(): return "IDFW Test Version of XC Functional"

    @property
    def abbr(self): return "IDFW"

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

   #def _make_unit(self, w, i, j):
   #    N = self._wfn.nmo()
   #    U = numpy.identity(N)        
   #    c = numpy.cos(w); s = numpy.sin(w)
   #    U[i,i] = U[j,j] = c
   #    U[i,j] = s
   #    U[j,i] =-s 
   #    return U

    def _make_U(self, w):
        a = 0.001 * 10000
        a = 0.0
        N = self._wfn.nmo()
        na= self._wfn.nalpha()
        A = numpy.zeros((N,N))
        s = numpy.sin(w/2.)
        v = a * s*s
        for i in range(N):
         #if i >= na:
            for j in range(i):
             #if ((i <na) and (j <na)):
                A[i,j] = v
                A[j,i] =-v
       #f = numpy.linalg.det(A)
       #print(f)
       #print(A); exit()
        U = scipy.linalg.expm(A)
       #print(U); exit()
        return U
        #c = numpy.cos(w); s = numpy.sin(w)
        #u = numpy.array([c, s, -s, c]).reshape(2,2)
        #if N==2: return u
        #else:                               
        #  U = numpy.identity(N)
        #  for i in range(N):
        #      for j in range(N-1):
        #        if i<j:
        #          u = self._make_unit(w, i, j)
        #          U = U @ u
        #  return U

    def load(self):
        "Load miscellanea"
        self._U = self._make_U(self._omega)
        return
   
    def gradient_P(self, x):
        "Analytical gradient"
        p, c = x.unpack()
        P = x.matrix()
        P2= P @ P

        U = self._U @ c # OK
       #U = c @ self._U

        Q1 = scipy.linalg.fractional_matrix_power(P, 2.*self._p - 1.0).real
        Q2 = scipy.linalg.fractional_matrix_power(P, 2.*self._q - 1.0).real

       #X = self._U @ numpy.diag(p**(2.*self._p)) @ self._U.T
       #Y = self._U @ numpy.diag(p**(2.*self._q)) @ self._U.T

       #X = c @ X @ c.T
       #Y = c @ Y @ c.T

        X = U @ numpy.diag(p**(2.*self._p)) @ U.T
        Y = U @ numpy.diag(p**(2.*self._q)) @ U.T

        J = oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(P2, ""))[0].to_array(dense=True)
        Jx= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(X, ""))[0].to_array(dense=True)
        Ky= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(Y, ""))[1].to_array(dense=True)

        Jx= c.T @ Jx @ c
        Ky= c.T @ Ky @ c

       #U = c @ self._U

        W1= U @ Jx @ U.T
        W2= U @ Ky @ U.T

        QW1= Q1 @ W1 + W1 @ Q1
        QW2= Q2 @ W2 + W2 @ Q2

       #print(W1);print(J); exit()

        gradient = 2.*(2.*self._p) * QW1 - (2.*self._q) * QW2

        gradient_H = J @ P
        gradient_H+= gradient_H.T
        gradient_H*= 4.0

        gradient-= gradient_H

       #print(gradient_H)
       #print(2.*(2.*self._p) * QW1)
       #exit()

        return Guess.create(matrix=gradient)
 
    def energy_P(self, x):
        "Exchange-correlation energy of OHF functional"
       #dim = self._wfn.nmo()
       #self._U     = numpy.identity(dim)
       #self._U.fill(numpy.sin(self._omega)**2)
       #for i in range(dim): self._U[i,i] = numpy.cos(self._omega)
       #self._U = self._make_U(self._omega)

        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T
       #P = x.matrix()

       #U = c.T @ self._U @ c

       #X = U @ numpy.diag(n**self._p) @ U.T
       #Y = U @ numpy.diag(n**self._q) @ U.T

       #X = c @ X @ c.T
       #Y = c @ Y @ c.T

        U = self._U @ c # OK
       #U = c @ self._U

        X = U @ numpy.diag(n**self._p) @ U.T
        Y = U @ numpy.diag(n**self._q) @ U.T

        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        Jx= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(X, ""))[0].to_array(dense=True)
        Ky= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(Y, ""))[1].to_array(dense=True)

       #Jx= c.T @ Jx @ c
       #Ky= c.T @ Ky @ c

       #Jx= U @ Jx @ U.T
       #Ky= U @ Ky @ U.T
       
        xc_energy =      -numpy.dot(Ky, Y).trace()
        xc_energy+= 2.0 * numpy.dot(Jx, X).trace()
        xc_energy-= 2.0 * numpy.dot(Jd, D).trace()
        return xc_energy

    def potential_P(self, x):
        "Exchange-Correlation Potential in MO basis"
        pref_p = 1./(2. * self._p)
        pref_q = 1./(2. * self._q)
       #pref_p = 1.0
       #pref_q = 0.5
        print(pref_p, pref_q)

        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T
        y = Guess.create(n, c, D, "matrix")
        V =  (self._gradient_D_numerical_unfolding_part_p(y).matrix()) * pref_p
        V+=  (self._gradient_D_numerical_unfolding_part_q(y).matrix()) * pref_q

        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        V-= Jd * 2.
        return Guess.create(matrix=V)  

    def _gradient_D_numerical_unfolding_part_p(self, x):
        "Numerical gradient with respect to P matrix"
        E0 = self.energy_D_part_p(x)
        D  = x.matrix()
        N  = D.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_D_part_p(i, j, D, E0)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)

    def _gradient_D_numerical_unfolding_part_q(self, x):
        "Numerical gradient with respect to P matrix"
        E0 = self.energy_D_part_q(x)
        D  = x.matrix()
        N  = D.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_D_part_q(i, j, D, E0)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)

    def _compute_deriv_D_part_p(self, n, m, D, E0): #ADD
        step = self._ffstep; h = 2.0*step
        D1= D.copy()
        D1[m,n] += step; D1[n,m] += step

        guess = Guess.create(matrix=D1)
        guess.update()
        E1 = self.energy_D_part_p(guess)
        der = (E1 - E0) / h
        return der 

    def _compute_deriv_D_part_q(self, n, m, D, E0): #ADD
        step = self._ffstep; h = 2.0*step
        D1= D.copy()
        D1[m,n] += step; D1[n,m] += step

        guess = Guess.create(matrix=D1)
        guess.update()
        E1 = self.energy_D_part_q(guess)
        der = (E1 - E0) / h
        return der 

    def energy_D_part_p(self, x):
        n, c = x.unpack()
       #n = p*p
        D = c @ numpy.diag(n) @ c.T

       #F = abs(numpy.outer(n,n))
       #X = self._U * F**self._p

       #self._U = numpy.identity(n.size)

       #X = self._U @ numpy.diag(n**self._p) @ self._U.T
       #X = c @ X @ c.T

        U = self._U @ c

        X = U @ numpy.diag(n**self._p) @ U.T

        Jx= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(X, ""))[0].to_array(dense=True)
       
        energy = 2.0 * numpy.dot(Jx, X).trace()
        return energy

    def energy_D_part_q(self, x):
        n, c = x.unpack()
       #n = p*p
        D = c @ numpy.diag(n) @ c.T

       #F = abs(numpy.outer(n,n))
       #Y = self._U * F**self._q

       #self._U = numpy.identity(n.size)

       #Y = self._U @ numpy.diag(n**self._q) @ self._U.T
       #Y = c @ Y @ c.T

        U = self._U @ c

        Y = U @ numpy.diag(n**self._q) @ U.T

        Ky= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(Y, ""))[1].to_array(dense=True)
       
        energy =      -numpy.dot(Ky, Y).trace()
        return energy

class IDF_JKPlus_XCFunctional(XCFunctional):
    """
"""
    BOUNDS  = [[0.00001,None],]
    TOL = 1.e-7
    MAXITER = 1000
    OPTIONS = {"disp": False, "maxiter": MAXITER, "ftol": TOL, "iprint": 0}

    def __init__(self, omega):
        super(IDF_JKPlus_XCFunctional, self).__init__()
        self._omega = omega

    def load(self):
        "Load miscellanea"
        n = self._wfn.nmo()
        self._DIM_N  = n
        self._DIM_N2 = n*n
        self._DIM_N4 = n*n*n*n
        self._N      = self._wfn.nalpha()*2.
       #self._u      = numpy.linalg.pinv(self._Ca.T @ self._S) # D_ao = u @ D_mo_scf @ u.T

       #L = numpy.dot(self._Ca.T, self._S)
       #R = numpy.dot(numpy.linalg.inv(numpy.dot(L.T, L)), L.T).T
       #self._u = R.T
        self._u = self._Ca.copy()

       #print(self._u-self._Ca)
       #self._u =self._Ca
        r = numpy.ones((n,n)) - numpy.identity(n)
        T = numpy.ones((n,n))
        for i in range(n):
            T[i,i] = 0.0
            for j in range(i): T[i,j] = -1.0
        t = numpy.einsum("ij,kl->ijkl",T,T)
        self._excl_ijkl = 1*1.*numpy.einsum("ij,kl,ik,jl,il,jk->ijkl",r,r,r,r,r,r) #* t
        print("Calculating AO ERIs in memory")
        self._eri_ao = numpy.array( psi4.core.MintsHelper(self._wfn.basisset()).ao_eri() )

        return

    @staticmethod
    def name(): return "IDF JKPlus XC Functional for closed-shell systems"

    @property
    def abbr(self): return "IDF-JKPlus"

    def energy_P(self, x):
        "Exchange-correlation energy"
        p, c = x.unpack(); p = abs(p)

        # JK-Only part
        f = self.fij(p*p)
        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);

        # Plus part
        C = self._u @ c
        eri_NO = numpy.einsum("ijkl,ia,jb,kc,ld->abcd",self._eri_ao,C,C,C,C)

        phases =-numpy.sign(eri_NO.round(9)) #.transpose(0,2,1,3)
        phases = numpy.einsum("ikjl->ijkl", phases)
       #phases = 0.5*(phases + phases.transpose(1,0,2,3)) # symmetrization

        zw=     p**(self._omega+1.0)
        U = numpy.einsum("i,j,k,l->ijkl",zw,zw,zw,zw) - numpy.einsum("i,j,k,l->ijkl",p,p,p,p)
        dg=-0.5*self._excl_ijkl*U
       #dg = self.antisymmetrize_4(dg)
       #dg = 0.5*(dg+numpy.einsum("ijkl->klij",dg))

       # phases =-numpy.sign(eri_NO.transpose(0,2,1,3)*dg)
       # phases = 0.5 * (phases + phases.transpose(1,0,2,3))
       # phases = 0.5 * (phases + phases.transpose(2,3,0,1))

       #dg = phases * dg

       #xc_energy += numpy.einsum("ijkl,ikjl", dg, eri_NO)

        t = dg * numpy.einsum("ikjl->ijkl", eri_NO)
        t = -abs(t) 
        xc_energy += t.sum()
       #print(t.sum())
        return xc_energy

    @staticmethod
    def antisymmetrize_4(t): 
        return 0.25*( ( t - t.transpose(1,0,2,3) ) + (-t.transpose(0,1,3,2) + t.transpose(1,0,3,2)) )
       #return 0.5*(t - t.transpose(1,0,2,3))

    def fij(self, n): 
       #f = n.copy() #/n.sum(); f*= self._N / 2.
        F = self.unfold_2(n, self._omega)

       # # modes
       # mode = 0   # antisym
       ##mode = 1   # antisym_a
       ##mode = 2   # mix

       # if 0:
       #    x_0 = (0.02,)
       #    if mode == 1:
       #       obj = self.obj_antisymmetry_a
       #    elif mode == 0:
       #       obj = self.obj_antisymmetry
       #    else:
       #       obj = self.obj_mix

       #    res = scipy.optimize.minimize(obj, x_0, args=(f,), tol=IDF_JK_XCFunctional.TOL, method='slsqp',
       #                options=IDF_JK_XCFunctional.OPTIONS, bounds=IDF_JK_XCFunctional.BOUNDS, constraints=())

       #    if res.success: 
       #       x = res.x 
       #    else: raise ValueError("Minimization has an error!")
       # else: x = self._omega

       # if 0:
       #    npoints = 1000
       #    a = 0.29
       #    F = numpy.zeros((self._DIM_N, self._DIM_N))
       #    x = numpy.linspace(0.0,10.0,npoints)
       #    w = numpy.exp(-a*x); w/= w.sum()
       #    for i in range(npoints):
       #        F += w[i] * self.unfold_2(n, x[i])
       # else:
       #    if mode == 0 or mode == 2:
       #      #scale = numpy.sqrt(self._DIM_N ) * 0.83
       #      #x/=scale
       #       F = self.unfold_2(n, x)
       #    else:
       #       F = self.unfold_exp(n, x*4)
       ##print(x)
        return F

    def unfold_2(self, f, x, eps=1.e-20):
        "Pseudo 2-rank unfolding function"
        K = f.size
        ff = numpy.outer(f,f) ** ((x+1.)/2.)
        s  = f**(x*(x+1.))
        d  = s[:,numpy.newaxis] + s[numpy.newaxis,:]
        d  = d**(1./(x+1.))
        pref = 2.**(1./(x+1.))
        dbs= abs(d)
    
        F = numpy.zeros(ff.shape)
        for i in range(K):
            F[i,i] = f[i]
            for j in range(i):
               #pref = 1.0
               #if i==j: pref = 2.0**(1./(x+1.))
                if dbs[i,j] > eps:
                   F[i,j] = ff[i,j] / d[i,j]  * pref
                else: F[i,j] = 0.0
                F[j,i] = F[i,j]
       #F *= pref
        return F
 
    def energy_D_add(self, x): #ADD
        "Exchange-correlation energy"
        n, c = x.unpack()
        p    = numpy.sqrt(abs(n))
        P    = c @ numpy.diag(p) @ c.T
      # P = self.generalized_density(n, c, 0.5)
        y    = Guess.create(p, c, P, "matrix")
        xc_energy = self.energy_P(y)
        return xc_energy


class IDF_JK_XCFunctional(XCFunctional):
    """
"""
    BOUNDS  = [[0.00001,None],]
    TOL = 1.e-7
    MAXITER = 1000
    OPTIONS = {"disp": False, "maxiter": MAXITER, "ftol": TOL, "iprint": 0}

    def __init__(self, omega):
        super(IDF_JK_XCFunctional, self).__init__()
        self._omega = omega

    def load(self):
        "Load miscellanea"
        n = self._wfn.nmo()
        self._DIM_N  = n
        self._DIM_N2 = n*n
        self._DIM_N4 = n*n*n*n
        self._N      = self._wfn.nalpha()*2.
        self._u      = numpy.linalg.pinv(self._Ca.T @ self._S) # D_ao = u @ D_mo_scf @ u.T
       #L = numpy.dot(self._Ca.T, self._S)
       #R = numpy.dot(numpy.linalg.inv(numpy.dot(L.T, L)), L.T).T
       #self._u = R.T

       #print(self._u-self._Ca)
       #self._u =self._Ca
        return

    @staticmethod
    def name(): return "IDF JK-Only XC Functional for closed-shell systems"

    @property
    def abbr(self): return "IDF-JK"

    def obj_antisymmetry(self, param, f):
        x = param[0]
        D = numpy.diag(f); I = numpy.identity(len(f))
        F2= self.unfold_2(f, x)
        g_H = numpy.einsum("ik,jl->ijkl",D,D)
        g_X = numpy.einsum("ij,il,jk->ijkl",F2,I,I)
        g_ijkl  = g_H - g_X
        g_ijlk  = g_ijkl.transpose(0,1,3,2)
        error = (g_ijkl + g_ijlk)**2
        return error.sum()

    def obj_antisymmetry_and_positivity(self, param, f):
        x = param[0]
        D = numpy.diag(f); I = numpy.identity(len(f))
        F2= self.unfold_2(f, x)
        g_H = numpy.einsum("ik,jl->ijkl",D,D)
        g_X = numpy.einsum("ij,il,jk->ijkl",F2,I,I)
        # antisymmetry
        g_ijkl  = g_H - g_X
        g_ijlk  = g_ijkl.transpose(0,1,3,2)
        error_antisymmetry  = ( (g_ijkl + g_ijlk)**2 ).sum()
        # positivity
        E, U = numpy.linalg.eigh(g_ijkl.reshape(self._DIM_N2,self._DIM_N2))
        error_positivity = -E.min()
        return error_antisymmetry, error_positivity


    def unfold_exp(self, f, a):
        npoints = 100
        F2= numpy.zeros((self._DIM_N, self._DIM_N))
        x = numpy.linspace(0.0,10.0,npoints)
        w = numpy.exp(-a*x*x); w/= w.sum()
        for i in range(npoints):
            F2+= w[i] * self.unfold_2(f, x[i])
        return F2

    def obj_antisymmetry_a(self, param, f):
        x = param[0]
        D = numpy.diag(f); I = numpy.identity(len(f))
        F2 = self.unfold_exp(f, x)
        g_H = numpy.einsum("ik,jl->ijkl",D,D)
        g_X = numpy.einsum("ij,il,jk->ijkl",F2,I,I)
        g_ijkl  = g_H - g_X
        g_ijlk  = g_ijkl.transpose(0,1,3,2)
        error = (g_ijkl + g_ijlk)**2
        return error.sum()

    def obj_mix(self, param, f):
        w_max = [6.0]
        Z_1    , Z_2     = self.obj_antisymmetry_and_positivity(param, f)
        Z_1_max, Z_2_max = self.obj_antisymmetry_and_positivity(w_max, f)
       #print(Z_1_max, Z_2_max)

        z_1 = Z_1/Z_1_max
        z_2 = Z_2/Z_2_max

        error = z_1 + z_2
        return error


    def fij(self, n): 
        f = n.copy() #/n.sum(); f*= self._N / 2.

        # modes
        mode = 0   # antisym
       #mode = 1   # antisym_a
       #mode = 2   # mix

        if 0:
           x_0 = (0.02,)
           if mode == 1:
              obj = self.obj_antisymmetry_a
           elif mode == 0:
              obj = self.obj_antisymmetry
           else:
              obj = self.obj_mix

           res = scipy.optimize.minimize(obj, x_0, args=(f,), tol=IDF_JK_XCFunctional.TOL, method='slsqp',
                       options=IDF_JK_XCFunctional.OPTIONS, bounds=IDF_JK_XCFunctional.BOUNDS, constraints=())

           if res.success: 
              x = res.x 
           else: raise ValueError("Minimization has an error!")
        else: x = self._omega

        if 0:
           npoints = 1000
           a = 0.29
           F = numpy.zeros((self._DIM_N, self._DIM_N))
           x = numpy.linspace(0.0,10.0,npoints)
           w = numpy.exp(-a*x); w/= w.sum()
           for i in range(npoints):
               F += w[i] * self.unfold_2(n, x[i])
        else:
           if mode == 0 or mode == 2:
             #scale = numpy.sqrt(self._DIM_N ) * 0.83
             #x/=scale
              F = self.unfold_2(n, x)
           else:
              F = self.unfold_exp(n, x*4)
       #print(x)
        return F

    def unfold_2(self, f, x, eps=1.e-20):
        "Pseudo 2-rank unfolding function"
        K = f.size
        ff = numpy.outer(f,f) ** ((x+1.)/2.)
        s  = f**(x*(x+1.))
        d  = s[:,numpy.newaxis] + s[numpy.newaxis,:]
        d  = d**(1./(x+1.))
        pref = 2.**(1./(x+1.))
        dbs= abs(d)
    
        F = numpy.zeros(ff.shape)
        for i in range(K):
            F[i,i] = f[i]
            for j in range(i):
               #pref = 1.0
               #if i==j: pref = 2.0**(1./(x+1.))
                if dbs[i,j] > eps:
                   F[i,j] = ff[i,j] / d[i,j]  * pref
                else: F[i,j] = 0.0
                F[j,i] = F[i,j]
       #F *= pref
        return F
 
    def energy_D_add(self, x): #ADD
        "Exchange-correlation energy"
        n, c = x.unpack()
        p    = numpy.sqrt(abs(n))
        P    = c @ numpy.diag(p) @ c.T
      # P = self.generalized_density(n, c, 0.5)
        y    = Guess.create(p, c, P, "matrix")
        xc_energy = self.energy_P(y)
        return xc_energy



class IDF_JK_Di_XCFunctional(XCFunctional):
    """
 Discrete
"""
    def __init__(self, ngrid, constant):
        super(IDF_JK_Di_XCFunctional, self).__init__()
        self._ngrid = ngrid
        self._omega_domain = numpy.linspace(0.0, 2.*numpy.pi, ngrid)
        self._weights = None
        self._constant_weights = constant
        if constant:
           data = numpy.mafromtxt('weights.dat').data
           self._omega_domain = data[:,0]
           self._weights = abs(data[:,1]); self._weights/=self._weights.sum()
           self._ngrid = len(data)
        self._Q = self._q(self._omega_domain)

    def _p(self, omega): return (3.+numpy.cos(omega/2.))/4.
    def _q(self, omega): return (3.-numpy.cos(omega/2.))/4.

    def load(self):
        "Load miscellanea"
        n = self._wfn.nmo()
        self._DIM_N  = n
        self._DIM_N2 = n*n
        self._DIM_N4 = n*n*n*n
        self._DIM_W  = self._ngrid
        self._N      = self._wfn.nalpha()*2.
        self._u      = numpy.linalg.pinv(self._Ca.T @ self._S) # D_ao = u @ D_mo_scf @ u.T
       #L = numpy.dot(self._Ca.T, self._S)
       #R = numpy.dot(numpy.linalg.inv(numpy.dot(L.T, L)), L.T).T
       #self._u = R.T

       #print(self._u-self._Ca)
       #self._u =self._Ca
        return

    @staticmethod
    def name(): return "IDF JK-Only Discrete XC Functional for closed-shell systems"

    @property
    def abbr(self): return "IDF-JK-Di"

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

    def weights(self, n, c):
        "Discretized spectral density - normalized to 1.0"
        n_ = n/n.sum(); n_*= self._N/2.
        D = c @ numpy.diag(n_) @ c.T

        f = numpy.zeros((self._DIM_W, self._DIM_N, self._DIM_N))
        for i in range(self._DIM_W):
            f[i] = c @ numpy.diag(n_**self._Q[i]) @ c.T
            
        F = numpy.einsum("aij,ajk->aik",f,f)                                                       
       #F = f
                                                                                             
        # Dimensions
        K = self._DIM_N
        W = self._DIM_W

        l = numpy.zeros(K*K*K*K + 1); l[-1] = 1.0
                                                                                             
        l[:  -1] = ( (1./2.)*(numpy.einsum("ik,jl->ijkl",D,D)+numpy.einsum("il,jk->ijkl",D,D)) ).ravel()
        H        = numpy.einsum("aik,ajl->aijkl",f, f)+numpy.einsum("ail,ajk->aijkl",f, f) \
                -(1./4.)*(numpy.einsum("aik,jl->aijkl", F, D) + numpy.einsum("ik,ajl->aijkl", D, F) + \
                          numpy.einsum("ail,jk->aijkl", F, D) + numpy.einsum("il,ajk->aijkl", D, F) )
                                                                                             
        h = numpy.zeros((W, K*K*K*K+1))
        for i in range(W):
            h[i,:-1] = H[i].ravel()
        h[:,-1] = 1.0
       #print(linalg.det(h))
                                                                                             
        # Compute w_i
        hi = numpy.linalg.pinv(h)
        w  = l @ hi

        self._weights = w
        if 1:
           out = open('weights.dat', 'w')                      
           for i in range(self._ngrid):
               out.write("%14.5f %13.5f\n" % (self._omega_domain[i],self._weights[i]))
           out.close()
          #print("Sum = ", self._weights.sum())
          #exit()

        return w

    def weights_old(self, n, c):
        "Discretized spectral density - normalized to 1.0"
        n_ = n/n.sum(); n_*= self._N/2.
        D = c @ numpy.diag(n_) @ c.T

        f = numpy.zeros((self._DIM_W, self._DIM_N, self._DIM_N))
        for i in range(self._DIM_W):
            f[i] = c @ numpy.diag(n_**self._Q[i]) @ c.T
            
        F = numpy.einsum("aij,ajk->aik",f,f)                                                       

        # Dimensions
        K = self._DIM_N
        W = self._DIM_W

        l = numpy.zeros(K*K + 1); l[-1] = 1.0
                                                                                             
        l[:  -1] = D.ravel()

        h = numpy.zeros((W, K*K+1))
        for i in range(W):
            h[i,:-1] = F[i].ravel()
        h[:,-1] = 1.0
       #print(linalg.det(h))
                                                                                             
        # Compute w_i
        hi = numpy.linalg.pinv(h)
        w  = l @ hi

        self._weights = w
        if 1:
           out = open('weights.dat', 'w')                      
           for i in range(self._ngrid):
               out.write("%14.5f %13.5f\n" % (self._omega_domain[i],self._weights[i]))
           out.close()
          #print("Sum = ", self._weights.sum())
          #exit()

        return w
   
    def energy_P(self, x):
        "Exchange-correlation energy"
        p, c = x.unpack()
        n = p*p;
       #D = c @ numpy.diag(n) @ c.T
       #P = x.matrix()

        # weights
        if not self._constant_weights:
           weights = self.weights(n, c)

        # XC energy
        f = numpy.zeros((len(n),len(n)))
        for i in range(self._ngrid):
            m = n**self._Q[i]
            f+= numpy.outer(m,m) * self._weights[i]

        psi_f = psi4.core.Matrix.from_array(f, "")
        psi_c = psi4.core.Matrix.from_array(c, "")
        xc_energy = oepdev.calculate_e_xc(self._wfn, self._ints, psi_f, psi_c);

        return xc_energy


class IDF_Alpha_XCFunctional(XCFunctional):
    """
 The IDF test version of functional for particular alpha value.

 Alpha value: Unitary transformation matrix for 4-rank tensor unfolding
              Alpha = 0 -> JK-Only functional
              Alpha must be larger than 0 and larger than the critical value
              (around 1.4) in order to satisfy the N-representability contition
              exactly.

 Interpolation: 1-D domain of omega from 0.0 to 2PI.
                Number of points for interpolation is controlled.
"""
    def __init__(self, alpha, ngrid):
        super(IDF_Alpha_XCFunctional, self).__init__()
        self._alpha = alpha
        self._ngrid = ngrid
        self._ngrid_half = int(ngrid/2)
        self._omega_domain = numpy.linspace(0.0, 2.*numpy.pi, ngrid)
        self._weights = None

    def _p(self, omega): return (3.+numpy.cos(omega/2.))/4.
    def _q(self, omega): return (3.-numpy.cos(omega/2.))/4.
    def _U(self, omega):
        R = self._wfn.nmo()
        a = self._alpha ** R
        na= self._wfn.nalpha()
        A = numpy.zeros((R,R))
        s = numpy.sin(omega/2.)
        v = a * s*s
        for i in range(R):
         #if i >= na:
            for j in range(i):
             #if ((i <na) and (j <na)):
                A[i,j] = v
                A[j,i] =-v
        U = scipy.linalg.expm(A)
       #print(U); exit()
        return U


    @staticmethod
    def name(): return "IDF-Alpha Test Version of XC Functional"

    @property
    def abbr(self): return "IDF-Alpha"

    @staticmethod
    def fij(n): 
        raise NotImplementedError
 
    def load(self):
        "Load miscellanea"
        n = self._wfn.nmo()
        self._DIM_N2 = n*n
        self._DIM_W  = self._ngrid_half
        self._N      = self._wfn.nalpha()*2.
        self._u      = numpy.linalg.pinv(self._Ca.T @ self._S) # D_ao = u @ D_mo_scf @ u.T
       #L = numpy.dot(self._Ca.T, self._S)
       #R = numpy.dot(numpy.linalg.inv(numpy.dot(L.T, L)), L.T).T
       #self._u = R.T

       #print(self._u-self._Ca)
       #self._u =self._Ca
        return
   
    def gradient_P(self, x):
        "Analytical gradient"
        p, c = x.unpack()
        P = x.matrix()
        P2= P @ P
        raise NotImplementedError

    def weights_P(self, x):
        "Discretized Spectral Density over Interpolating Domain"
        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T

        # allocate
        F = numpy.ones ((self._DIM_W, self._DIM_N2+1))
        g = numpy.zeros(              self._DIM_N2+1 )
       #print(F.shape); exit()

        # populate g vector
        d = (self._N - 1.0) * n
        g[:-1] = numpy.diag(d).ravel()
       #g[:-1] = D.ravel() * (self._N -1.0)
        g[ -1] = 1./2.

        # populate F matrix
        Al = []; Bl = []
        for i in range(self._ngrid_half):
            omega = self._omega_domain[i]
            p     = self._p(omega)
            q     = self._q(omega)

            U = self._U(omega) @ c

            A = U @ numpy.diag(n**p) @ U.T
            B = U @ numpy.diag(n**q) @ U.T

            Al.append(A.copy())
            Bl.append(B.copy())

            A = c.T @ A @ c
            B = c.T @ B @ c

            f = 2.0 * (A * A.trace() + B * B.trace()) - (A@A + B@B)

            F[i,:-1] = f.ravel()

        # compute weights
        Fi = numpy.linalg.pinv(F)
        w  = g @ Fi

        # glue with the mirror weights
        wall = numpy.concatenate((w, w[::-1]))

       #wall.fill(0.0); wall[0] = wall[-1] = 0.5
       #wall.fill(0.0); wall[0] = 1.0            
        return wall, Al, Bl

    def energy_P(self, x):
        "Exchange-correlation energy"
        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T

        weights, Al, Bl = self.weights_P(x)
        if 0:
           out = open('weights.dat', 'w')                      
           for i in range(self._ngrid):
               out.write("%14.3d %13.5f\n" % (i+1,weights[i]))
           out.close()
           exit()

        # transform to non-orthogonal AO basis (from orthogonal SCF-MO basis)
       #self._u = numpy.identity(len(n))
        D_ao = self._u @ D @ self._u.T
        for i in range(self._ngrid_half):
            A_ao = self._u @ Al[i] @ self._u.T
            B_ao = self._u @ Bl[i] @ self._u.T
           #if i==0:
           #   print(A_ao);
            Al[i] = psi4.core.Matrix.from_array(A_ao, "")
            Bl[i] = psi4.core.Matrix.from_array(B_ao, "")

        D_ao_psi = psi4.core.Matrix.from_array(D_ao, "")

        xc_energy = oepdev.calculate_idf_alpha_xc_energy(self._wfn, weights[:self._ngrid_half], D_ao_psi, Al, Bl)
        return xc_energy

    def potential_P(self, x):
        "Exchange-Correlation Potential in MO basis"
        p, c = x.unpack()
        n = p*p
        D = c @ numpy.diag(n) @ c.T

        weights, Al, Bl = self.weights_P(x)

        Jd= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(D, ""))[0].to_array(dense=True)
        V =-Jd * 2.

        for i in range(self._ngrid):
            omega = self._omega_domain[i]
            p     = self._p(omega)
            q     = self._q(omega)
            Uw    = self._U(omega)

            weight= weights[i]

            pref_p = 1./(2.*p)
            pref_q = 1./(2.*q)

            y = Guess.create(n, c, D, "matrix")
            V+=  (self._gradient_D_numerical_unfolding_part_p(y, p, omega, Uw).matrix()) * pref_p * weight
            V+=  (self._gradient_D_numerical_unfolding_part_q(y, q, omega, Uw).matrix()) * pref_q * weight

        self._weights = weights
        out = open('weights.dat', 'w')
        for i in range(self._ngrid):
            out.write("%14.3d %13.5f\n" % (i+1,weights[i]))
        out.close()

        return Guess.create(matrix=V)  

    def _gradient_D_numerical_unfolding_part_p(self, x, p, omega, Uw):
        "Numerical gradient with respect to P matrix"
        E0 = self.energy_D_part_p(x, p, omega, Uw)
        D  = x.matrix()
        N  = D.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_D_part_p(i, j, D, E0, p, omega, Uw)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)

    def _gradient_D_numerical_unfolding_part_q(self, x, q, omega, Uw):
        "Numerical gradient with respect to P matrix"
        E0 = self.energy_D_part_q(x, q, omega, Uw)
        D  = x.matrix()
        N  = D.shape[0]
        Der = numpy.zeros((N,N), numpy.float64)
        for i in range(N):
            for j in range(i+1):
                d = self._compute_deriv_D_part_q(i, j, D, E0, q, omega, Uw)
                Der[i, j] = d
                if i!=j: Der[j, i] = d
        return Guess.create(matrix=Der)

    def _compute_deriv_D_part_p(self, n, m, D, E0, p, omega, Uw): #ADD
        step = self._ffstep; h = 2.0*step
        D1= D.copy()
        D1[m,n] += step; D1[n,m] += step

        guess = Guess.create(matrix=D1)
        guess.update()
        E1 = self.energy_D_part_p(guess, p, omega, Uw)
        der = (E1 - E0) / h
        return der 

    def _compute_deriv_D_part_q(self, n, m, D, E0, q, omega, Uw): #ADD
        step = self._ffstep; h = 2.0*step
        D1= D.copy()
        D1[m,n] += step; D1[n,m] += step

        guess = Guess.create(matrix=D1)
        guess.update()
        E1 = self.energy_D_part_q(guess, q, omega, Uw)
        der = (E1 - E0) / h
        return der 

    def energy_D_part_p(self, x, p, omega, Uw):
        n, c = x.unpack()
        U = Uw @ c

        A = U @ numpy.diag(n**p) @ U.T

        Ja= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(A, ""))[0].to_array(dense=True)
       
        energy = 2.0 * numpy.dot(Ja, A).trace()
        return energy

    def energy_D_part_q(self, x, q, omega, Uw):
        n, c = x.unpack()
        U = Uw @ c

        B = U @ numpy.diag(n**q) @ U.T

        Kb= oepdev.calculate_JK_r(self._wfn, self._ints, psi4.core.Matrix.from_array(B, ""))[1].to_array(dense=True)
       
        energy =      -numpy.dot(Kb, B).trace()
        return energy






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

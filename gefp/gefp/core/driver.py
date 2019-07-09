#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 GEFP Driver Routines.
 Bartosz Błasiak, Gundelfingen, 7 May 2019
"""

import os
import sys
import math
import numpy
import psi4
from ..density.dmft import DMFT
from ..density.functional import XCFunctional
from ..basis.optimize import *

__all__ = ["dmft_solver", "gdf_basisset_fitter"]

class PadeApproximant_2D:
    """\
 Pade Approximant of Function of 2 Variables:

 f(x, y) =        \sum_{n=0,m=0}            a_{nm} x^n y^m   \ 
            1.0 + \sum_{n,m}^{n,m != (0,0)} b_{nm} x^n y^m
"""
    def __init__(self):
        "Initialize"
        self.a = []
        self.b = []
    def add_a(self, nx, ny, a):
        "Add a_nm coefficient to the Pade approximant"
        self.a.append((nx, ny, a))
    def add_b(self, nx, ny, b):
        "Add b_nm coefficient to the Pade approximant"
        assert(nx + ny != 0)
        self.b.append((nx, ny, b))
    def value(self, x, y):
        "Compute value of the function f(x,y) from Pade Approximant"
        f = 0.0
        d1= 0.0
        d2= 1.0
        for key in self.a:
            nx, ny, a = key
            d1 += a * x**nx * y**ny
        for key in self.b:
            nx, ny, b = key
            d2 += b * x**nx * y**ny
        f = d1/d2
        return f

class Entry:
    def __init__(self, pade, description_short, description_full):
        self.pade = pade
        self.description_short= description_short
        self.description_full = description_full

class UniversalSurface:
    #
    par_descr_fci_sto3g_1 = """\
 Universal Surface at FCI/STO-3G level.
 Systems: 2el (H2), 4el (H4), 10el (H2O)
"""
    par_fulld_fci_sto3g_1 = """\
 function used for fitting: g(x,y)
	g(x,y) = (a0 + a1*x + a2*y + a3*x*y + a4*y*y + a5*y*y*x + a6*y*y*y)/(1.0 + b1*x + b2*y + b3*x*y + b4*y*y + b5*y*y*x + b6*y*y*y)

 degrees of freedom    (FIT_NDF)                        : 47
 rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.0217961
 variance of residuals (reduced chisquare) = WSSR/ndf   : 0.00047507
 
 Final set of parameters            Asymptotic Standard Error
 =======================            ==========================
 a0              = 0.173323         +/- 0.007845     (4.526%)
 a1              = -6.85743         +/- 3.626        (52.87%)
 a2              = 0.732897         +/- 0.5255       (71.7%)
 a3              = -8.10918         +/- 3.951        (48.73%)
 a4              = -0.272454        +/- 0.5683       (208.6%)
 a5              = 1.73152          +/- 3.307        (191%)
 a6              = -0.0591462       +/- 0.08416      (142.3%)
 b1              = -68.04           +/- 38.64        (56.8%)
 b2              = -3.57443         +/- 3.66         (102.4%)
 b3              = 23.0961          +/- 43.09        (186.5%)
 b4              = -3.80028         +/- 3.203        (84.29%)
 b5              = 30.2229          +/- 27.85        (92.14%)
 b6              = -0.269245        +/- 0.328        (121.8%)

"""
    par_pade_fci_sto3g_1 = PadeApproximant_2D()
    par_pade_fci_sto3g_1.add_a(0, 0,   0.173323)  # a0
    par_pade_fci_sto3g_1.add_a(1, 0,  -6.85743 )  # a1
    par_pade_fci_sto3g_1.add_a(0, 1,   0.732897)  # a2
    par_pade_fci_sto3g_1.add_a(1, 1,  -8.10918 )  # a3
    par_pade_fci_sto3g_1.add_a(0, 2,  -0.272454)  # a4
    par_pade_fci_sto3g_1.add_a(1, 2,   1.73152 )  # a5
    par_pade_fci_sto3g_1.add_a(0, 3, -0.0591462)  # a6
    #
    par_pade_fci_sto3g_1.add_b(1, 0,-68.04    )   # b1
    par_pade_fci_sto3g_1.add_b(0, 1, -3.57443 )   # b2
    par_pade_fci_sto3g_1.add_b(1, 1, 23.0961  )   # b3
    par_pade_fci_sto3g_1.add_b(0, 2, -3.80028 )   # b4
    par_pade_fci_sto3g_1.add_b(1, 2, 30.2229  )   # b5
    par_pade_fci_sto3g_1.add_b(0, 3, -0.269245)   # b6


    #
    par_descr_fci_sto3g_2 = """\
 Universal Surface at FCI/STO-3G level.
 Systems: 2el (H2), 4el (H4), 10el (H2O)
"""
    par_fulld_fci_sto3g_2 = """\

function used for fitting: g(x,y)
        g(x,y) = (a0 + a1*x + a2*y + a3*x*y + a4*y*y + a5*y*y*y + a6*y*y*y*y + a7*y**5 + a8*x*x*y)/(1.0 + b1*x + b2*y + b3*x*y + b4*y*y + b5*y*x*x + b6*x*x + b7*x*y*y + b8*x*x*y*y + b9*x*x*x*y)

degrees of freedom    (FIT_NDF)                        : 42
rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.0131158
variance of residuals (reduced chisquare) = WSSR/ndf   : 0.000172024


Final set of parameters            Asymptotic Standard Error
=======================            ==========================
a0              = 0.180128         +/- 0.004829     (2.681%)
a1              = -7.77423         +/- 1.696        (21.81%)
a2              = 1.87765          +/- 0.647        (34.46%)
a3              = -6.73461         +/- 2.739        (40.67%)
a4              = 2.75464          +/- 1.005        (36.47%)
a5              = 1.60223          +/- 0.5967       (37.24%)
a6              = 0.372731         +/- 0.1429       (38.33%)
a7              = 0.0299822        +/- 0.01182      (39.44%)
b1              = -68.707          +/- 14.46        (21.04%)
b2              = 1.35046          +/- 0.3228       (23.9%)
b3              = -112.231         +/- 19.2         (17.1%)
b4              = 0.342321         +/- 0.2692       (78.64%)
a8              = 12.9053          +/- 14.26        (110.5%)
b5              = -890.429         +/- 482.9        (54.23%)
b6              = -1109.35         +/- 351.1        (31.65%)
b7              = -45.5445         +/- 8.899        (19.54%)
b8              = -32.6005         +/- 159.2        (488.4%)
b9              = 1703.16          +/- 1828         (107.4%)
"""
    par_pade_fci_sto3g_2 = PadeApproximant_2D()
    par_pade_fci_sto3g_2.add_a(0, 0, 0.180128)
    par_pade_fci_sto3g_2.add_a(1, 0, -7.77423)
    par_pade_fci_sto3g_2.add_a(0, 1, 1.87765)
    par_pade_fci_sto3g_2.add_a(1, 1, -6.73461)
    par_pade_fci_sto3g_2.add_a(0, 2, 2.75464)
    par_pade_fci_sto3g_2.add_a(0, 3, 1.60223)
    par_pade_fci_sto3g_2.add_a(0, 4, 0.372731)
    par_pade_fci_sto3g_2.add_a(0, 5, 0.0299822)
    par_pade_fci_sto3g_2.add_a(2, 1, 12.9053)
    #
    par_pade_fci_sto3g_2.add_b(1, 0, -68.707)
    par_pade_fci_sto3g_2.add_b(0, 1, 1.35046)
    par_pade_fci_sto3g_2.add_b(1, 1, -112.231)
    par_pade_fci_sto3g_2.add_b(0, 2, 0.342321)
    par_pade_fci_sto3g_2.add_b(2, 1, -890.429)
    par_pade_fci_sto3g_2.add_b(2, 0, -1109.35)
    par_pade_fci_sto3g_2.add_b(1, 2, -45.5445)
    par_pade_fci_sto3g_2.add_b(2, 2, -32.6005)
    par_pade_fci_sto3g_2.add_b(3, 1, 1703.16)
   
    # parameterizations
    par_fci_sto3g_1 = Entry(par_pade_fci_sto3g_1, par_descr_fci_sto3g_1, par_fulld_fci_sto3g_1)
    par_fci_sto3g_2 = Entry(par_pade_fci_sto3g_2, par_descr_fci_sto3g_2, par_fulld_fci_sto3g_2)


    
def dmft_solver(wfn, xc_functional='MBB', 
                g_0=0.001, g=0.0001, verbose=True, guess='current', step_mode='search', algorithm='proj-p',
                gradient_mode='exact', conv=0.000005, conv_int=0.000010, maxit=300,
                with_restart=False, parameterization=None, t=None, kmax=30,
                return_density=False):
    "Compute the electron density by using the DMFT method."

    do_medi = True if 'medi' in xc_functional.lower() else False

    # Interpolation XC Functional
    if do_medi:
       if algorithm.lower() != 'proj-p':
          print(" Warning: For MEDI functionals algorithm must be 'PROJ-P'. Ignoring current option 'algorithm'.")
          algorithm = 'proj-p'

       # compute t
       if t is None:
          surface  = UniversalSurface()                                                         
          func_mbb = XCFunctional.create('MBB')
          dmft_mbb = DMFT.create(wfn, func_mbb, 
                                 algorithm='proj-p', 
                                 guess='current', 
                                 step_mode='search')
          dmft_mbb.run(maxit=maxit, conv=conv, g_0=g_0, verbose=verbose)
                                                                                               
          N          = int(wfn.nalpha())
          E_mbb      = dmft_mbb.E
          D_mbb      = dmft_mbb.D
          d_mbb      = dmft_mbb.dipole
          i_d, i_n   = dmft_mbb.scalar_correlation
          #I_d, I_n   = dmft_mbb.matrix_correlation
          #I_d        = numpy.dot(I_d, D_mbb).trace() * 2.0
          #I_n        = numpy.dot(I_n, D_mbb).trace() * 2.0
          x          = math.log(i_d/N + 1.0)
          y          = math.log(i_n/N * 2.0) if i_n > 0.0 else -1.e10
                                                                                               
          if   parameterization=='fci/sto-3g.1':
               par = surface.par_fci_sto3g_1
          elif parameterization=='fci/sto-3g.2':
               par = surface.par_fci_sto3g_2
          else: raise ValueError("Incorret Universal Surface Specification!")
          t = par.pade.value(x, y)
          if verbose: 
             print(par.description_short)
             print(" MBB Electron Correlation: i_d=%14.4f in=%14.4f t=%14.4f" % (i_d, i_n, t))
          del func_mbb, dmft_mbb
          psi4.core.clean()

       else: pass
                                                                                            
       # Run DMFT with interpolation functional
       if '1' in xc_functional.lower(): coeff = {'a0': t}
       else:                            coeff = {'t' : t}
       func_int = XCFunctional.create(xc_functional, coeff=coeff, kmax=kmax)

    # Conventional XC Functionals
    else:
       func_int = XCFunctional.create(xc_functional)

    # Compute density
    dmft_int = DMFT.create(wfn, func_int, 
                           algorithm=algorithm, 
                           guess=guess, 
                           step_mode=step_mode)

    # All XC Functionals
    if   gradient_mode.lower() == 'num':
         dmft_int.set_gradient_mode(num=True)
         dmft_int.run(maxit=maxit, conv=conv    , g_0=g_0, g=g, verbose=verbose)
    # Only MEDI XC Functionals
    elif gradient_mode.lower() == 'exact' and do_medi: 
         dmft_int.set_gradient_mode(approx=True)                   
         dmft_int.run(maxit=maxit, conv=conv_int, g_0=g_0, g=g, verbose=verbose)
         dmft_int.set_gradient_mode(exact=True)
         dmft_int.run(maxit=maxit, conv=conv    , g_0=g_0, g=g, verbose=verbose, restart=True)
    # Only non-MEDI XC Functionals
    elif gradient_mode.lower() == 'exact' and not do_medi:
         if with_restart:
             dmft_int.set_gradient_mode(approx=True)
             dmft_int.run(maxit=maxit, conv=conv_int, g_0=g_0, g=g, verbose=verbose)
             dmft_int.set_gradient_mode(exact=True)
             dmft_int.run(maxit=maxit, conv=conv    , g_0=g_0, g=g, verbose=verbose, restart=True)
         else:
             dmft_int.run(maxit=maxit, conv=conv    , g_0=g_0, g=g, verbose=verbose)
    # Only non-MEDI XC Functionals
    elif gradient_mode.lower() == 'approx' and not do_medi:
         dmft_int.set_gradient_mode(approx=True)
         dmft_int.run(maxit=maxit, conv=conv    , g_0=g_0, g=g, verbose=verbose)
    #
    else: pass
        
    # 
    E_int = dmft_int.E
    D_int = dmft_int.D
    N_int = dmft_int.N
    C_int = dmft_int.C

    # Return
    if return_density: 
        return E_int, func_int, dmft_int
    else: 
        del func_int, dmft_int
        psi4.core.clean()
        return E_int


def gdf_basisset_fitter(mol, oep_type, basis="6-31G",
                        basis_xpl="6-311++G**", basis_int="AUG-CC-PVDZ-JKFIT"):
    "Fit the auxiliary basis set for GDF-OEP purposes"

    # extract the intermediate basis set
    basis_gdf_int = psi4.core.BasisSet.build(mol, "BASIS", basis_int, fitrole='ORBITAL', puream=-1)
    basis_gdf_xpl = psi4.core.BasisSet.build(mol, "BASIS", basis_xpl, fitrole='ORBITAL', puream=-1)
    
    # prepare the system
    e_hf, w_hf = psi4.energy('hf/%s' % basis, molecule=mol, return_wfn=True)
    w_hf.set_basisset("BASIS_DF_INT", basis_gdf_int)
    
    # fit the auxiliary basis
    dfbasis = DFBasis(w_hf.molecule())
    oep     = OEP.create(oep_type, w_hf, dfbasis)
    opt     = DFBasisOptimizer(oep)
    
    success = opt.fit()
    
    # fetch the optimal basis
    basis_gdf_aux = opt.oep.dfbasis.basis
    
    # compute error for intermediate basis
    Z_aux = opt.compute_error(basis_gdf_aux)
    Z_xpl = opt.compute_error(basis_gdf_xpl)
    Z_int = opt.compute_error(basis_gdf_int)
    print(" Error for auxiliary    basis = %14.6f  Size: %6d" % (Z_aux, basis_gdf_aux.nbf()))
    print(" Error for example      basis = %14.6f  Size: %6d" % (Z_xpl, basis_gdf_xpl.nbf()))
    print(" Error for intermediate basis = %14.6f  Size: %6d" % (Z_int, basis_gdf_int.nbf()))
    return

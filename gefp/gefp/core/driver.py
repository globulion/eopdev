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
from ..basis.optimize import OEP, DFBasis, DFBasisOptimizer, oep_ao_basis_set_optimizer
from ..basis.parameters import StandardizedInput

__all__ = ["gdf_basisset_optimizer"]


def gdf_basisset_optimizer(mol, oep_type, 
                           out='oepfit.gbs', maxiter=1000, tolerance=1.0e-9, method='slsqp',
                           opt_global=False, temperature=500, stepsize=500,
                           take_step=None, accept_test=None,
                           basis=None, 
                           basis_int="AUG-CC-PVDZ-JKFIT", basis_xpl="6-311++G**",
                           templ_file='templ.dat', param_file='param.dat', 
                           bounds_file=None, constraints=(),
                           exp_lower_bound=None, exp_upper_bound=None, 
                           ctr_lower_bound=None, ctr_upper_bound=None,
                           old_interface=False, use_standardized_input='standard'):
    """
 ---------------------------------------------------------------------------------------------
 Fit the auxiliary basis set for GDF-OEP purposes.

 ---------------------------------------------------------------------------------------------

 Usage:
  
   basis_aux = gdf_basisset_optimizer(molecule, oep_type, 
                                      out             = 'oepfit.gbs',
                                      basis           =  None,               
                                      basis_xpl       = "6-311++G**", 
                                      basis_int       = "AUG-CC-PVDZ-JKFIT",
                                      templ_file      = "templ.dat",
                                      param_file      = "param.dat",
                                      bounds_file     = None,
                                      constraints     = (),
                                      maxiter         = 1000,
                                      exp_lower_bound = None,
                                      exp_upper_bound = None,
                                      ctr_lower_bound = None,
                                      ctr_upper_bound = None)


   where

    o basis_aux       - the optimalized auxiliary basis set object of type DFBasis                                 
    o molecule        - psi4.core.Molecule object
    o oep_type        - type of OEP overlap-like matrix element. Available: `pauli`, `ct`.
    o out             - output file to save the optimized basis set in Psi4 format. No output if out=None.
    o basis           - primary basis set of calculation. If not given, the global basis is used (keyword BASIS).
    o basis_int       - intermediate basis set for double GDF calculations (should be RI or JKFIT)
    o basis_xpl       - exemplary basis set to compare (any basis set)
    o templ_file      - template for auxiliary basis set structure
    o param_file      - list of parameters to optimize
    o bounds_file     - list of parameter bounds. If not given, exponents only are assumed 
                        with bounds (exp_lower_bound, exp_upper_bound)
    o constraints     - parameter constraints. If provided, need to be in format 
                        of scipy.optimize.minimize constraints
    o maxiter         - maximum number of optimization iterations
    o exp_lower_bound - custom value of default exponent lower bound 
    o exp_upper_bound - custom value of default exponent upper bound 
    o ctr_lower_bound - custom value of default contraction coefficient lower bound 
    o ctr_upper_bound - custom value of default contraction coefficient lower bound 


   Notes:

    o at least two input files are necessary:
      - 'templ.dat' - template of basis set
      - 'param.dat' - file with parameters to optimize
      - 'bounds.dat'- optional file with parameter bounds. Each bound is specified 
                      by a symbol `E` for exponents and `C` for contraction coefficients.
                      If only symbol is provided, then defaults are taken. If you want to provide
                      custom values, use the syntax `X:min,max` without spaces.

    o example of template file and parameter file contents:
 
      Example for basis for H atom: (2s, 1p) -> [1s, 1p] contraction scheme
      There are 5 parameters to optimize.

      >>>
      cartesian
      ****
      H     0
      S   2   1.00
      %20.10f             %20.10f
      %20.10f             %20.10f
      P   1   1.00
      %20.10f             1.0000000
      <<<

      Parameter file then has to have these entries, e.g.,

      >>>
      # maybe some title can be placed or other comments
      20.0  0.3  # H 1s orbital: 1st GTO (starting values)
      10.3  0.7  # H 1s orbital: 2nd GTO (starting values)
      # some other comments can be provided
      2.4        # H 2p orbital: only exponent is to specify (starting value)
      <<< 

      If you want to specify your own bounds, provide an additional input file:

      >>>
      # maybe some title can be placed or other comments
      E             C:-0.4,0.5  # H 1s orbital: 1st GTO (starting values)
      E:0.03,400.0  C           # H 1s orbital: 2nd GTO (starting values)
      # some other comments can be provided
      E                         # H 2p orbital: only exponent is to specify (starting value)
      <<< 
     
      In the above example, first and third exponent and second contraction coefficient bounds 
      are taken as defaults whereas the rest is provided with custom bounds.

    o Note that you can provide comments after '#' sign. 
    
    o Remember to place the parameters in exact same order as present in template
      from left to right, up to down direction. 
  ---------------------------------------------------------------------------------------------
                                                       Created      : Dresden     ,  9 Jul 2019
                                                       Last Revision: Gundelfingen, 12 Jul 2019
"""

    # extract the intermediate basis set
    puream = psi4.core.get_global_option("PUREAM")
    basis_gdf_int = psi4.core.BasisSet.build(mol, "BASIS", basis_int, fitrole='ORBITAL', puream=puream, quiet=True)
    basis_gdf_xpl = psi4.core.BasisSet.build(mol, "BASIS", basis_xpl, fitrole='ORBITAL', puream=puream, quiet=True)
    
    # prepare the system
    if basis is None: basis = psi4.core.get_global_option("BASIS")
    e_hf, w_hf = psi4.energy('hf/%s' % basis, molecule=mol, return_wfn=True)
    w_hf.set_basisset("BASIS_DF_INT", basis_gdf_int)

    # are the custom defaults given?
    if exp_lower_bound is not None: DFBasis.exp_lower_bound = exp_lower_bound
    if exp_upper_bound is not None: DFBasis.exp_upper_bound = exp_upper_bound
    if ctr_lower_bound is not None: DFBasis.ctr_lower_bound = ctr_lower_bound
    if ctr_upper_bound is not None: DFBasis.ctr_upper_bound = ctr_upper_bound
 
    # fit the auxiliary basis
    if old_interface:
       dfbasis = DFBasis(w_hf.molecule(), templ_file=templ_file, param_file=param_file,                         
                                          bounds_file=bounds_file, constraints=constraints)
                                                                                                                
       OEP.read_vints = True
       oep     = OEP.create(oep_type, w_hf, dfbasis)
       #oep.compute()
       oep.compute_and_save_V()
       opt     = DFBasisOptimizer(oep)
                                                                                                                
       success = opt.fit(maxiter, tolerance, method, opt_global, temperature, stepsize, take_step, accept_test)
       
       # fetch the optimal basis
       basis_gdf_aux = opt.oep.dfbasis.basis
       
       # compute error for intermediate basis
       Z_aux = opt.compute_error(basis_gdf_aux, rms=True)
       Z_xpl = opt.compute_error(basis_gdf_xpl, rms=True)
       Z_int = opt.compute_error(basis_gdf_int, rms=True)
                                                                                                                
       # print the test results
       print(" Error for auxiliary    basis = %14.6f  Size: %6d" % (Z_aux, basis_gdf_aux.nbf()))
       print(" Error for example      basis = %14.6f  Size: %6d" % (Z_xpl, basis_gdf_xpl.nbf()))
       print(" Error for intermediate basis = %14.6f  Size: %6d" % (Z_int, basis_gdf_int.nbf()))
       print(opt.oep.dfbasis.print())
                                                                                                                
       # save
       if out is not None: opt.oep.dfbasis.save(out)

       dfbasis_opt = opt.oep.dfbasis

    else:

       # Updated procedure

       target = "OCC" if oep_type.lower() == "efp2-rep" else "VIR"
       standardized_input = None
       if use_standardized_input: 
          standardized_input = StandardizedInput(mol, oep_type, use_standardized_input)

       dfbasis_opt, err, GB = oep_ao_basis_set_optimizer(\
           w_hf, interm=basis_gdf_int, test=basis_gdf_int, exemplary=basis_gdf_xpl, target=target, cpp=True, more_info=True,
              templ_file=templ_file, param_file=param_file, bound_file=bounds_file, constraints=constraints, outname=out,
              opt_global=opt_global, standardized_input=standardized_input) 

    return dfbasis_opt

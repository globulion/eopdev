#!-*- coding: utf-8 *-*
#
# @BEGIN LICENSE
#
# oepdev by Bartosz Błasiak (email: blasiak.bartosz@gmail.com).
# oepdev copy, use and redistribution not allowed without written 
# consent of the Project Administrator.
#
# oepdev: A plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util

def run_oepdev(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    oepdev can be called via :py:func:`~driver.energy`. 

    >>> energy('oepdev')

    """
    # setup
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None: 
       print " [1] Computing HF/SCF for the dimer."
       ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # grab molecule aggregate and all fragments (now only for dimer)
    print " [2] Preparing basis sets for the monomers."
    molecule      = ref_wfn.molecule()

    # case when OEP build is requested
    if psi4.core.get_global_option("OEPDEV_TARGET") == "OEP_BUILD":
           
       basis_df_oep = psi4.core.BasisSet.build(molecule  , "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"),
                                               puream=ref_wfn.basisset().has_puream())
       ref_wfn.set_basisset("BASIS_DF_OEP", basis_df_oep)

    # case when solver tast on wavefunction union is requested
    else:
      
       molecule_A     = molecule.extract_subsets(1)                                                               
       molecule_B     = molecule.extract_subsets(2)
                                                                                                                 
       basis_A        = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("BASIS"),
                                                 puream=ref_wfn.basisset().has_puream())
       basis_B        = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("BASIS"),
                                                 puream=ref_wfn.basisset().has_puream())
       basis_df_oep_A = psi4.core.BasisSet.build(molecule_A, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"),
                                                 puream=ref_wfn.basisset().has_puream())
       basis_df_oep_B = psi4.core.BasisSet.build(molecule_B, "BASIS", psi4.core.get_global_option("DF_BASIS_OEP"),
                                                 puream=ref_wfn.basisset().has_puream())
                                                                                                                 
       ref_wfn.set_basisset("BASIS_1"       , basis_A       )
       ref_wfn.set_basisset("BASIS_2"       , basis_B       )
       ref_wfn.set_basisset("BASIS_DF_OEP_1", basis_df_oep_A)
       ref_wfn.set_basisset("BASIS_DF_OEP_2", basis_df_oep_B)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    print " [3] Running OEPDev plugin."
    oepdev_wfn = psi4.core.plugin('oepdev.so', ref_wfn)

    return oepdev_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['oepdev'] = run_oepdev

# Auxiliary routines
def test(dummy=True):
    r"An empty test."
    if dummy: psi4.core.print_out("\n\n  ==> Running Empty Test <==\n\n")
    pass

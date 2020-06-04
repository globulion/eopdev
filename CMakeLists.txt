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
cmake_policy(SET CMP0048 NEW)

# The Project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(oepdev 
	       VERSION 1.0.2 
	       LANGUAGES CXX)
set (OEPDEV_AUTHORS         "Bartosz Błasiak")
set (OEPDEV_DESCRIPTION     "Quantum Chemistry for Extended Molecular Aggregates")
set (OEPDEV_URL             "")
set (OEPDEV_LICENSE         "Private")
find_package(psi4 1.2 REQUIRED)

# External Numerical Library
set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
if( NOT EIGEN3_INCLUDE_DIR ) 
 message( FATAL_ERROR "
 OEPDev uses Eigen3 library for some numerical calculations that are not supported by Psi4 internal libraries.
 Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of Eigen3 library on your computer.
 Eigen3 can be downloaded for free from http://eigen.tuxfamily.org/ and does not require installation.
 Directory pointed by EIGEN3_INCLUDE_DIR should contain the directory called 'Eigen' (which should not be in the path).
 After settling with Eigen3, type:
 export EIGEN3_INCLUDE_DIR=<your path to Eigen3> 
 or add the above line to your .bashrc file and type 'source'.
")
endif()
include_directories ("${EIGEN3_INCLUDE_DIR}")

# Testing
include (CTest)

# Configure tests
include ("cmake/CTestTestfile.cmake.in")

include(psi4OptionsTools)

# Compile
add_psi4_plugin(oepdev main.cc 
	               # libutil
                       oepdev/libutil/util.cc
                       oepdev/libutil/diis.cc
                       oepdev/libutil/cphf.cc
                       oepdev/libutil/gram_schmidt.cc
                       oepdev/libutil/davidson_liu.cc
                       oepdev/libutil/wavefunction_union.cc
                       oepdev/libutil/integrals_iter.cc
                       oepdev/libutil/unitary_optimizer.cc
                       oepdev/libutil/scf_perturb.cc
                       oepdev/libutil/cis_base.cc
                       oepdev/libutil/cis_rhf.cc
                       oepdev/libutil/cis_rhf_explicit.cc
                       oepdev/libutil/cis_rhf_dl.cc
                       oepdev/libutil/cis_uhf.cc
                       oepdev/libutil/cis_uhf_explicit.cc
                       oepdev/libutil/cis_uhf_dl.cc
		       oepdev/libutil/cis_rhf_direct.cc
		       # lib3d
                       oepdev/lib3d/space3d.cc
                       oepdev/lib3d/esp.cc
                       oepdev/lib3d/dmtp_base.cc
                       oepdev/lib3d/dmtp_camm.cc
		       # libpsi
                       oepdev/libpsi/potential.cc
                       oepdev/libpsi/integral.cc
		       # libints
                       oepdev/libints/recurr.cc
                       oepdev/libints/eri.cc
		       # liboep
                       oepdev/liboep/oep_base.cc   
                       oepdev/liboep/oep_energy_coul.cc   
                       oepdev/liboep/oep_energy_pauli.cc   
                       oepdev/liboep/oep_energy_ct.cc   
                       oepdev/liboep/oep_coupling_eet.cc   
                       oepdev/liboep/oep_gdf.cc
		       # libgefp  
                       oepdev/libgefp/gefp.cc
                       oepdev/libgefp/gefp_efp2.cc
                       oepdev/libgefp/gefp_oep_efp2.cc
                       oepdev/libgefp/gefp_polar_base.cc
                       oepdev/libgefp/gefp_polar_uniform_base.cc
                       oepdev/libgefp/gefp_polar_uniform_field_1.cc
                       oepdev/libgefp/gefp_polar_uniform_field_2.cc
                       oepdev/libgefp/gefp_polar_nonuniform_base.cc
                       oepdev/libgefp/gefp_polar_nonuniform_field_1.cc
                       oepdev/libgefp/gefp_polar_nonuniform_field_2.cc
                       oepdev/libgefp/gefp_polar_nonuniform_field_1_grad_1.cc
                       oepdev/libgefp/gefp_polar_nonuniform_field_2_grad_1.cc
                       oepdev/libgefp/gefp_polar_abinitio.cc
                       oepdev/libgefp/gefp_polar_ffabinitio.cc
		       # libsolver
                       oepdev/libsolver/solver_base.cc
                       oepdev/libsolver/solver_energy_coul.cc
                       oepdev/libsolver/solver_energy_pauli.cc
                       oepdev/libsolver/solver_energy_ct.cc
                       oepdev/libsolver/solver_coupling_eet.cc
                       oepdev/libsolver/ti_data.cc
		       # libtest
                       oepdev/libtest/test.cc
		       oepdev/libtest/basic.cc
                       oepdev/libtest/cis_rhf_explicit.cc
                       oepdev/libtest/cis_rhf_dl.cc
                       oepdev/libtest/cis_uhf_explicit.cc
                       oepdev/libtest/cis_uhf_dl.cc
		       oepdev/libtest/cphf.cc
		       oepdev/libtest/camm.cc
		       oepdev/libtest/dmtp_energy.cc
		       oepdev/libtest/efp2_energy.cc
		       oepdev/libtest/eri_2_2.cc
		       oepdev/libtest/eri_3_1.cc
		       oepdev/libtest/eri_1_1.cc
		       oepdev/libtest/scf_perturb.cc
		       oepdev/libtest/dmatpol.cc
		       oepdev/libtest/dmatpolX.cc
		       oepdev/libtest/unitary_optimizer.cc
		       oepdev/libtest/unitary_optimizer_4_2.cc
                       oepdev/libtest/esp_solver.cc
                       oepdev/libtest/test_custom.cc
		       oepdev/libtest/points_collection3d.cc
		       oepdev/libtest/ct_energy_benchmark_ol.cc
		       oepdev/libtest/rep_energy_benchmark_hs.cc
		       oepdev/libtest/rep_energy_benchmark_dds.cc
		       oepdev/libtest/rep_energy_benchmark_efp2.cc
		       oepdev/libtest/rep_energy_benchmark_ol.cc
		       oepdev/libtest/rep_energy_benchmark_murrell_etal.cc
		       oepdev/libtest/rep_energy_oep_based_murrell_etal.cc
		       # exports
		       oepdev/export/export_dmtp.cc
		       oepdev/export/export_cphf.cc
		       oepdev/export/export_solver.cc
		       oepdev/export/export_util.cc
                       oepdev/export/export_oep.cc
                       oepdev/export/export_gefp.cc
               )

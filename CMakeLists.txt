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

# The Project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(oepdev CXX)
set (OEPDEV_VERSION_MAJOR   "1")
set (OEPDEV_VERSION_MINOR   "0")
set (OEPDEV_VERSION_PATCH   "0")
set (OEPDEV_VERSION "${OEPDEV_VERSION_MAJOR}.${BQC_VERSION_MINOR}.${BQC_VERSION_PATCH}")
find_package(psi4 1.0 REQUIRED)

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

# Compile
add_psi4_plugin(oepdev main.cc 
                       oepdev/libutil/util.cc
                       oepdev/libutil/diis.cc
                       oepdev/libutil/cphf.cc
                       oepdev/libutil/wavefunction_union.cc
                       oepdev/libutil/integrals_iter.cc
                       oepdev/libutil/space3d.cc
                       oepdev/libutil/esp.cc
                       oepdev/libutil/solver.cc
                       oepdev/libutil/unitary_optimizer.cc
                       oepdev/libpsi/potential.cc
                       oepdev/libpsi/integral.cc
                       oepdev/liboep/oep.cc   
                       oepdev/libints/recurr.cc
                       oepdev/libints/eri.cc
                       oepdev/libtest/test.cc
                       oepdev/libgefp/gefp.cc
               )

# Configure tests
add_test (NAME BASIC   COMMAND "psi4" "tests/oepdev/basic.inp"   "tests/oepdev/basic.out"  )
add_test (NAME CPHF    COMMAND "psi4" "tests/oepdev/cphf.inp"    "tests/oepdev/cphf.out"   )
add_test (NAME DMATPOL COMMAND "psi4" "tests/oepdev/dmatpol.inp" "tests/oepdev/dmatpol.out")
add_test (NAME ERI_1_1 COMMAND "psi4" "tests/oepdev/eri_1_1.inp" "tests/oepdev/eri_1_1.out")
add_test (NAME ERI_2_2 COMMAND "psi4" "tests/oepdev/eri_2_2.inp" "tests/oepdev/eri_2_2.out")
add_test (NAME ERI_3_1 COMMAND "psi4" "tests/oepdev/eri_3_1.inp" "tests/oepdev/eri_3_1.out")
add_test (NAME UNITARY_OPTIMIZER COMMAND "psi4" "tests/oepdev/unitaryOptimizer.inp" "tests/oepdev/unitaryOptimizer.out")
add_test (NAME UNITARY_OPTIMIZER_4_2 COMMAND "psi4" "tests/oepdev/unitaryOptimizer_4_2.inp" "tests/oepdev/unitaryOptimizer_4_2.out")

set_tests_properties (BASIC   PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.0000"    )
set_tests_properties (CPHF    PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.00"      )
set_tests_properties (DMATPOL PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.00"      )
set_tests_properties (ERI_1_1 PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.00000000")
set_tests_properties (ERI_2_2 PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.00000000")
set_tests_properties (ERI_3_1 PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.00000000")
set_tests_properties (UNITARY_OPTIMIZER PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.000")
set_tests_properties (UNITARY_OPTIMIZER_4_2 PROPERTIES PASS_REGULAR_EXPRESSION " Test result= 0.000")

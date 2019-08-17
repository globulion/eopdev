#!/usr/bin/python3
from gefp.core.utilities import wavefunction_union_from_dimer
from gefp.core.utilities import wavefunction_union_from_dfi_solver
from gefp.density.dfi import DFI_JK as DFI

import oepdev
import numpy
import psi4
import sys

# tentative test molecule resembling etylene dimer from Fujimoto's work R=3.0 Angstrom
inp = ("""
0 1
C     -0.031871   -0.000007    0.001312
H     -0.053949   -0.102551    1.069139
C      1.188324    0.000003   -0.695389
H     -0.962724   -0.100044   -0.523396
H      1.209504   -0.099215   -1.763690
H      2.119240   -0.102685   -0.171502
--
0 1
C     -0.031871    3.000007    0.001312
H     -0.053949    3.102551    1.069139
C      1.188324    2.999997   -0.695389
H     -0.962724    3.100044   -0.523396
H      1.209504    3.099215   -1.763690
H      2.119240    3.102685   -0.171502

units angstrom
symmetry c1
noreorient
nocom
""")

# read molecule
dimer = psi4.geometry(inp)

# run dfi?
run_dfi = True

# set Psi4 options
psi4.set_options({"scf_type"       : "df"    ,
                  "basis"          : "6-31G*",
                  "e_convergence"  : 1e-9    ,
                  "puream"         : False   ,
                  "print"          : 0       ,
                  "excited_state_a":-1       ,
                  "excited_state_b":-1       ,
                  "trcamm_symmetrize": False})

# set Psi4 output
psi4.core.set_output_file(sys.argv[0].replace('.py','.log'), True)

# construct wavefunction union 
if run_dfi:
   dfi = DFI(dimer)
   dfi.run()
   un = wavefunction_union_from_dfi_solver(dfi)
else:
   un = wavefunction_union_from_dimer(dimer)


# perform 4-index transformations in dimer spaces
un.transform_integrals()

# run TrCAMM and TDFI-TI
solver = oepdev.OEPDevSolver.build("EET COUPLING CONSTANT", un)
solver.compute_benchmark("FUJIMOTO_TI_CIS")

#!/usr/bin/python3
"""
 Analyse ethylene dimer from Fujimoto JCP 2012 (TDFI-TI work)
 Equilibrium geometry: MP2/6-31G* (g09)
 Source: data@pauli.ch.pwr.wroc.pl:/home/data/2015.EET/logs/C2H4/opt/c2h4_2@mp2_6-31d_d2h.log

 Tasks:
  - compute TDFI-TI(CIS) and TrCAMM EET coupling between first bright states of ethylene in a symmetric dimer

 BB, Gundelfingen 22 Aug 2019
"""
from gefp.core.utilities import wavefunction_union_from_dimer
from gefp.core.utilities import wavefunction_union_from_dfi_solver
from gefp.density.dfi import DFI_J as DFI

import oepdev
import numpy
import psi4
import sys

# ethylene dimer from Fujimoto's JCP 2012 paper (TDFI-TI). 
inp = ("""
0 1
C
C  1  rcc
H  1  rch  2  ahcc
H  2  rch  1  ahcc  3  d1
H  2  rch  1  ahcc  3  d2
H  1  rch  2  ahcc  5  d1

units angstrom
symmetry c1
no_reorient
no_com
--
0 1
C  1  ro   2  a1    4  d3
C  7  rcc  1  a1    2  d1
H  7  rch  8  ahcc  2  d3
H  8  rch  7  ahcc  9  d1
H  8  rch  7  ahcc  9  d2
H  7  rch  8  ahcc  11 d1

units angstrom
symmetry c1
no_reorient
no_com

rcc = 1.33627
rch = 1.08510
ro  = 4.16896
ahcc= 121.693
a1  = 90.0000
d1  = 0.00000
d2  = 180.000
d3  = 89.906
""")

def good_dimer(xyz,n):
    'Format dimer and return psi4.core.Molecule from xyz string'
    xyz = xyz.split('\n')
    xyz.pop()
    xyz.insert(n+1,'--')
    xyz.insert(n+2,xyz[0])
    xyz.insert(0,'\n')
    xyz.append('units angstrom')
    xyz.append('symmetry c1')
    xyz.append('no_reorient')
    xyz.append('nocom')
    xyz = '\n'.join(xyz)
    xyz = psi4.geometry(xyz)
    xyz.update_geometry()
    return xyz


# read molecule
dimer = psi4.geometry(inp)
dimer.ro = 3.0 # Angstroms
dimer.update_geometry()
dimer = good_dimer(dimer.save_string_xyz(), 6)
dimer.print_cluster()

# run dfi?
run_dfi = 0

# set Psi4 options
psi4.set_options({"scf_type"       : "df"    ,
                  "basis"          : "6-31G*",
                  "e_convergence"  : 1e-9    ,
                  "puream"         : False   ,
                  "print"          : 0       ,
                  "excited_state_a":-1       ,
                  "excited_state_b":-1       ,
                  "trcamm_symmetrize": True , 
                  "ti_cis_scf_fock_matrix": False, 
                  "ti_cis_print_fock_matrix": False,})

# set Psi4 output
psi4.core.set_output_file(sys.argv[0].replace('.py','.log'), True)

# construct wavefunction union 
if run_dfi:
   dfi = DFI(dimer)
   dfi.run()
   un = wavefunction_union_from_dfi_solver(dfi)
else:
   un = wavefunction_union_from_dimer(dimer)

# run TrCAMM and TDFI-TI
solver = oepdev.OEPDevSolver.build("EET COUPLING CONSTANT", un)
solver.compute_benchmark("FUJIMOTO_TI_CIS")

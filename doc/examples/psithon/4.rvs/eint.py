#!/usr/bin/python3
"""
 Counterpoise corrected interaction energies. Benchmark.

 Usage: [basis] [xyz.psi]
"""

import psi4, gefp, sys, os

psi4.set_options({"scf_type"        : "direct"           ,
                  "guess"           : "auto"             ,
                  "df_scf_guess"    : True               ,
                  "e_convergence"   : 1e-10              ,
                  "d_convergence"   : 1e-10              ,
                  "basis"           : sys.argv[1]        ,
                  "puream"          : False              ,
                  "print"           : 1                  ,})
#psi4.set_output_file("eint.log", False)

mol = gefp.core.utilities.psi_molecule_from_file(sys.argv[2])
gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])

e = psi4.energy('scf', bsse_type=['nocp', 'cp', 'vmfc'], molecule=mol)

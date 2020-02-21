#!/usr/bin/python3
"""
 Test for reproducing results from Steven & Fink, CPL 139, pp.: 15-22 (1987)
 
 Usage: [basis] [xyz.psi]
"""

import psi4, gefp, sys

psi4.set_options({"scf_type"        : "direct"   ,
                  "guess"           : "auto"     ,
                  "df_scf_guess"    : True       ,
                  "e_convergence"   : 1e-10     ,
                  "d_convergence"   : 1e-10     ,
                   "basis"          : sys.argv[1],
                  "puream"          : False      ,
                  "print"           : 1          ,})
psi4.set_output_file("test.log", False)

mol = gefp.core.utilities.psi_molecule_from_file(sys.argv[2])
gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])

e_hf, wfn = psi4.energy('scf', molecule=mol, return_wfn=True)


sol = gefp.density.rvs.RVS(wfn) 
sol.run_dimer(conver=1e-8)
psi4.core.print_out(str(sol))

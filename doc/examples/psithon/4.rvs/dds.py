#!/usr/bin/python3
import psi4, gefp, sys
"""
 DDS interaction energies. Benchmark for CEX - run in MCBS mode.

 Usage: [basis] [xyz.psi]
"""

psi4.set_options({"scf_type"        : "direct"   ,
                  "guess"           : "auto"     ,
                  "df_scf_guess"    : True       , 
                  "e_convergence"   : 1e-10      ,
                  "d_convergence"   : 1e-10      ,
                  "basis"           : sys.argv[1], 
                  "puream"          : False      ,
                  "print"           : 1          ,})
psi4.set_output_file("dds.log", False)

mol = gefp.core.utilities.psi_molecule_from_file(sys.argv[2])
gefp.math.matrix.move_atom_rotate_molecule(mol, [-30., -83., 43.1])

dds = gefp.density.partitioning.DensityDecomposition(mol, method='hf', 
     acbs=False, jk_type='direct', no_cutoff=0.000, xc_scale=1.0, l_dds=False,
                       cc_relax=False, verbose=False, n_eps=5.0E-5)
dds.compute(polar_approx=False)
v = dds.vars
print(dds)

conv = psi4.constants.hartree2kcalmol

CEX = v["e_cou_t"] + v["e_exr_t"]
print(" CEX = %16.6f [kcal/mol]" % (CEX*conv))

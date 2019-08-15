#!/usr/bin/python3
"""
 Attempts to recover plot in Fig.4 from Fujimoto and Yang JCP 129,
 054102 (2008). In fact, Fujimoto used B3LYP/6-31G(d) with Coulomb
 potential.

 Plot dfi_jk.pdf shows the artifacts at short distances if exchange
 integrals are introduced in DFI-HF/6-31(d). If it does not work at
 HF level it probably won't work at DFT level.
"""

import psi4
import sys
import numpy

sys.path.append('..')

from gefp.density.dfi import SCF, DFI_JK as DFI

psi4.set_memory('2 GB')

h4o2 = ("""
0 1
 h
 o    1 oh2     
 h    2 ho3         1 hoh3      

oh2   =     0.967887
ho3   =     0.976947
hoh3  =     104.053
units angstrom
symmetry c1
no_reorient
nocom
--
0 1
 o    3 {0}         2 oho4          1 dih4   
 h    4 ho5         3 hoh5          2 dih5   
 h    4 ho6         3 hoh6          2 dih6   
 
oh4   =     1.919163
oho4  =     161.059
dih4  =    -179.048
ho5   =     0.970444
hoh5  =      95.003
dih5  =      51.351
ho6   =     0.970429
hoh6  =      95.409
dih6  =     -53.167
units angstrom
symmetry c1
no_reorient
nocom
""")

Rs=[1.919163]
Rs.extend(numpy.arange(1.0,2.5,0.1))
Rs.extend(numpy.arange(3.,10.,1.))
Rs.sort()

def header(s):
    header = str('\n'+len(s)*'-'+'\n'+s+'\n'+len(s)*'-'+'\n')
    psi4.core.print_out(header)

header("PES scan for the following values of R") 
psi4.core.print_out(", ".join(['{:5.3f}'.format(v) for v in Rs])+'\n\n')

# set the global options
psi4.set_options({'basis'    : '6-31G*',
                  'scf_type' : 'direct',
                  'guess'    : 'core'})

# run the program quietly
psi4.core.set_output_file(sys.argv[0].replace('.py','.log'), True)

def dimer(xyz,n):
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
    return xyz

E_DFI = {}
E_SCF = {}

# perform scan 
for R in Rs:
    geom = psi4.geometry(h4o2.format(R))
    geom.update_geometry()
    geom = geom.save_string_xyz()
    geom = dimer(geom,3)
    geom.update_geometry()
    geom.print_cluster()
    #geom.reset_point_group('c1')
    #geom.fix_orientation(True)
    #geom.print_cluster()
    dfi = DFI(geom)
    E_DFI[R] = dfi.run()
    psi4.core.clean()
    E_SCF[R], wfn = psi4.energy('scf', bsse_type='cp', return_wfn=True, molecule=geom)
    psi4.core.clean()

# generate plot
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline

Table_DFI = numpy.array(sorted(E_DFI.items()))
print(Table_DFI)
R_dfi = Table_DFI[:,0]
E_dfi = Table_DFI[:,1]
R_dfi_i = numpy.linspace(R_dfi.min(),R_dfi.max(),500)
E_dfi_i = make_interp_spline(R_dfi, E_dfi, k=3)

Table_SCF = numpy.array(sorted(E_SCF.items()))
print(Table_SCF)
R_scf = Table_SCF[:,0]
E_scf = Table_SCF[:,1]
R_scf_i = numpy.linspace(R_scf.min(),R_scf.max(),500)
E_scf_i = make_interp_spline(R_scf, E_scf, k=3)

au2kcal = psi4.constants.hartree2kcalmol
header("{:>14s} {:>14s} {:>14s} [A \ kcal/mol]\n".format('R', 'E(RHF)', 'E(DFI)'))
for i in range(len(Table_SCF)):
    psi4.core.print_out("{:>14.3f} {:>14.3f} {:>14.3f}\n".format(Table_SCF[i][0], Table_SCF[i][1]*au2kcal, Table_DFI[i][1]*au2kcal))
header("{:>14s} {:>14s} {:>14s} [A \ au]\n".format('R', 'E(RHF)', 'E(DFI)'))
for i in range(len(Table_SCF)):
    psi4.core.print_out("{:>14.3f} {:>14.3f} {:>14.3f}\n".format(Table_SCF[i][0], Table_SCF[i][1], Table_DFI[i][1]))

fig = plt.figure()
plt.xlabel('R$_{H-O}$ [A]')
plt.ylabel('E [a.u.]')
plt.plot(R_dfi, E_dfi, 'bo', markersize=12)
plt.plot(R_scf, E_scf, 'ro', markersize=12)
plt.plot(R_dfi_i, E_dfi_i(R_dfi_i), 'b-', linewidth=2, label="DFI")
plt.plot(R_scf_i, E_scf_i(R_scf_i), 'r-', linewidth=2, label="SCF")
plt.legend(loc=1)
plt.show()
#fig.savefig( 'sys.argv[0].replace('.py','.pdf') )
fig.savefig('dfi_jk.pdf')


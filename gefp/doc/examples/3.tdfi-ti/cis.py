#!/usr/bin/python3
"""
 Analyse ethylene dimer from Fujimoto JCP 2012 (TDFI-TI work)
 Equilibrium geometry: MP2/6-31G* (g09)
 Source: data@pauli.ch.pwr.wroc.pl:/home/data/2015.EET/logs/C2H4/opt/c2h4_2@mp2_6-31d_d2h.log

 Tasks:
  - compute CIS energies and amplitudes for isolated monomer taken from equilibrium geometry
  - compute CIS energies for dimer in equilibrium geometry
  - compute CIS energies and amplitudes for dimer in infinite separation of monomers
  - perform t1 scan (R_min = 3.0, R_max = 9.0 Angs) by translating dimer maintaining D2H symmetry. Compute CIS energies
    Writes "energies.dat" file with excitation energies (eV) in each scan for first 10 states 

 By comparing the CIS amplitudes of the monomer and in infinitely stretched dimer, the scaling factor is
   CIS(dimer) = CIS(monomer) / sqrt(2)

 From scan we can determine dimer excited states for computation of the exact couplings as 1/2 * (E - E')

 BB, Gundelfingen 22 Aug 2019
"""
import oepdev
import numpy
import psi4
import gefp
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
--
0 1
C  1  ro   2  a1    4  d3
C  7  rcc  1  a1    2  d1
H  7  rch  8  ahcc  2  d3
H  8  rch  7  ahcc  9  d1
H  8  rch  7  ahcc  9  d2
H  7  rch  8  ahcc  11 d1

rcc = 1.33627
rch = 1.08510
ro  = 4.16896
ahcc= 121.693
a1  = 90.0000
d1  = 0.00000
d2  = 180.000
d3  = 89.906

units angstrom
symmetry c1
noreorient
nocom
""")

class CIS_WFN:
  def __init__(self, mol, idx=None):
      e, w = psi4.energy('scf', return_wfn=True, molecule=mol)
      cis  = oepdev.CISComputer.build("RESTRICTED", w, psi4.core.get_options(), "RHF")
      cis.compute()
      cis.clear_dpd()
      psi4.core.clean()
      if idx is None: idx = numpy.arange(cis.nstates())
      wfn  = gefp.density.ci.CIS_CIWavefunction(w, cis.E().to_array(dense=True)[idx], cis.U().to_array(dense=True)[:,idx])
      self.cis = cis
      self.wfn = wfn

  def print_amplitudes(self, i, thr=0.1):
      na = self.wfn.ref_wfn.nalpha()
      nb = self.wfn.ref_wfn.nbeta()
      nmo= self.wfn.ref_wfn.nmo()
      na_occ = na
      nb_occ = nb
      na_vir = nmo - na
      nb_vir = nmo - nb
      C = self.wfn.ci_c[:,i]
      print(" Alpha Amplitudes:")
      for i in range(na_occ):
          for a in range(na_vir):
              ia = na_vir * i + a
              c = C[ia]
              if abs(c)>thr:
                 print(" Ta %2d --> %2d  = %14.6f" % (i+1, a+1, c))
      print(" Beta Amplitudes:")
      off = na_occ*na_vir
      for i in range(nb_occ):
          for a in range(nb_vir):
              ia = off + na_vir * i + a
              c = C[ia]
              if abs(c)>thr:
                 print(" Tb %2d --> %2d  = %14.6f" % (i+1, a+1, c))



# set Psi4 options
psi4.set_options({"scf_type"       : "df"    ,
                  "basis"          : "6-31G*",
                  "e_convergence"  : 1e-9    ,
                  "puream"         : False   ,
                  "print"          : 1       })

# set Psi4 output
psi4.core.set_output_file(sys.argv[0].replace('.py','.log'), True)

# print only first 10 excited states
N = 10

# CIS for monomer
mol = psi4.geometry(inp).extract_subsets(1); mol.update_geometry()
m = CIS_WFN(mol, idx=numpy.arange(N))
print(" Monomer")
for i in range(N):
    print("E  = %14.2f f = %14.3f c = %14.5f" % (m.wfn.ci_e[i]*psi4.constants.hartree2ev, 
                                                 m.cis.oscillator_strength(i), m.cis.U_homo_lumo(i, 0, 0)[0]))
print()
print(" CIS amplitudes for state 2")
m.print_amplitudes(1)

# CIS dimer: equilibrium geometry
mol = psi4.geometry(inp)
mol.update_geometry()
d = CIS_WFN(mol, idx=numpy.arange(N))
print(" Dimer (equilibrium)")
for i in range(N):
    print("E  = %14.4f f = %14.3f c = %14.5f" % (d.wfn.ci_e[i]*psi4.constants.hartree2ev, 
                                                 d.cis.oscillator_strength(i),d.cis.U_homo_lumo(i, 0, 0)[0]))
print()


# CIS dimer: large separation (monomers not interacting)
mol = psi4.geometry(inp)
mol.ro = 13.0
mol.update_geometry()
d = CIS_WFN(mol, idx=numpy.arange(N))
print(" Dimer (equilibrium)")
for i in range(N):
    print("E  = %14.4f f = %14.3f c = %14.5f" % (d.wfn.ci_e[i]*psi4.constants.hartree2ev, 
                                                 d.cis.oscillator_strength(i),d.cis.U_homo_lumo(i, 0, 0)[0]))
print()
print(" CIS amplitudes for state 3")
d.print_amplitudes(2)
print(" CIS amplitudes for state 4")
d.print_amplitudes(3)


# CIS dimer: scan
out = open('energies.dat','w')
for ro in numpy.linspace(3.0, 9.0, 100):
    dimer = psi4.geometry(inp); dimer.ro = ro; dimer.update_geometry()
    s = CIS_WFN(dimer, idx=numpy.arange(N))
    f = numpy.zeros(N)
    for i in range(N):
        f[i] = s.cis.oscillator_strength(i)
    print("r = %14.2f" % ro + 10*"%6.3f" % tuple(f))
    e = s.wfn.ci_e * psi4.constants.hartree2ev; del s
    log = "%14.5E" % ro
    log+= len(e)*"%14.6E" % tuple(e)
    out.write(log + '\n')
out.close()

#e, c, s = gefp.math.matrix.rearrange_eigenpairs(u=c_R.wfn.ci_c, u_ref=c_0.wfn.ci_c, n=c_R.wfn.ci_e, return_sim=True)
#print(s)

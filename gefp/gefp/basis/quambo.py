#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 QUAMBO module.
 Bartosz Błasiak, Gundelfingen, 10 July 2020
"""
import psi4
import numpy
from ..math.matrix import matrix_power

__all__ = ["QUAMBO"]

class QUAMBO:
  """
 Quasiatomic Minimal Basis Set Molecular Orbitals

 Usage:

 solve = QUAMBO(molecule, psi4.core.get_options(), acbs=False)
 solve.compute()
 c_vir = solve.mo(space='vir')
 e_all = solve.eps(space='all')
 quambo_b_nonorthogonal = solve.quambo(type='nonorthogonal', spin='beta')

 Notes:
  o first argument to QUAMBO constructor can also be Wavefunction object
    with the SCF solution of a molecule. In this case, SCF is not run 
    for the molecule.

 References: 

  [1] W. C. Lu, C. Z. Wang, M.W. Schmidt, L. Bytautas, K. M. Ho, K. Reudenberg, 
      J. Chem. Phys. 120, 2629 (2004) [original QUAMBO paper]
  [2] P. Xu, M. S. Gordon, J. Chem. Phys. 139, 194104 (2013) [application 
      of QUAMBO in EFP2 CT term]
"""
  def __init__(self, mol, options, acbs=False):
      self._options = options
      if isinstance(options, dict): 
         psi4.core.set_options(options)
         self._options = psi4.core.get_options()
      self.acbs = acbs

      # ---> Namespace of protected constant attributes <--- #

      # numbers of minimal basis functions of free atoms
      self._nbas_atom_mini= { "H": 1, "He": 1, 
                             "Li": 5, "Be": 5, "B": 5, "C": 5, "N": 5, "O": 5, "F": 5, "Ne": 5,
                             "Na": 9, "Mg": 9,"Al": 9,"Si": 9, "P": 9, "S": 9,"Cl": 9, "Ar": 9, 
                             "K" :13, "Ca":13, }

      # numbers of unpaired electrons in free atoms
      self._unpe_atom     = { "H": 1, "He": 0, 
                             "Li": 1, "Be": 0, "B": 1, "C": 4, "N": 3, "O": 2, "F": 1, "Ne": 0, # assume C* state
                             "Na": 1, "Mg": 0,"Al": 1,"Si": 4, "P": 3, "S": 2,"Cl": 1, "Ar": 0,  
                             "K" : 1, "Ca": 0, }


      # ---> Namespace of protected calculables <--- #

      # AO Overlap Matrix
      self._Sao = None
      # QUAMBO (Alpha, non-orthogonal) 
      self._quambo_a_nonorthogonal = None
      # QUAMBO (Alpha, orthogonal) 
      self._quambo_a_orthogonal    = None
      # QUAMBO (Beta, non-orthogonal) 
      self._quambo_b_nonorthogonal = None
      # QUAMBO (Beta, orthogonal) 
      self._quambo_b_orthogonal    = None
      # Virtual Valence Molecular Orbitals (Alpha, VVO)
      self._c_a_mini_vir = None
      # Virtual Valence Molecular Orbitals (Beta, VVO)
      self._c_b_mini_vir = None
      # VVO Energies (Alpha)
      self._e_a_mini_vir = None
      # VVO Energies (Beta)
      self._e_b_mini_vir = None
      # All Molecular orbitals (Alpha, OCC + VVO)
      self._c_a_mini = None 
      # All Molecular orbitals (Beta, OCC + VVO)
      self._c_b_mini = None
      # Energies of All Molecular Orbitals (Alpha)
      self._e_a_mini = None
      # Energies of All Molecular Orbitals (Beta)
      self._e_b_mini = None

      # Determine the molecule
      if isinstance(mol, psi4.core.Molecule):
         self._mol = mol
         # Run SCF of a molecule
         en, self._wfn = psi4.energy('scf', molecule=self._mol, return_wfn = True); psi4.core.clean()

      elif isinstance(mol, psi4.core.Wavefunction):
         self._mol = mol.molecule()
         self._wfn = mol

      else: raise TypeError("Wrong data type provided to QUAMBO constructor first argument")


  def quambo(self, spin="alpha", type="orthogonal"):
      "QUAMBO minimal basis set"
      if type.lower().startswith("non"):
         if spin.lower().startswith("alpha"): return self._quambo_a_nonorthogonal
         else                               : return self._quambo_b_nonorthogonal
      elif type.lower().startswith("ort"):
         if spin.lower().startswith("alpha"): return self._quambo_a_orthogonal
         else                               : return self._quambo_b_orthogonal
      else: raise ValueError("Wrong options for QUAMBO type!")

  def mo(self, spin="alpha", space="all"):
      "Molecular orbitals in minimal QUAMBO basis"
      if   space.lower() == "all": 
        if spin.lower().startswith("a"): return self._c_a_mini
        else                           : return self._c_b_mini
      elif space.lower() == "occ": raise NotImplementedError
      elif space.lower() == "vir": 
        if spin.lower().startswith("a"): return self._c_a_mini_vir
        else                           : return self._c_b_mini_vir
      else:
        raise ValueError("Wrong orbital space chosen")

  def eps(self, spin="alpha", space="all"):
      "Energies of molecular orbitals in minimal QUAMBO basis"
      if   space.lower() == "all": 
        if spin.lower().startswith("a"): return self._e_a_mini
        else                           : return self._e_b_mini
      elif space.lower() == "occ": raise NotImplementedError
      elif space.lower() == "vir": 
        if spin.lower().startswith("a"): return self._e_a_mini_vir
        else                           : return self._e_b_mini_vir
      else:
        raise ValueError("Wrong orbital space chosen")

  def overlap(self, spin="alpha", type="nonorthogonal"):
      "Compute overlap matrix between QUAMBOs"
      if type.lower().startswith("non"):
         if spin.lower().startswith("alpha"): q = self._quambo_a_nonorthogonal
         else                               : q = self._quambo_b_nonorthogonal
      elif type.lower().startswith("ort"):
         if spin.lower().startswith("alpha"): q = self._quambo_a_orthogonal
         else                               : q = self._quambo_b_orthogonal
      else: raise ValueError("Wrong options for type!")
      S = q.T @ self._Sao @ q
      return S

  def compute(self):
      "Calculate QUAMBO expansion coefficients -a- and transformation matrix -T-"

      psi4.core.print_out("\n ===> Computing QUAMBOs <===\n\n")

      # [1] Read the orbitals and Fock matrix
      ca_occ = self._wfn.Ca_subset("AO","OCC").to_array(dense=True)
      ca_vir = self._wfn.Ca_subset("AO","VIR").to_array(dense=True)
      cb_occ = self._wfn.Cb_subset("AO","OCC").to_array(dense=True)
      cb_vir = self._wfn.Cb_subset("AO","VIR").to_array(dense=True)

      bfs = self._wfn.basisset()

      eps_a_occ = self._wfn.epsilon_a_subset("MO", "OCC").to_array(dense=True)
      eps_a_vir = self._wfn.epsilon_a_subset("MO", "VIR").to_array(dense=True)
      eps_a_all = self._wfn.epsilon_a_subset("MO", "ALL").to_array(dense=True)
      eps_b_occ = self._wfn.epsilon_b_subset("MO", "OCC").to_array(dense=True)
      eps_b_vir = self._wfn.epsilon_b_subset("MO", "VIR").to_array(dense=True)
      eps_b_all = self._wfn.epsilon_b_subset("MO", "ALL").to_array(dense=True)
     
      Fa = self._wfn.Fa().to_array(dense=True)
      Fb = self._wfn.Fb().to_array(dense=True)

      Sao= self._wfn.S ().to_array(dense=True); self._Sao = Sao
      #X  = matrix_power(Sao, -0.5); self._X = X
      #Y  = matrix_power(Sao,  0.5); self._Y = Y

      # sizing
      ndocc = self._wfn.doccpi()[0]
      nsocc = self._wfn.soccpi()[0]

      naocc = ca_occ.shape[1]
      navir = ca_vir.shape[1]
      nbocc = cb_occ.shape[1]
      nbvir = cb_vir.shape[1]

      naocc_mini= self._wfn.nalpha()          ; self._naocc_mini = naocc_mini
      nbocc_mini= self._wfn.nbeta()           ; self._nbocc_mini = nbocc_mini
      nbas_mini = self._calculate_nbas_mini() ; self. _nbas_mini =  nbas_mini
      navir_mini= nbas_mini - naocc_mini      ; self._navir_mini = navir_mini
      nbvir_mini= nbas_mini - nbocc_mini      ; self._nbvir_mini = nbvir_mini

      psi4.core.print_out("\n QUAMBO orbital analysis:\n")
      psi4.core.print_out  (" Na_occ = %5d    Na_occ_mini = %5d\n" % (naocc,naocc_mini))
      psi4.core.print_out  (" Nb_occ = %5d    Nb_occ_mini = %5d\n" % (nbocc,nbocc_mini))
      psi4.core.print_out  (" Na_vir = %5d    Na_vir_mini = %5d\n" % (navir,navir_mini))
      psi4.core.print_out  (" Nb_vir = %5d    Nb_vir_mini = %5d\n" % (nbvir,nbvir_mini))
      psi4.core.print_out  (" There are %5d minimal basis molecular orbitals (QUAMBOs) for alpha and beta spin\n" % nbas_mini)
      psi4.core.print_out  (" Looking for %5d virtual ALPHA valence orbitals (VVOs)\n"   % navir_mini)
      psi4.core.print_out  (" Looking for %5d virtual BETA  valence orbitals (VVOs)\n\n" % nbvir_mini)
      # 
      #Fa = X @ Fa @ X
      #Fb = X @ Fb @ X
      #ca_occ = Y @ ca_occ
      #cb_occ = Y @ cb_occ
      #ca_vir = Y @ ca_vir
      #cb_vir = Y @ cb_vir


      # [2] Run ROHF for all free atoms
      reference_stash = self._options.get_str("REFERENCE")
      mints = psi4.core.MintsHelper(bfs)
      psi4.set_options({"reference":"rohf"})

      a_occ_a_set = numpy.zeros((naocc, 0)) ; a_occ_b_set = numpy.zeros((nbocc, 0))
      a_vir_a_set = numpy.zeros((navir, 0)) ; a_vir_b_set = numpy.zeros((nbvir, 0))

      for atom in self._atomize():
          en_a, wfn_a = psi4.energy('hf', molecule=atom, return_wfn = True); psi4.core.clean()
         
          # A* orbitals (free-atom minimal basis occupied valence+core orbitals)
         #ca_a  = wfn_a.Ca_subset("AO","ALL").to_array(dense=True)[:,:self._nbas_atom_mini[atom.symbol(0)]]
          ca_a  = wfn_a.Ca_subset("AO","OCC").to_array(dense=True) 

          # a* coefficients (projections of A* onto molecule's occupied and virtual orbitals)
          if self.acbs:
              a_occ_a = ca_occ.T  @ Sao @ ca_a  ; a_occ_b = cb_occ.T  @ Sao @ ca_a 
              a_vir_a = ca_vir.T  @ Sao @ ca_a  ; a_vir_b = cb_vir.T  @ Sao @ ca_a
          else:
              # metric 
              bfs_a = wfn_a.basisset()
              S_a = mints.ao_overlap(bfs, bfs_a).to_array(dense=True)

              a_occ_a = ca_occ.T  @ S_a @ ca_a  ; a_occ_b = cb_occ.T  @ S_a @ ca_a 
              a_vir_a = ca_vir.T  @ S_a @ ca_a  ; a_vir_b = cb_vir.T  @ S_a @ ca_a

         #atom.print_out()

          a_occ_a_set = numpy.hstack((a_occ_a_set.copy(), a_occ_a.copy()))
          a_vir_a_set = numpy.hstack((a_vir_a_set.copy(), a_vir_a.copy()))
          a_occ_b_set = numpy.hstack((a_occ_b_set.copy(), a_occ_b.copy()))
          a_vir_b_set = numpy.hstack((a_vir_b_set.copy(), a_vir_b.copy()))


      # [3] Compute QUAMBOs 
      quambo_a_orthogonal, quambo_a_nonorthogonal, e_a_mini_vir, c_a_mini_vir, e_a_mini, c_a_mini = \
         self._compute_quambo(ca_occ, ca_vir, eps_a_occ, eps_a_vir, Fa, a_occ_a_set, a_vir_a_set, naocc_mini, navir_mini, "ALPHA")

      quambo_b_orthogonal, quambo_b_nonorthogonal, e_b_mini_vir, c_b_mini_vir, e_b_mini, c_b_mini = \
         self._compute_quambo(cb_occ, cb_vir, eps_b_occ, eps_b_vir, Fb, a_occ_b_set, a_vir_b_set, nbocc_mini, nbvir_mini, "BETA")
 
     
      # [4] Save
      self._quambo_a_nonorthogonal = quambo_a_nonorthogonal
      self._quambo_a_orthogonal    = quambo_a_orthogonal 
      self._quambo_b_nonorthogonal = quambo_b_nonorthogonal
      self._quambo_b_orthogonal    = quambo_b_orthogonal 

      self._c_a_mini_vir = c_a_mini_vir
      self._c_b_mini_vir = c_b_mini_vir

      self._e_a_mini_vir = e_a_mini_vir
      self._e_b_mini_vir = e_b_mini_vir

      self._c_a_mini = c_a_mini
      self._c_b_mini = c_b_mini

      self._e_a_mini = e_a_mini
      self._e_b_mini = e_b_mini

      # [5] restore Psi4 options state to original
      psi4.set_options({"reference":reference_stash})

      psi4.core.print_out(" @QUAMBO: Done.\n")


  # ---> protected interface <--- #

  def _compute_quambo(self, 
                            c_occ,      # Occupied MOs of the molecule
                            c_vir,      # Virtual MOs of the molecule
                            eps_occ,    # Occupied MO energies of the molecule
                            eps_vir,    # Virtual MO energies of the molecule
                            F,          # Fock matrix in non-orthogonal AO basis of the molecule
                            a_occ_set,  # The a_{nj}^* coefficients from Ref. [1] (projections of free-atom orbitals on c_occ)
                            a_vir_set,  # The a_{vj}^* coefficients from Ref. [1] (projections of free-atom orbitals on c_vir)
                            nocc_mini,  # Number of occupied orbitals in the minimal QUAMBO basis
                            nvir_mini,  # Number of virtual orbitals (VVOs) in the minimal QUAMBO basis
                            label):     # Label of the orbitals
      "Intrinsic routine to compute QUAMBOs from A* and < Occ(Vir) | A* > projections"

      # [3] Compute B matrix
      B = a_vir_set @ a_vir_set.T
     
      # [4] Diagonalize to find transformation matrix T 
      E, U = numpy.linalg.eigh(B)
      E = E[::-1]
      U = U[:,::-1]
      T = U[:,:nvir_mini].copy()
      Et= E[:nvir_mini].copy()
     
      # [5] Calculate R matrix 
      R = T @ T.T

      # [6] Calculate normalization
      Dj = (a_occ_set**2).sum(axis=0) + numpy.diag(a_vir_set.T @ R @ a_vir_set)
      Dj[numpy.where(Dj<0.0)] = 1.0
      Djm = 1./numpy.sqrt(Dj)
     
      # [7] Calculate non-orthogonal QUAMBOs 
      a_occ_set_nonorthogonal = numpy.einsum("j,nj->nj", Djm, a_occ_set)         # a_nj from Ref. [1]
      a_vir_set_nonorthogonal = numpy.einsum("j,nj->nj", Djm, R @ a_vir_set)     # a_vj from Ref. [1]
      quambo_nonorthogonal = c_occ @ a_occ_set_nonorthogonal 
      quambo_nonorthogonal+= c_vir @ a_vir_set_nonorthogonal

      # [8] Calculate orthogonal QUAMBOs
      S = a_occ_set_nonorthogonal.T @ a_occ_set_nonorthogonal 
      S+= a_vir_set_nonorthogonal.T @ a_vir_set_nonorthogonal
      Sm12 = matrix_power(S, -0.5)
      a_occ_set_orthogonal = a_occ_set_nonorthogonal @ Sm12                      # a_nj' from Ref. [1]
      a_vir_set_orthogonal = a_vir_set_nonorthogonal @ Sm12                      # a_vj' from Ref. [1]
      quambo_orthogonal = c_occ @ a_occ_set_orthogonal 
      quambo_orthogonal+= c_vir @ a_vir_set_orthogonal

      # [9] Diagonalize Fock matrix in QUAMBO basis
      quambo = quambo_orthogonal
      F_quambo = quambo.T @ F @ quambo
      e_quambo, c_quambo = numpy.linalg.eigh(F_quambo)

      # [9a] Test whether occupied orbital energies are correctly obtained
      e_quambo_occ = e_quambo[:nocc_mini]
      error = self._error(e_quambo_occ, eps_occ)
      assert (error < self._options.get_double("QUAMBO_EPS_THRESHOLD_CHECK")), \
                "Error in QUAMBO calculations! Error=%f" % error

      # [10] Compute virtual minimal MOs (VMO)
      e_quambo_vir = e_quambo[nocc_mini:]
      vir_mini_mo = c_quambo[:,nocc_mini:]
      vir_mini = quambo @ vir_mini_mo

      psi4.core.print_out("\n Eigenvalues of %6s Fock Operator in QUAMBO Basis\n\n" % label)
      psi4.core.print_out  (" Occupied\n")
      psi4.core.print_out  (" --------\n\n")
      psi4.core.print_out  ("          Full HF       QUAMBO        Error\n")
      for i in range(nocc_mini):
          psi4.core.print_out(" %3d %13.6f %13.6f %13.6f\n" % (i+1, eps_occ[i], e_quambo_occ[i], e_quambo_occ[i]-eps_occ[i]))
      psi4.core.print_out("\n Virtual\n")
      psi4.core.print_out  (" -------\n\n")
      psi4.core.print_out  ("          QUAMBO\n")
      for i in range(nvir_mini):
          psi4.core.print_out(" %3d %13.6f\n" % (i+1, e_quambo_vir[i]))
      psi4.core.print_out("\n")
      
      # [10a] Test for VMOs --> this test is wrong. vir_mini_test is not the same as vir_mini. VVOs are vir_mini, not vir_mini_test
      #vir_mini_test = c_vir @ T
      #error =  abs(vir_mini).sum() - abs(vir_mini_test).sum()
      #print(error)
      #print(vir_mini.T)
      #print(vir_mini_test.T)

      return quambo_orthogonal, quambo_nonorthogonal, e_quambo_vir, vir_mini, e_quambo, c_quambo

  def _atomize(self):
      "Create a list of psi molecule objects containing free atoms"
      atoms = []
      end = "symmetry c1\nunits bohr\nno_reorient\nno_com\n"

      if self.acbs:

         I = 0                                                                          
         for i in range(self._mol.natom()):
          coord_log = ""
          for j in range(self._mol.natom()):
           x = self._mol.x(j)
           y = self._mol.y(j)
           z = self._mol.z(j)
           if i==j:
             s = self._mol.symbol(j)
             try: m = self._unpe_atom[s] + 1
             except KeyError:
               raise KeyError("Please add antry in _atomize for atom %s and rerun" % s)
                                                                                        
             charge = 0
             multiplicity = m
                                                                                        
           else:
             s = "@"+self._mol.symbol(j)
           coord_log+= "%3s %14.8f %14.8f %14.8f\n" % (s, x, y, z)
                                                                                        
          log = "\n%i %i\n" % (charge, multiplicity)
          log+= coord_log
          log+= end 
                                                                                        
          atom = psi4.geometry(log)
          atoms.append(atom)
                                                                                        
          I += 1

      else:

         for i in range(self._mol.natom()):
             s = self._mol.symbol(i)
             x = self._mol.x(i)
             y = self._mol.y(i)
             z = self._mol.z(i)
             try: m = self._unpe_atom[s] + 1
             except KeyError:
               raise KeyError("Please add antry in _atomize for atom %s and rerun" % s)
    
             charge = 0
             multiplicity = m
             log = "\n%i %i\n" % (charge, multiplicity)
             log+= "%3s %14.8f %14.8f %14.8f\n" % (s, x, y, z) + end
    
             atom = psi4.geometry(log)
             atoms.append(atom)
         
      return atoms

  def _calculate_nbas_mini(self):
      # number of basis functions in minimal AO basis per atom
      nbf = 0
      for i in range(self._mol.natom()):
          atom = self._mol.symbol(i)
          try: nb = self._nbas_atom_mini[atom]
          except KeyError:
            raise KeyError("Please add antry in _calculate_nbas_mini for atom %s and rerun" % atom)
          nbf += nb
      return nbf

  def _error(self, a, b): return ((a-b)**2).sum()

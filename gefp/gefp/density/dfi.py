#!/usr/bin/python3
"""
 Demonstrates the use of Psi4 from Python level.
 Useful notes:
  o Use psi4.core module for most of the work
  o Useful modules within psi4.core: 
     - MintsHelper
     - Molecule
     - BasisSet
     - ExternalPotential
     others
  o Psi4 defines its own matrix type (psi4.core.Matrix). 
    Extracting numpy.array is easy:
       numpy_array = numpy.asarray(psi4_matrix)
    Creating Psi4 matrix from array is also easy:
       psi4_matrix = psi4.core.Matrix.from_array(numpy_array)     
  o To compute 1-el potential matrix for a set of charges 
    use ExternalPotential (charge positions are to be provided in Angstroms)
    unless charges are just nuclei within the basis set (in this case use of ao_potential method
    of MintsHelper is easier).  
  o ao_potential method of MintsHelper is limited only for nuclei within the same basis set
    (the nuclei are taken from the first basis set axis, for example:
      mints = MintsHelper(basis_X)         
      mints.ao_potential()                 -> nuclei taken from basis of mints object (basis_X)
      mints.ao_potential(basis_1, basis_2) -> nuclei taken from basis_1
  o Psi4 has efficient and easy to use method of defining fragments within a molecule (use '--' separator).
    Defining ghost atoms and extracting fragment i in the multimer-centred basis set is also very straighforward
    (method extract_subsets(...) of psi4.core.Molecule)
"""
__all__ = ['SCF', 'DFI']

import psi4
import numpy
import abc
import oepdev

MAX_NBF = 128


class SCF:
  """
 ---------------------------------------------------------------------------------------------------------------
                              Self-Consistent Field (SCF) Procedure for Hartree-Fock Model
 ---------------------------------------------------------------------------------------------------------------

 Demo for RHF-SCF method (closed shells). Implements SCF algorithm
 with primitive damping of the AO Fock matrix. 

 Usage:
  scf = SCF(molecule)
  scf.run(maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True)

 The above example runs SCF on 'molecule' psi4.core.Molecule object 
 starting from core Hamiltonian as guess (guess=None) 
 and convergence 1.0E-7 A.U. in total energy with 30 maximum iterations
 (10 of which are performed by damping of the Fock matrix with damping coefficient of 0.01).
 The SCF iterations are printed to standard output (verbose=True).
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, May 4th 2018
"""
  def __init__(self, mol):
      "Initialize BasisSet, Wavefunction and JK objects"
      # Basis set
      self._bfs = psi4.core.BasisSet.build(mol, "BASIS", psi4.core.get_global_option("BASIS"), puream=-1)
      # Wavefunction
      self._wfn = psi4.core.Wavefunction(mol, self._bfs)
      # Number of alpha electrons
      self._ndocc = self._wfn.nalpha()
      # Integral calculator
      self._mints = psi4.core.MintsHelper(self._bfs)
      # JK object
      self._jk = psi4.core.JK.build(self._bfs, jk_type="Direct")
      self._jk.set_memory(int(5e8))
      self._jk.initialize()
      ### Accessors
      # nuclear repulsion energy
      self.e_nuc = mol.nuclear_repulsion_energy()
      # Total Energy
      self.E = None
      # Density Matrix
      self.D = None
      # LCAO-MO coeffs (occ)
      self.Co= None
      # LCAO-MO coeffs (occ+vir)
      self.C = None
      # Fock matrix 
      self.F = None
      # Orbital energies
      self.eps = None
      # Hcore matrix
      self.H = None
      # Overlap integrals and orthogonalizer
      self.S = numpy.asarray( self._mints.ao_overlap() )
      self.X = self._orthogonalizer(self.S)
      return

  def run(self, maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, v_ext=None):
      "Solve SCF (public interface)"
      # Form Hcore                    
      T = self._mints.ao_kinetic()
      V = self._mints.ao_potential()
      H = T.clone()
      H.add(V)
      H = numpy.asarray(H)
      if v_ext is not None:
         H += numpy.asarray(v_ext)
      self.H = H.copy()
      # Getermine guess
      if guess is None:
         guess = H
      else:
         guess = numpy.asarray(guess)
      # Run SCF
      self._run(guess, maxit, conv, damp, ndamp, verbose)
      return

  # --- protected --- #

  def _run(self, guess, maxit, conv, damp, ndamp, verbose):
      "Solve SCF (protected interface)"
      # First step: Guess density matrix
      F    = numpy.dot(self.X, numpy.dot(guess, self.X)) 
      E, C = numpy.linalg.eigh(F)
      C    = numpy.dot(self.X, C)
      idx = numpy.argsort(E)
      E = E[  idx]
      C = C[:,idx]
      Co= C[:,:self._ndocc].copy()
      D = numpy.dot(Co,Co.T)
      
      niter = 0
      e_old = 1e8
      e_new = 1e7
      F_old = guess.copy()
      while (abs(e_old - e_new) > conv):
        niter += 1
        # form Fock matrix
        self._jk.C_clear()
        self._jk.C_left_add(psi4.core.Matrix.from_array(Co, "C matrix"))
        self._jk.compute()
        F_new = self.H + 2.0 * numpy.asarray(self._jk.J()[0]) - numpy.asarray(self._jk.K()[0])
        if niter < ndamp: 
           F = damp * F_old + (1.0 - damp) * F_new
        else:             
           F = F_new
        F_old = F.copy()
        # compute total energy
        e_old = e_new
        e_new = numpy.trace( numpy.dot(D, self.H + F) ) + self.e_nuc
        if verbose:
            print (" @SCF Iter {:02} E = {:14.8f}".format(niter, e_new))
        # transform Fock matrix to orthogonal AO basis           
        F = numpy.dot(self.X, numpy.dot(F, self.X))
        # diagonalize the Fock matrix
        E, C = numpy.linalg.eigh(F)
        # convert LCAO-MO coefficiets to non-orthogonal AO basis
        C = numpy.dot(self.X, C)
        # form density matrix
        idx = numpy.argsort(E) 
        E = E[  idx]
        C = C[:,idx]
        Co= C[:,:self._ndocc].copy()
        D = numpy.dot(Co,Co.T)
        # save
        self.D = D.copy()
        self.E = e_new
        self.F = F_old.copy()
        self.C = C.copy()
        self.Co= Co.copy()
        self.eps = E.copy()
        if niter > maxit: break
      return

  def _orthogonalizer(self, S):
      "Form orthogonalizer"
      L, U = numpy.linalg.eigh(S)
      L    = numpy.diag(1./numpy.sqrt(L))
      X    = numpy.dot(U, numpy.dot(L, U.T))
      return X


class DFI(abc.ABC):
  """
 ---------------------------------------------------------------------------------------------------------------
                                  Density Fragment Interaction (DFI) Method
 ---------------------------------------------------------------------------------------------------------------

 Demo for SCF-DFI method (closed shells).

 Usage:
  dfi = DFI(fragment_1, fragment_2, [...])  # OR: dfi = DFI(fragments)
  dfi.run(maxit=30, conv=1.0e-5, verbose_scf=False, conv_scf=1.0e-5, maxit_scf=100, damp_scf=0.14, ndamp_scf=0)

 Notes: 
  o fragmet_i is a psi4.core.Molecule wit one Psi4 Fragment
  o fragments is a psi4.core Molecule with multiple Psi4 Fragments ('--' separator in input)
  o SCF of unperturbed molecule is solved by Psi4, while the subsequent SCF's in DFI iterations are solved
    by SCF class instances of this Demo.
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, May 4th 2018
"""
  def __init__(self, *frags):
      "Initialize BasisSet and fragment lists"
      abc.ABC.__init__(self)
      self._j_only = False
      # --- Fragments as psi4.core.Molecule objects
      # Handle Psi4 Fragments within one Psi4 Molecule
      if len(frags) == 1: 
         nfrag = frags[0].nfragments() 
         assert nfrag > 1, " You must provide at least two fragments ('--' separator)!"
         self._nfrag = nfrag
         self._frags = []
         for i in range(frags[0].nfragments()):
             frags[0].fix_orientation(True)
             frags[0].print_out()
             f = frags[0].extract_subsets(i+1)
             f.print_out()
             f.reset_point_group('c1')
             f.fix_orientation(True)
             self._frags.append(f)
      # Handle separate Psi4 Molecules
      else:
         self._nfrag = len(frags)
         self._frags = frags
      #
      # --- Fragment lists                                                   Status every DFI iteration
      self._ens   = []  # Total energy of a fragment                         Updated
      self._ndocc = []  # number of doubly-occupied MO's                     Const
      self._nbfs  = []  # number of basis functions                          Const
      self._bfs   = []  # basis set                                          Const
      self._mints = []  # Integral calculator                                Running calculator
      self._jks   = []  # JK calculator                                      Running calculator
      self._Xs    = []  # AO orthogonalizer (S^-1/2)                         Const
      self._eps   = []  # Fock matrix eigenvalues                            Updated
      self._Cs    = []  # LCAO-MO coefficients (occupied)                    Updated
      self._Csall = []  # LCAO-MO coefficients (occ+vir)                     Updated
      self._H     = []  # Hcore matrix in original basis                     Const
      self._Hs    = []  # Hcore matrix in original basis (effective)         Updated
      self._Fs    = []  # Fock matrix in original basis                      Updated
      self._Ds    = []  # One-particle density matrix in original basis      Updated
      return

  @staticmethod
  def create(*frags, j_only=False):
      if j_only: return DFI_J (*frags)
      else:      return DFI_JK(*frags)

  def run(self, maxit=100, conv=1.0e-5, verbose_scf=False, conv_scf=1.0e-5, maxit_scf=100, damp_scf=0.14, ndamp_scf=0):
      """Runs DFI iterations"""
      self._run_init()
      self._run_iter(maxit, conv, conv_scf, maxit_scf, damp_scf, ndamp_scf, verbose_scf)
      return sum(self._ens) - self.en_0 

  # --- protected --- #
  def _run_init(self):
      "Solve unperturbed SCF equations"
      for frag in self._frags:
          frag.set_point_group(psi4.core.PointGroup('c1'))
          frag.reset_point_group('c1')
          frag.fix_orientation(True)
          frag.print_out()
          en, wfn = psi4.energy('scf', return_wfn=True, molecule=frag)
          bfs     = wfn.basisset()
          assert bfs.nbf() <= MAX_NBF, " Maximum number of basis functions %d exceeded!" % MAX_NBF
          mints   = psi4.core.MintsHelper(bfs)
          jk      = psi4.core.JK.build(bfs, jk_type="Direct")
          jk.set_memory(int(5e8)) # 4GB 
          jk.initialize()
          X       = self._orthogonalizer(numpy.asarray(wfn.S()))
          eps     = numpy.asarray(wfn.epsilon_a())
          C       = numpy.asarray(wfn.Ca())
          Co      = numpy.asarray(wfn.Ca_subset("AO","OCC"))
          H       = numpy.asarray(wfn.H())
          F       = numpy.asarray(wfn.Fa())
          D       = numpy.asarray(wfn.Da())
          #
          self._ens  .append(en)             # Total energy of a fragment
          self._ndocc.append(wfn.nalpha())   # number of doubly-occupied MO's
          self._nbfs .append(bfs.nbf())      # number of basis functions
          self._bfs  .append(bfs)            # basis set
          self._mints.append(mints)          # Integral calculator
          self._jks  .append(jk)             # JK calculator
          self._Xs   .append(X)              # AO orthogonalizer (S^-1/2)
          self._eps  .append(eps)            # Eigenvalues of Fock matrix
          self._Cs   .append(Co)             # LCAO-MO coefficients (occ)
          self._Csall.append(C)              # LCAO-MO coefficients (occ+vir)
          self._H    .append(H)              # Hcore matrix in original basis
          self._Hs   .append(H)              # Effective Hcore matrix in original basis
          self._Fs   .append(F)              # Fock matrix in original basis
          self._Ds   .append(D)              # One-particle density matrix in original basis
      # compute interfragment nuclear repulsion energy
      self.enuc = self._interfragmentNuclearRepulsionEnergy()
      # 
      self.en_0 = None
      return

  def _run_iter(self, maxit, conv, conv_scf, maxit_scf, damp_scf, ndamp_scf, verbose_scf):
      "DFI iterations (protected interface)"
      niter = 0
      sum_new = sum(self._ens)
      sum_old = 1e8
      self._print_iter(0)
      self.en_0 = sum(self._ens)
      while (abs(sum_new - sum_old) > conv):
          niter += 1
          for i in range(self._nfrag): self._update_frag(i, conv_scf, maxit_scf, damp_scf, ndamp_scf, verbose_scf)
          self._print_iter(niter)
          if niter > maxit: break
          sum_old = sum_new 
          sum_new = sum(self._ens)
      self._finalize()
      return

  def _print_iter(self, n):
      "Prints the DFI iterates"
      print (" @DFI Iter {:02} E = {:14.8f}".format(n, sum(self._ens)))
      log = self._nfrag * "%15.6f" % tuple(self._ens)
      print (log)
      return

  def _finalize(self):
      "Finalize the DFI process"
      e_int = sum(self._ens) - self.en_0 #+ self.enuc
      print (" DFI interaction energy = {:14.6f}".format(e_int))
      return

  def _update_frag(self, i, conv_scf, maxit_scf, damp_scf, ndamp_scf, verbose_scf):
      "Update the i-th fragment"
      # compute effective potential matrix of i-th fragment in the presence of other fragments 
      V = numpy.zeros(self._H[i].shape)
      for j in range(self._nfrag):
          if j!=i:
             V += self._eval_vnuc(i, j)
             V += self._eval_vel (i, j)
      # solve SCF for i-th fragment
      scf = SCF(self._frags[i])
      scf.run(guess=self._Fs[i].copy(), v_ext=V, conv=conv_scf, maxit=maxit_scf, damp=damp_scf, ndamp=ndamp_scf, verbose=verbose_scf)
      # save
      self._Hs[i]  = self._H[i] + V
      self._ens[i] = scf.E
      self._Fs [i] = scf.F.copy()
      self._Cs [i] = scf.Co.copy()
      self._Csall[i]= scf.C.copy()
      self._Ds [i] = scf.D.copy()
      self._eps[i] = scf.eps.copy()
      return

  def _eval_vnuc(self, i, j):
      "Nuclear contribution from j-th fragment to i-th fragment"
      ep = psi4.core.ExternalPotential()
      BohrToAngstrom = 0.5291772086
      for a in range(self._frags[j].natom()):
          q = self._frags[j].Z(a)
          x = self._frags[j].x(a) * BohrToAngstrom
          y = self._frags[j].y(a) * BohrToAngstrom
          z = self._frags[j].z(a) * BohrToAngstrom
          ep.addCharge(q, x, y, z)
      V = numpy.asarray(ep.computePotentialMatrix(self._bfs[i]))
      return V

  def _eval_vel(self, i, j):
      "Electronic contribution from j-th fragment to i-th fragment"
      ## OLD code: memory consuming
      #eri_J_ij = numpy.asarray(\
      #              self._mints[i].ao_eri(self._bfs[j], self._bfs[j],
      #                                    self._bfs[i], self._bfs[i]))
      #ni = self._nbfs[i]
      #nj = self._nbfs[j]
      #J = numpy.dot(self._Ds[j].ravel(), eri_J_ij.reshape(nj*nj,ni*ni))
      #V = 2.0*J
      #if not self._j_only:
      #   eri_K_ij = numpy.asarray(\
      #                 self._mints[i].ao_eri(self._bfs[j], self._bfs[i],
      #                                       self._bfs[j], self._bfs[i]))
      #                                                                                        
      #   K = numpy.dot(self._Ds[j].ravel(), eri_K_ij.transpose(0,2,1,3).reshape(nj*nj,ni*ni))
      #   V-= K
      #return V.reshape(ni, ni)
      f_aabb = psi4.core.IntegralFactory(self._bfs[i], self._bfs[i], self._bfs[j], self._bfs[j])
      D      = psi4.core.Matrix.from_array(self._Ds[j])
      if not self._j_only:
         f_abab = psi4.core.IntegralFactory(self._bfs[i], self._bfs[j], self._bfs[i], self._bfs[j])
         V      = oepdev.calculate_DFI_Vel_JK(f_aabb, f_abab, D)
      else:
         V      = oepdev.calculate_DFI_Vel_J(f_aabb, D)
      return V.to_array(dense=True)

  def _orthogonalizer(self, S):
      "Compute the orthogonalizer matrix X"
      L, U = numpy.linalg.eigh(S)
      L    = numpy.diag(1./numpy.sqrt(L))
      X    = numpy.dot(U, numpy.dot(L, U.T))
      return X

  def _interfragmentNuclearRepulsionEnergy(self):
      "Compute the nuclear repulsion between all fragments"
      E = 0.0
      for ni in range(self._nfrag):
          for nj in range(self._nfrag):
              if ni!=nj:
                 for i in range(self._frags[ni].natom()):
                     Z_i = self._frags[ni].Z(i)
                     x_i = self._frags[ni].x(i)
                     y_i = self._frags[ni].y(i)
                     z_i = self._frags[ni].z(i)
                     for j in range(self._frags[nj].natom()):
                         Z_j = self._frags[nj].Z(j)
                         x_j = self._frags[nj].x(j)
                         y_j = self._frags[nj].y(j)
                         z_j = self._frags[nj].z(j)
                         r2_ij = (x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2
                         E += Z_i * Z_j / numpy.sqrt(r2_ij)
      return E

class DFI_J(DFI):
  def __init__(self, *frags):
      DFI.__init__(self, *frags)
      self._j_only = True

class DFI_JK(DFI):
  def __init__(self, *frags):
      DFI.__init__(self, *frags)
      self._j_only = False

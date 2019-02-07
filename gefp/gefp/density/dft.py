#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DFT module.
 Bartosz Błasiak, Gundelfingen, Feb 2019
"""

import sys
import math
import numpy
import numpy.linalg
import psi4
import oepdev

__all__ = ["DFT"]


class DFT:
  """
 ---------------------------------------------------------------------------------------------------------------
 ---------------------------------------------------------------------------------------------------------------

 Demo for DFT method (closed shells). Implements DFT density matrix projected gradient algorithm (TODO)

 Usage:
  dft = DFT(wfn)
  e = dft.static_energy(dmft='MBB')
  dft.run(maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, V_ext=None)

 The above example runs DFT on 'molecule' psi4.core.Molecule object 
 starting from core Hamiltonian as guess (guess=None) 
 and convergence 1.0E-7 A.U. in total energy with 30 maximum iterations
 (10 of which are performed by damping of the Fock matrix with damping coefficient of 0.01).
 The projected gradient iterations are printed to standard output (verbose=True).
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, May 4th 2018
"""
  def __init__(self, wfn):
      "Initialize BasisSet, Wavefunction and JK objects"
      # Wavefunction
      self._wfn = wfn#.c1_deep_copy(wfn.basisset())
      # Molecule
      self._mol = self._wfn.molecule()
      # Basis set
      self._bfs = self._wfn.basisset()
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
      self.e_nuc = self._mol.nuclear_repulsion_energy()
      # Total Energy
      self.E = self._wfn.energy()
      # Density Matrix
      self.D = self._wfn.Da()
      # LCAO-MO coeffs
      self.C = None #self._wfn.Ca()
      # Fock matrix 
      self.F = None #self._wfn.Fa()
      # Hcore matrix
      self.H = self._wfn.H()
      # Overlap integrals and orthogonalizer
      self.S = numpy.asarray( self._mints.ao_overlap() )
      self.X = self._orthogonalizer(self.S)
      # External potential
      self.V_ext = None
      return

  def static_energy(self, dmft='MBB'):
      """\
 Compute static total energy of the system based on the current 1-particle density matrix
 and a chosen DMFT functional. Available functionals:

   o 'HF'  - relaxed HF functional
   o 'CHF' - corrected HF functional
   o 'MBB' - Muller-Buijse-Baerends functional (default)
   o 'GU'  - Goedecker-Urmigar functional
   o 'BB1' - Blasiak 1 functional
   o 'BB2' - Blasiak 2 functional

"""
      # Natural Orbitals from current OPDM
      D = self._wfn.Da().to_array(dense=True)
      S = self._wfn.S ().to_array(dense=True)
      n, c = self.natural_orbitals(D, orthogonalize_first=S, order='descending', no_cutoff=0.0, renormalize=True)
      # Nuclear Repulsion
      E_N = self.e_nuc
      # 1-electron energy
      E_1 = self._compute_hcore_energy(n, c)
      # Hartree 2-electron energy
      E_H = self._compute_hartree_energy(n, c)
      # Exchange-Correlation energy
      E_XC = self._compute_XC_energy(n, c, dmft)
      # Total energy
      E = E_N + E_1 + E_H + E_XC
      return E

  def natural_orbitals(self, D, orthogonalize_first=None, order='descending', original_ao_mo=True, renormalize=False, no_cutoff=False):
      "Compute the Natural Orbitals from a given ODPM"
      if orthogonalize_first is not None:
         S = orthogonalize_first
         D_ = self._orthogonalize_OPDM(D, S)
      else:
         D_ = D
      n, U = numpy.linalg.eigh(D_)
      n[numpy.where(n<0.0)] = 0.0
      if original_ao_mo:
         assert orthogonalize_first is not None
         U = numpy.dot(self._orthogonalizer(orthogonalize_first), U)
      if no_cutoff is False: no_cutoff = self.no_cutoff
      if no_cutoff != 0.0:
         ids = numpy.where(n>=self.no_cutoff)
         n = n[ids]
         U =(U.T[ids]).T
      if order=='ascending': 
         pass
      elif order=='descending':
         n = n[  ::-1]
         U = U[:,::-1]
      else: raise ValueError("Incorrect order of NO orbitals. Possible only ascending or descending.")
      if renormalize is True:
         d = numpy.round(n.sum()) - n.sum()
         d/= numpy.float64(n.size)
         n+= d
         n[numpy.where(n<0.0)] = 0.0
      return n, U

  def _orthogonalize_OPDM(self, D, S):
      "Transforms the one-particle density matrix to orthogonal space"
      Y = self._deorthogonalizer(S)
      return numpy.dot(Y, numpy.dot(D, Y.T))

  def _deorthogonalizer(self, S):
      "Compute the deorthogonalizer matrix from the overlap matrix"
      s, u = numpy.linalg.eig(S)
      s = numpy.sqrt(s)
      Y = numpy.dot(u, numpy.dot(numpy.diag(s), u.T))
      return Y

  def _orthogonalizer(self, S):
      "Compute the orthogonalizer matrix from the overlap matrix"
      s, u = numpy.linalg.eig(S)
      sm = 1.0/numpy.sqrt(s)
      X = numpy.dot(u, numpy.dot(numpy.diag(sm), u.T))
      return X



  def _compute_hcore_energy(self, n, c):
      raise NotImplementedError
  def _compute_hartree_energy(self, n, c):
      raise NotImplementedError
  def _compute_XC_energy(self, dmft, n, c):
      raise NotImplementedError

  def run(self, maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, V_ext=None):
      "Solve SCF (public interface)"

      raise NotImplementedError

      # prepare H_core
      if guess is None:
         # Form Hcore                    
         T = self._mints.ao_kinetic()
         V = self._mints.ao_potential()
         guess = T.clone()
         guess.add(V)
         guess = numpy.asarray(guess)
      else: guess = numpy.asarray(guess)

      # external one-electron potential
      self.V_ext = V_ext
      if self.V_ext is not None:
         guess  += self.V_ext
         self.H += self.V_ext

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
      C = C[:,idx]
      C = C[:,:self._ndocc]
      D = numpy.dot(C,C.T)
      
      niter = 0
      e_old = 1e8
      e_new = 1e7
      F_old = guess.copy()
      while (abs(e_old - e_new) > conv):
        niter += 1
        # form Fock matrix
        self._jk.C_clear()
        self._jk.C_left_add(psi4.core.Matrix.from_array(C, "C matrix"))
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
        if verbose: print(" @SCF Iter %02d. E = %14.8f" % (niter, e_new))
        # transform Fock matrix to orthogonal AO basis           
        F = numpy.dot(self.X, numpy.dot(F, self.X))
        # diagonalize the Fock matrix
        E, C = numpy.linalg.eigh(F)
        # convert LCAO-MO coefficiets to non-orthogonal AO basis
        C = numpy.dot(self.X, C)
        # form density matrix
        idx = numpy.argsort(E) 
        C = C[:,idx]
        C = C[:,:self._ndocc]
        D = numpy.dot(C,C.T)
        # save
        self.D = D.copy()
        self.E = e_new
        self.F = F_old.copy()
        self.C = C.copy()
        if niter > maxit: break
      return

  def _orthogonalizer(self, S):
      "Form orthogonalizer"
      L, U = numpy.linalg.eigh(S)
      L    = numpy.diag(1./numpy.sqrt(L))
      X    = numpy.dot(U, numpy.dot(L, U.T))
      return X

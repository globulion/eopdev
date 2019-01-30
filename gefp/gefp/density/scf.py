#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 SCF module.
 Bartosz Błasiak, Gundelfingen, Jan 2019
"""

import sys
import math
import numpy
import numpy.linalg
import psi4
import oepdev

__all__ = ["SCF"]


class SCF:
  """
 ---------------------------------------------------------------------------------------------------------------
                              Self-Consistent Field (SCF) Procedure for Hartree-Fock Model
 ---------------------------------------------------------------------------------------------------------------

 Demo for RHF-SCF method (closed shells). Implements SCF algorithm
 with primitive damping of the AO Fock matrix. 

 Usage:
  scf = SCF(wfn)
  scf.run(maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, V_ext=None)

 The above example runs SCF on 'molecule' psi4.core.Molecule object 
 starting from core Hamiltonian as guess (guess=None) 
 and convergence 1.0E-7 A.U. in total energy with 30 maximum iterations
 (10 of which are performed by damping of the Fock matrix with damping coefficient of 0.01).
 The SCF iterations are printed to standard output (verbose=True).
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, May 4th 2018
"""
  def __init__(self, wfn):
      "Initialize BasisSet, Wavefunction and JK objects"
      # Wavefunction
      self._wfn = wfn.c1_deep_copy(wfn.basisset())
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
      self.C = self._wfn.Ca()
      # Fock matrix 
      self.F = self._wfn.Fa()
      # Hcore matrix
      self.H = self._wfn.H()
      # Overlap integrals and orthogonalizer
      self.S = numpy.asarray( self._mints.ao_overlap() )
      self.X = self._orthogonalizer(self.S)
      # External potential
      self.V_ext = None
      return

  def run(self, maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, V_ext=None):
      "Solve SCF (public interface)"
      if guess is None:
         # Form Hcore                    
         T = self._mints.ao_kinetic()
         V = self._mints.ao_potential()
         guess = T.clone()
         guess.add(V)
         guess = numpy.asarray(guess)
      else: guess = numpy.asarray(guess)

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

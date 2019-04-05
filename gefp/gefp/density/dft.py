#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT module.
 Bartosz Błasiak, Gundelfingen, Feb 2019
"""

import os
import sys
import math
import numpy
import numpy.linalg
import scipy.optimize
import psi4
import oepdev
from . import functional

__all__ = ["DMFT"]


class DMFT:
  """
 ---------------------------------------------------------------------------------------------------------------
                                              The DMFT method
 ---------------------------------------------------------------------------------------------------------------

 Demo for DMFT method (closed shells). Implements DMFT density matrix projected gradient algorithm (TODO)

 Usage:
  dmft = DMFT(wfn)
  e = dmft.static_energy(dmft='MBB')
  dmft.run(maxit=30, conv=1.0e-7, guess=None, damp=0.01, ndamp=10, verbose=True, V_ext=None)

 The above example runs DMFT on 'molecule' psi4.core.Molecule object 
 starting from core Hamiltonian as guess (guess=None) 
 and convergence 1.0E-7 A.U. in total energy with 30 maximum iterations
 (10 of which are performed by damping of the Fock matrix with damping coefficient of 0.01).
 The projected gradient iterations are printed to standard output (verbose=True).
 ---------------------------------------------------------------------------------------------------------------
                                                                       Last Revision: Gundelfingen, Feb 7th 2019
"""
  def __init__(self, wfn, guess='current', V_ext=None):
      "Initialize BasisSet, Wavefunction and JK objects"
      # Wavefunction
      self._wfn = wfn#.c1_deep_copy(wfn.basisset())
      # Molecule
      self._mol = self._wfn.molecule()
      # Basis set
      self._bfs = self._wfn.basisset()
      # Integral calculator
      self._mints = psi4.core.MintsHelper(self._bfs)
      # JK object
      self._jk = psi4.core.JK.build(self._bfs, jk_type="Direct")
      self._jk.set_memory(int(5e8))
      self._jk.initialize()
      # Nuclear repulsion energy
      self.e_nuc = self._mol.nuclear_repulsion_energy()
      # Number of electron pairs
      self.Np = self._wfn.nalpha()

      ### Constant AO matrices
      # Overlap integrals and orthogonalizer
      self.S = self._wfn.S().to_array(dense=True) #numpy.asarray( self._mints.ao_overlap() )
      self.X = self._orthogonalizer(self.S)
      self.Y = numpy.linalg.inv(self.X)
      # Hcore matrix
      self.H = self._wfn.H().to_array(dense=True)
      #self.H = numpy.linalg.multi_dot([self.X, self.H, self.X])
      # External potential
      self.V_ext = V_ext
      if V_ext is not None and guess == 'current':
          raise ValueError(" External potential can only be set for 'HCore' guess!")
      #if V_ext is not None:
      #   self.V_ext = numpy.linalg.multi_dot([self.X, self.V_ext, self.X])

      ### Current OPDM and total energy
      if guess == 'current':  # Guess based on current density matrix from the input wavefunction
         assert self.V_ext is None
         self.E = self._wfn.energy()
         self.D = self._wfn.Da().to_array(dense=True)
         #self.D = self._orthogonalize_OPDM(self.D, self.S)
      elif guess == 'hcore':  # Guess based on one-electron Hamiltonian
         if self.V_ext is not None: self.H += self.V_ext
         e, c = numpy.linalg.eigh(self.H)
         c = c[:,::-1]
         D = numpy.zeros((self._bfs.nbf(), self._bfs.nbf()), numpy.float64)
         for i in range(self.Np):
             D += numpy.outer(c[:,i], c[:,i])
         self.D = D
         self.E = 2.0*e.sum() + self.e_nuc
      else:
         raise ValueError("Only 'Current' or 'HCore' ODPM's are supported as starting points.")
      # Natural Orbitals
      self.n, self.c = self.natural_orbitals(self.D, orthogonalize_first=self.S, 
                            order='descending', no_cutoff=0.0, renormalize=False, original_ao_mo=True)
      return

  def run(self, dmft, algorithm='sd-par',
              maxit=30, conv=1.0e-7, g_0=0.0001, step=0.000001, verbose=True, **kwargs):
      "Solve SCF-DMFT equations (public interface)"
      if algorithm.lower() == 'sd-par':
         success = self._run(dmft, maxit, conv, g_0, step, verbose, self.n, self.c, **kwargs)
      elif algorithm.lower() == 'sd-proj-r':
         success = self._run_proj_R(dmft, maxit, conv, g_0, step, verbose, self.n, self.c, **kwargs)
      elif algorithm.lower() == 'sd-proj-t':
         success = self._run_proj_T(dmft, maxit, conv, g_0, step, verbose, self.n, self.c, **kwargs)
      else: raise ValueError(" Unknown algorithm %s" % algorithm)
      return success

  def static_energy(self, dmft='MBB', **kwargs):
      """\
 Compute static total energy of the system based on the current 1-particle density matrix
 and a chosen exchange-correlation density functional. Available functionals:

   o 'HF'  - relaxed HF functional
   o 'CHF' - corrected HF functional
   o 'MBB' - Muller-Buijse-Baerends functional (default)
   o 'GU'  - Goedecker-Urmigar functional
   o 'BB1' - Blasiak 1 functional
   o 'BB2' - Blasiak 2 functional (requires GAMMA)
   o 'BBT' - Blasiak: Taylor Series of Unknown Generator functional

"""
      # Nuclear Repulsion
      E_N = self.e_nuc
      # 1-electron energy
      E_1 = self._compute_hcore_energy()
      # Hartree 2-electron energy
      E_H = self._compute_hartree_energy()
      # Exchange-Correlation energy
      E_XC = self._compute_XC_energy(dmft, **kwargs)
      # Total energy
      E = E_N + E_1 + E_H + E_XC
      return E

  # --- Public interface (expert) --- #

  def fij(self, n, dmft, **kwargs):
      "Compute Exchange-Correlation Coefficients"
      if   dmft.lower() == 'mbb'           : f = self._fij_mbb(n)
      elif dmft.lower() == 'mbb_pc'        : f = self._fij_mbb_PC(n)  # same as bbc1
      elif dmft.lower() == 'mbb0'          : f = self._fij_mbb_0(n)
      elif dmft.lower() == 'hf'            : f = self._fij_hf (n)
      elif dmft.lower() == 'chf'           : f = self._fij_chf(n)
      elif dmft.lower() == 'cga'           : f = self._fij_cga(n)
      elif dmft.lower() == 'gu'            : f = self._fij_gu (n)
      elif dmft.lower() == 'bbc1'          : f = self._fij_bbc1(n)
      elif dmft.lower() == 'bbc2'          : f = self._fij_bbc2(n)
      elif dmft.lower() == 'bb1'           : f = self._fij_bb1(n)
      elif dmft.lower() == 'bb2'           : f = self._fij_bb2(n, kwargs["gamma"])
      elif dmft.lower() == 'bb3'           : f = self._fij_bb3(n)
      elif dmft.lower() == 'bbp'           : f = self._fij_bbp(n)
      elif dmft.lower() == 'bbb'           : f = self._fij_bbb(n, kwargs["kmax"], kwargs["scale"])
      elif dmft.lower() == 'bby'           : f = self._fij_bby(n, kwargs["coeffs"], kwargs["kmax"])
      elif dmft.lower() == 'bbz'           : f = self._fij_bbz(n, kwargs["coeffs"], kwargs["kmax"], kwargs["pc"])
      elif dmft.lower() == 'bbi'           : f = self._fij_bbZ(n, kwargs["kmax"], kwargs["pc"], kwargs["ao"])
      elif dmft.lower() == 'bb_opt'        : f = self._fij_bbopt(n, kwargs["kmax"], kwargs["coeffs"])
      elif dmft.lower() == 'bb_opt_2'      : f = self._fij_bbopt_2(n, kwargs["kmax"], kwargs["coeffs"])
      elif dmft.lower() == 'bb_pade_p'     : f = self._fij_bbpade(n, kwargs["kmax"], coeffs=[4.19679,-164.543,5.92545,94.7776,0.0,-0.475193,-31.6913,-226.502,17.1221,91.8458,0.0,-1.81797,-49.7973])
     #elif dmft.lower() == 'bb_opt_2_p'    : f = self._fij_bbopt_2(n, kwargs["kmax"], coeffs=[-0.14881973,5.16901966,-23.44068276,-0.66111539,1.32165395])
      elif dmft.lower() == 'bb_opt_2_p'    : f = self._fij_bbopt_2(n, kwargs["kmax"], coeffs=[ 1.14633566,7.31708483,-23.28681062,-0.80297864,1.27633697])
      elif dmft.lower() == 'bbq2'          : f = self._fij_bbq2(n)
      elif dmft.lower() == 'bbh'           : f = self._fij_bbh(n)
      elif dmft.lower() == 'bbh2'          : f = self._fij_bbh2(n)
      elif dmft.lower() == 'pow'           : f = self._fij_power(n, kwargs["pow_exp"])
      elif dmft.lower() == 'bbt'           : f = self._fij_bbt(n, kwargs["bbt_coeffs"], kwargs["bbt_kmax"], kwargs["bbt_hf"], kwargs["bbt_log"])
      elif dmft.lower() == 'bbg'           : f = self._fij_bbg(n, kwargs["coeffs"], kwargs["kmax"])
      elif dmft.lower() == 'bbx'           : f = self._fij_BBX(n, kwargs["coeffs"], kwargs["kmax"])
      else: raise ValueError("The functional %s is not available. Misspelling error?" % dmft.lower())
      return f

  def fij_1(self, n, m, dmft, step=None, **kwargs):
      "Compute first derivatives of EX coefficients wrt the m-th NO occupancy"
      if   dmft.lower() == 'hf'       : f_m = self._fij_1_hf(n, m)
      elif dmft.lower() == 'mbb'      : f_m = self._fij_1_numerical(n, m, step, self._fij_mbb, **kwargs)
      elif dmft.lower() == 'bb_pade'  : f_m = self._fij_1_numerical(n, m, step, self._fij_bbpade, **kwargs)
      elif dmft.lower() == 'gu'       : f_m = self._fij_1_numerical(n, m, step, self._fij_gu)
      elif dmft.lower() == 'chf'      : f_m = self._fij_1_numerical(n, m, step, self._fij_chf)
      else: raise NotImplementedError
      return f_m

  def Kij(self, c):
      "Exchange (ij|ji) integral matrix"
      C = psi4.core.Matrix.from_array(c, "Natural Orbitals LCAO matrix")
      #eri_K_ij = oepdev.calculate_Kij(self._wfn.c1_deep_copy(self._wfn.basisset()), C).to_array(dense=True)
      eri_K_ij = oepdev.calculate_JK(self._wfn, C)[1].to_array(dense=True)
      psi4.core.clean()

      #bfs = self._wfn.basisset()
      #mints = psi4.core.MintsHelper(bfs)
      #eri = numpy.asarray(mints.ao_eri(bfs, bfs,
      #                                 bfs, bfs))
      #eri_K_ij = numpy.einsum("ijkl,ia,jb,kb,la->ab", eri, c, c, c, c) 
      return eri_K_ij



  # ===> protected <=== #

  # ----- Steepest-Descents Algorithm ----- #

  def _run(self, dmft, maxit, conv, g_0, step, verbose, n, C, **kwargs):
      "Run DMFT-SCF with the Steepest-Descents minimization and density matrix projection algorithm (protected interface)"

      nn = len(n)
      iteration = 0
      success = False

      # Exchange-Correlation coefficients
      fij = self.fij(n, dmft, **kwargs)

      # Starting energy
      self.E = self.static_energy(dmft, **kwargs)
      E_old = self.E
      if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_old))

      # First iteration
      iteration += 1
      C_ = numpy.dot(self.Y, C)
      x_old_2 = numpy.hstack([n, C_.ravel()])
      gradient_2 = self._gradient(n, fij, C, dmft, step, **kwargs)
      x_old_1 = x_old_2 - g_0 * gradient_2

      n_old_1 = x_old_1[:nn]
      C_old_1_= x_old_1[nn:].reshape(nn,nn) ; C_old_1 = numpy.dot(self.X, C_old_1_)
      n_old_1, C_old_1 = self._density_matrix_projection(n_old_1, C_old_1)
      self.D = numpy.linalg.multi_dot([C_old_1, numpy.diag(n_old_1), C_old_1.T])
      self.n = n_old_1.copy()
      self.c = C_old_1.copy()
      E_new = self.static_energy(dmft, **kwargs)
      if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_new))

      # Prepare for further iterations
      fij_old_1 = self.fij(n_old_1, dmft, **kwargs)

      n_old_2 = n.copy()
      C_old_2 = C.copy()
      fij_old_2 = fij.copy()

      # Continue iterations with changing step of SD
      stop = False
      iteration += 1
      while stop is False:
          # New density matrix and XC coefficients
          n_new, C_new = self._st_step(n_old_1, fij_old_1, C_old_1, n_old_2, fij_old_2, C_old_2, dmft, step, **kwargs)
          n_new, C_new = self._density_matrix_projection(n_new, C_new)
          fij_new = self.fij(n_new, dmft, **kwargs)
          self.D = numpy.linalg.multi_dot([C_new, numpy.diag(n_new), C_new.T])
          self.n = n_new.copy()
          self.c = C_new.copy()

          # compute current energy
          E_new = self.static_energy(dmft, **kwargs)
          if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_new))

          # check if converged?
          if abs(E_new-E_old) < conv: 
             stop = True
             success = True

          # check if max iterations were exceeded
          if iteration >= maxit: 
             stop = True
             success = False

          # prepare for next iteration
          iteration += 1
          E_old = E_new
          self.E = E_old

          n_old_2 = n_old_1.copy()
          fij_old_2 = fij_old_1.copy()
          C_old_2 = C_old_1.copy()

          n_old_1 = n_new.copy()
          fij_old_1 = fij_new.copy()
          C_old_1 = C_new.copy()

      if verbose and success:
         print(" DMFT-SCF iterations converged.")
         print(" Final Energy = %14.8f" % self.E)
      if verbose and not success:
         print(" DMFT-SCF iterations did not converge.")
      return success

  def _grad_n(self, n, fij, C, dmft, step, **kwargs):
      "Energy gradient wrt NO occupation numbers"
      nn = len(n)

      # 1-electron contributon
      Hmm = numpy.linalg.multi_dot([C.T, self.H, C]).diagonal()
      grad = 2.0 * Hmm

      # 2-electron contribution
      D = numpy.dot(C, numpy.dot(2.0*numpy.diag(n), C.T))
      I = numpy.identity(nn, numpy.float64)

      # J-type
      self._jk.C_clear()
      #self._jk.set_do_J(True)
      #self._jk.set_do_K(False)
      self._jk.C_left_add(psi4.core.Matrix.from_array(self.D, ""))
      self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      self._jk.compute()
      J = self._jk.J()[0].to_array(dense=True)
      grad += 4.0 * numpy.linalg.multi_dot([C.T, J, C]).diagonal()

      # K-type
      Kij = self.Kij(self.c)
      for m in range(nn):
          fij_m = self.fij_1(n, m, dmft, step, **kwargs)
          grad[m] -=(numpy.dot(Kij, fij_m)).trace()

      #self._jk.set_do_K(True)
      return grad

  def _grad_C(self, n, fij, C):
      "Energy gradient wrt LCAO-NO wavefunction coefficients"
      nn = len(n)

      # 1-electron contributon
      Ham = numpy.dot(self.H, C)
      grad = Ham * n

      # 2-electron contribution
      I = numpy.identity(nn, numpy.float64)

      # J-type
      self._jk.C_clear()                                           
      #self._jk.set_do_J(True)
      #self._jk.set_do_K(False)

      self._jk.C_left_add(psi4.core.Matrix.from_array(self.D, ""))
      self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      self._jk.compute()
      J = self._jk.J()[0].to_array(dense=True)
      J = 2.0 * numpy.dot(J, C)
      for m in range(nn): grad[:,m] += n[m] * J[:,m]

      #for m in range(nn):
      #    Am = self._A(n, C, m)
      #    self._jk.C_left_add(psi4.core.Matrix.from_array(Am, ""))
      #    self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      #self._jk.compute()

      #for m in range(nn):
      #    Jm = self._jk.J()[m].to_array(dense=True)
      #    grad[:,m] += numpy.dot(Jm, C[:,m])

      # K-type
      self._jk.C_clear()                                           
      #self._jk.set_do_J(False)
      #self._jk.set_do_K(True)

      for m in range(nn):
          Bm = self._B(fij, C, m)
          self._jk.C_left_add(psi4.core.Matrix.from_array(Bm, ""))
          self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      self._jk.compute()

      for m in range(nn):
          Km = self._jk.K()[m].to_array(dense=True)
          grad[:,m] -= numpy.dot(Km, C[:,m])

      # finalize
      grad *= 4.0
      #self._jk.set_do_J(True)
      #self._jk.set_do_K(True)

      # transform gradient to orthogonal AO basis
      grad = numpy.dot(self.X, grad)
      return grad

  def _A(self, n, C, m):
      Am = 2.0 * self.D * n[m]
      return Am
  def _B(self, fij, C, m):
      Bm = numpy.linalg.multi_dot([C, numpy.diag(fij[:,m]), C.T])
      return Bm

  def _gradient(self, n, fij, C, dmft, step, **kwargs):
      "Gradient vector in the parameter space of n x C.ravel(). Returns gradient in orthogonal AO basis"
      dE_n = self._grad_n(n, fij, C, dmft, step, **kwargs)
      dE_C = self._grad_C(n, fij, C)

      gradient = numpy.hstack([dE_n, dE_C.ravel()])
      return gradient

  def _st_step(self, n_old_1, fij_old_1, C_old_1, 
                     n_old_2, fij_old_2, C_old_2,
                     dmft, step, **kwargs):
      "Next Steepest-Descents Step"
      nn = len(n_old_1)

      # transform C to orthogonal AO basis
      C_old_1_ = numpy.dot(self.Y, C_old_1)
      C_old_2_ = numpy.dot(self.Y, C_old_2)

      x_old_1= numpy.hstack([n_old_1, C_old_1_.ravel()])
      x_old_2= numpy.hstack([n_old_2, C_old_2_.ravel()])

      gradient_1 = self._gradient(n_old_1, fij_old_1, C_old_1, dmft, step, **kwargs)
      gradient_2 = self._gradient(n_old_2, fij_old_2, C_old_2, dmft, step, **kwargs)

      norm = numpy.linalg.norm(gradient_1 - gradient_2)
      g = numpy.dot(x_old_1 - x_old_2, gradient_1 - gradient_2) / norm**2
      #print(" G = ", g)
      g = abs(g)
      x_new = x_old_1 - g * gradient_1
      n_new = x_new[:nn]
      C_new_= x_new[nn:].reshape(nn,nn)

      # transform C back to original AO basis
      C_new = numpy.dot(self.X, C_new_)
      return n_new, C_new

  # ----- Density Projection Steepest Descents ----- #

  def _make_P(self, n, C):
      return numpy.linalg.multi_dot([C, numpy.diag(numpy.sqrt(n)), C.T])
  def _make_D(self, n, C):
      return numpy.linalg.multi_dot([C, numpy.diag(n), C.T])

  # ----- ST-Proj-R Algorithm ----- #

  def _run_proj_R(self, dmft, maxit, conv, g_0, step, verbose, n, C, **kwargs):
      "Run DMFT-SCF with the Steepest-Descents minimization and density matrix projection algorithm (protected interface)"

      nn = len(n)
      iteration = 0
      success = False

      # Exchange-Correlation coefficients
      fij = self.fij(n, dmft, **kwargs)

      # Starting energy
      self.E = self.static_energy(dmft, **kwargs)
      E_old = self.E
      if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_old))

      # First iteration
      iteration += 1
      C_ = numpy.dot(self.Y, C)
      x_old_2 = numpy.hstack([n, C_.ravel()])
      gradient_2 = self._gradient(n, fij, C, dmft, step, **kwargs)
      x_old_1 = x_old_2 - g_0 * gradient_2

      n_old_1 = x_old_1[:nn]
      C_old_1_= x_old_1[nn:].reshape(nn,nn) ; C_old_1 = numpy.dot(self.X, C_old_1_)
      n_old_1, C_old_1 = self._density_matrix_projection(n_old_1, C_old_1)
      self.D = numpy.linalg.multi_dot([C_old_1, numpy.diag(n_old_1), C_old_1.T])
      self.n = n_old_1.copy()
      self.c = C_old_1.copy()
      E_new = self.static_energy(dmft, **kwargs)
      if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_new))

      # Prepare for further iterations
      fij_old_1 = self.fij(n_old_1, dmft, **kwargs)

      n_old_2 = n.copy()
      C_old_2 = C.copy()
      fij_old_2 = fij.copy()

      # Continue iterations with changing step of SD
      stop = False
      iteration += 1
      while stop is False:
          # New density matrix and XC coefficients
          n_new, C_new = self._st_step(n_old_1, fij_old_1, C_old_1, n_old_2, fij_old_2, C_old_2, dmft, step, **kwargs)
          n_new, C_new = self._density_matrix_projection(n_new, C_new)
          fij_new = self.fij(n_new, dmft, **kwargs)
          self.D = numpy.linalg.multi_dot([C_new, numpy.diag(n_new), C_new.T])
          self.n = n_new.copy()
          self.c = C_new.copy()

          # compute current energy
          E_new = self.static_energy(dmft, **kwargs)
          if verbose: print(" @DMFT-SCF Iter %2d. E = %14.8f" % (iteration, E_new))

          # check if converged?
          if abs(E_new-E_old) < conv: 
             stop = True
             success = True

          # check if max iterations were exceeded
          if iteration >= maxit: 
             stop = True
             success = False

          # prepare for next iteration
          iteration += 1
          E_old = E_new
          self.E = E_old

          n_old_2 = n_old_1.copy()
          fij_old_2 = fij_old_1.copy()
          C_old_2 = C_old_1.copy()

          n_old_1 = n_new.copy()
          fij_old_1 = fij_new.copy()
          C_old_1 = C_new.copy()

      if verbose and success:
         print(" DMFT-SCF iterations converged.")
         print(" Final Energy = %14.8f" % self.E)
      if verbose and not success:
         print(" DMFT-SCF iterations did not converge.")
      return success


  def _gradient_proj_R(self, n, fij, C, dmft, step, **kwargs):
      "Gradient vector in the parameter space of n x C.ravel(). Returns gradient in orthogonal AO basis"
      dE_D = self._grad_D(n, fij, C, dmft, step, **kwargs)
      gradient = dE_D
      return gradient

  def _gradient_proj_P(self, n, c, dmft, **kwargs):
      # 1-electron contribution + 2-electron Hartree contribution
      if dmft.lower()!= 'mbb': raise NotImplementedError

      P  = self._make_P(n, c)
      PS = numpy.dot(P, self.S)
      D  = numpy.linalg.multi_dot([P, self.S, P])
      I  = numpy.identity(len(P))

      self._jk.C_clear()
      self._jk.C_left_add(psi4.core.Matrix.from_array(D, ""))
      self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      self._jk.compute()
      J = self._jk.J()[0].to_array(dense=True)

      Jh = self.H + 2.0*J
      grad = numpy.dot(Jh, PS)
      grad+= grad.T
      

      # 2-electron contribution: Exchange-Correlation
      #f  = self.fij(n, dmft, **kwargs)
      #dP = self._grad_P_xc(f, n, c)
      dP = self._jk.K()[0].to_array(dense=True)
      grad -= dP
      # 
      Si = numpy.linalg.inv(self.S)
      grad = 2.0 * numpy.linalg.multi_dot([Si, grad, Si])
      return grad

  def _st_step_proj_P(self, n_old_1, fij_old_1, C_old_1, 
                            n_old_2, fij_old_2, C_old_2,
                            dmft, step, **kwargs):
      raise NotImplementedError
      "Next Steepest-Descents Step"
      nn = len(n_old_1)

      # transform C to orthogonal AO basis
      C_old_1_ = numpy.dot(self.Y, C_old_1)
      C_old_2_ = numpy.dot(self.Y, C_old_2)

      P_old_1_ = self._make_P(n_old_1, C_old_1_)
      P_old_2_ = self._make_P(n_old_2, C_old_2_)

      x_old_1= P_old_1_.ravel()
      x_old_2= P_old_2_.ravel()

      gradient_1 = self._gradient_proj_P(n_old_1, fij_old_1, C_old_1, dmft, step, **kwargs)
      gradient_2 = self._gradient_proj_P(n_old_2, fij_old_2, C_old_2, dmft, step, **kwargs)

      norm = numpy.linalg.norm(gradient_1 - gradient_2)
      g = numpy.dot(x_old_1 - x_old_2, gradient_1 - gradient_2) / norm**2
      #print(" G = ", g)
      g = abs(g)
      x_new = x_old_1 - g * gradient_1
      n_new = x_new[:nn]
      C_new_= x_new[nn:].reshape(nn,nn)

      # transform C back to original AO basis
      C_new = numpy.dot(self.X, C_new_)

      return n_new, C_new



  # ----- Density-Matrix Projection Algorithm ----- #

  def _a(self, a, mu):
      "Projected occupation numbers"
      a_ = a.copy();
      for i in range(len(a)):
          u = a[i] + mu
          if   u <= 0.0: a_[i] = 0.0
          elif u >= 1.0: a_[i] = 1.0
          else: a_[i] = u
      return a_

  def _find_mu(self, n):
      "Search for mu"
      mu = 0.0
      def obj(mu, x):
          u = self._a(x, mu)
          Z = (u.sum() - self.Np)**2
          return Z

      R = scipy.optimize.minimize(obj, mu, args=(n,))
      mu = R.x
      return mu

  def _density_matrix_projection(self, n, C, mu_tol=0.00000000001, mu_iter=1000):
      "Find n_new and C_new such that new density matrix is N-representable"
      # compute pre-density matrix
      preD = numpy.linalg.multi_dot([C, numpy.diag(n), C.T]) # cannot be here self.D because it is pre-density matrix!
      D_ = self._orthogonalize_OPDM(preD, self.S)
      a, b = numpy.linalg.eigh(D_)
      #print(" Init sum = %14.6f" % a.sum())

      mu = self._find_mu(a)

      # compute the projected density matrix
      n_new = self._a(a, mu)
      #print(" Nsum = ", n_new.sum())
      C_new_= b

      # sort (descending order)
      idx = numpy.argsort(n_new)[::-1]
      n_new = n_new [  idx]
      C_new_= C_new_[:,idx]
      C_new = numpy.dot(self.X, C_new_)
      return n_new, C_new

  def _density_matrix_projection_old(self, n, C, mu_tol=0.00000000001, mu_iter=1000):
      "Find n_new and C_new such that new density matrix is N-representable"
      # compute pre-density matrix
      preD = numpy.linalg.multi_dot([C, numpy.diag(n), C.T]) # cannot be here self.D because it is pre-density matrix!
      D_ = self._orthogonalize_OPDM(preD, self.S)
      a, b = numpy.linalg.eigh(D_)
      # 
      Zk = lambda mu, x: (x.sum() - self.Np)**2
      print(" Init sum = %14.6f" % a.sum())

      # start secant search for root mu such that Zk = 0
      stop = False
      mu_0=-0.0
      mu_1= 0.02
      a_0 = self._a(a, mu_0)
      a_1 = self._a(a, mu_1)
      Z_0 = Zk(mu_0, a_0)
      Z_1 = Zk(mu_1, a_1)

      niter = 0
      while stop is False:
          dZ = Z_1 - Z_0
          #if abs(dZ) > mu_tol: mu_k = mu_1 - Z_1 * (mu_1 - mu_0) / dZ
          #else: break
          mu_k = mu_1 - Z_1 * (mu_1 - mu_0) / dZ
          a_k = self._a(a, mu_k)
          Z_k = Zk(mu_k, a_k)
          #if Z_k < mu_tol or niter > mu_iter: stop = True
          if Z_k < mu_tol: stop = True
          mu_0 = mu_1  ; Z_0 = Z_1
          mu_1 = mu_k  ; Z_1 = Z_k
          a_1 = a_k
          niter += 1

      # compute the projected density matrix
      n_new = a_k
      print(" Nsum = ", n_new.sum())
      C_new_= b

      # sort (descending order)
      idx = numpy.argsort(n_new)[::-1]
      n_new = n_new [  idx]
      C_new_= C_new_[:,idx]
      C_new = numpy.dot(self.X, C_new_)
      return n_new, C_new



  # ----- Energy Components ----- #

  def _compute_hcore_energy(self):
      "1-Electron Energy"
      H = self.H.copy()
      D = self.D.copy()
      E = 2.0 * self.compute_1el_energy(D, H)
      return E

  def _compute_hartree_energy(self):
      "2-Electron Energy: Hartree"
      D = self.D.copy()
      E = 2.0 * self.compute_2el_energy(D, D, type='j')
      return E

  def _compute_XC_energy(self, dmft, **kwargs):
      "2-Electron Energy: Exchange-Correlation"
      E = None
      MO_list = ('bb1', 'bb2', 'bb3', 'pow', 'gu', 'cga', 'chf', 'mbb_pc', 'bbc1', 'bbc2', 'mbb0', 'bb_pade_p',
              'bbh', 'bbh2', 'bbt', 'bbg', 'bbb', 'bbx', 'bbq2', 'bby', 'bbp', 'bbz', 'bbi', 'bb_opt', 'bb_opt_2', 'bb_opt_2_p')

      # Summation in AO space using JK object
      if dmft.lower() not in MO_list:
         if   dmft.lower() == 'mbb':  D = self._generalized_density_matrix(numpy.sqrt(self.n), self.c) 
         elif dmft.lower() == 'hf':   D = self._generalized_density_matrix(self.n, self.c)
         E = -1.0 * self.compute_2el_energy(D, D, type='k')
      # summation in MO space
      else:
         fij = self.fij(self.n, dmft, **kwargs)
         Kij = self.Kij(self.c)
         E = -1.0 * numpy.dot(fij, Kij).trace()
      return E


  # ----- DMFT Functionals: Exchange-Correlation Coefficients ----- #

  def _pc(self, n):
      "Phase correction according to BBC1 functional"
      m = n.copy(); m.fill(0.0)
      m[numpy.where(n>=0.5)] =  1.0  # strong
      m[numpy.where(n< 0.5)] = -1.0  # weak
      pc = numpy.ones((len(n), len(n)))
      for i in range(len(n)):
          for j in range(len(n)):
              if i!=j:
                if m[i] < 0 and m[j] < 0:
                  pc[i,j] = -1.0
      return pc
  def _fij_hf(self, n):
      "The Hartree-Fock (HF) Functional"
      return numpy.outer(n, n)
  def _fij_mbb(self, n):
      "The Muller-Buijse-Baerends (MBB) Functional"
      ns = numpy.sqrt(n)
      return numpy.outer(ns, ns)
  def _fij_mbb_0(self, n):
      "The MBB Functional with Zero Exchange (MMB0)"
      f = numpy.diag(n)
      return f
  def _fij_bbc1(self, n):
      "The BBC1 Functional"
      f = self._fij_mbb_PC(n)
      return f
  def _fij_mbb_PC(self, n):  # remove
      f = self._fij_mbb(n) * self._pc(n)
      return f
  def _fij_bbc2(self, n):
      "The BBC2 Functional"
      f = self._fij_mbb_PC(n)
      f_hf = self._fij_hf(n)
      m = n.copy(); m.fill(0.0)
      m[numpy.where(n>=0.5)] =  1.0  # strong
      m[numpy.where(n< 0.5)] = -1.0  # weak
      for i in range(len(n)):
          for j in range(len(n)):
              if i!=j:
                 if m[i] > 0 and m[j] > 0:
                    f[i,j] = f_hf[i,j]
      return f
  def _fij_power(self, n, p):
      "The Power Functional"
      nn = numpy.outer(n, n)
      return nn**p
  def _fij_chf(self, n):
      "The Corrected Hartree-Fock (CHF) Functional"
      nm = n.copy()
      nm[numpy.where(nm<0.0)] = 0.0
      f = self._fij_hf(nm)
      ns = nm * (1.0 - nm)
      ns[numpy.where(ns<0.0)] = 0.0
      f += numpy.sqrt(numpy.outer(ns, ns))
      return f
  def _fij_cga(self, n):
      "The CGA Functional"
      f = self._fij_hf(n)
      ns = n * (2.0 - n)
      f += numpy.sqrt(numpy.outer(ns, ns))
      f/= 2.0
      return f
  def _fij_gu(self, n):
      "The GU Functional"
      nm = n.copy()
      nm[numpy.where(nm<0.0)] = 0.0
      f = self._fij_mbb(nm)
      ns = nm - nm*nm
      f -= numpy.diag(ns)
      return f
  def _fij_bbbk(self, n, k):
      "The BBB-(k) Functional"
      W = 1.0/float(k+1)
      K = float(k*(k+1))
      U = float(k+1)
      f = self._fij_mbb(n)
      f = 2.0**W * ( f**U )
      de= ( n[:,numpy.newaxis]**K + n[numpy.newaxis,:]**K )**W
      for i in range(len(n)):
          for j in range(len(n)):
              de_ij = de[i,j]
              if de_ij < 1.0e-20: f[i,j] = 0.0
              else: f[i,j] = f[i,j] / de_ij
      return f
  def _fij_bbb(self, n, a0, kmax):
      "The MBB-MBB0 Interpolation Functional with Exponential Decay (MBB/ED)"
      f = a0 * self._fij_mbb(n)
      a_sum = 0.0
      if a0<1:
         for k in range(1,kmax+1):                   
             ak = a0 * math.exp(k*math.log(1.0 - a0))
             f += ak * self._fij_bbbk(n, k)
             a_sum += ak
      #else: f = self._fij_mbb(n)
      #print( " Sum of a: %13.4f" % a_sum)
      return f
  def _fij_bbb2(self, n, a0, b0, kmax):
      "The MBB-MBB0 Interpolation Functional with Damped Oscillatory Exponential Decay (MBB/DED)"
      return NotImplementedError


  # other functionals

  def _fij_bbh(self, n):
      f = 2.0 * self._fij_hf(n)
      f/= n[:,numpy.newaxis] + n[numpy.newaxis,:]
      return f
  def _fij_bbh2(self, n):
      f = math.sqrt(2.0) * self._fij_hf(n)
      f/= numpy.sqrt(n[:,numpy.newaxis]**2 + n[numpy.newaxis,:]**2)
      return f
  def _fij_bb1(self, n, eps=1.0e-13):
      N = n.size
      f = numpy.zeros((N,N), numpy.float64)
      for i in range(N):
          ni = n[i]
          for j in range(N):
              nj = n[j]
              if   ni < eps or nj < eps:  f[i,j] = 0.0
              elif abs(ni - nj) < eps: f[i,j] = ni
              else: f[i,j] = (ni - nj) / (math.log(ni) - math.log(nj))
      return f
  def _fij_bb2(self, n, gamma, eps=1.0e-20):
      print("\n")
      print( " Sum of n = %14.8f" % n.sum())
      g = 1.0 - float(gamma)
      g_th = n.sum() / math.sqrt(self._fij_mbb(n).sum())
      s = self._fij_mbb(n).sum()
      s2= self._fij_hf (n).sum()
      s3= numpy.sqrt(n).sum()

      i_N = (n*(1.0 - n)).sum()
      i_D = numpy.sqrt(n*(1.0 - n)).sum() / 2.0 - i_N
      i_T = i_D + i_N

      #s = 1.0 - (g_th + math.sqrt(g_th)) / 2.0
      print ( " GAMMA = %14.6f" % g_th) 
      print ( " S     = %14.6f" % s) 
      print ( " S2    = %14.6f" % s2) 
      print ( " S3    = %14.6f" % s3) 
      print ( " I_DYNAMIC = %14.6f" % i_D)
      print ( " I_NONDYNA = %14.6f" % i_N)
      print ( " I_TOTAL   = %14.6f" % i_T)

      f = g * self._fij_mbb(n) + (1.0 - g) * self._fij_bb1(n, eps)
      return f
  def _fij_bbq2(self, n):
      N = n.sum()
      g = N / math.sqrt(self._fij_mbb(n).sum())
      A = 1./math.pi**2
      B = 0.5
      A = 0.103393; B = 0.452121
      a_0 = 1.0 + (g - 1.0) * (B + N*A)
      a_1 = 1.0 - a_0
      f = a_0 * self._fij_mbb(n) + a_1 * numpy.diag(n)
      return f
  def _fij_bb3(self, n, eps=1.0e-20):
      B = -0.07751677
      C = +0.008
      D = 0.00
      A = 1.0 - B - C - D
      f1 = self._fij_mbb(n)
      f2 = self._fij_bb1(n, eps)
      f3 = self._fij_bbh(n)
      f4 = self._fij_bbh2(n)
      f  = A * f1 + B * f2 + C * f3
      return f

  def _fij_bbg(self, n, coeffs, kmax):
      "Taylor expansion of unknown Fundamental Generator"
      N_terms = kmax + 1
      assert(len(coeffs) + 1 == N_terms), " Number of Coefficients not consistent with k_max"

      A = 1.0 - sum(coeffs)       # MBB term (k=0)
      Bs = coeffs                 # k=1,2,.. terms
                                                        
      # MBB term (k=0)
      f  = A * self._fij_mbb(n)
      # post-MBB terms (k=1,2,...)
      for k in range(1, kmax+1):
          f += Bs[k-1] * self._fij_bbbk(n, k)
      return f

  def _fij_bbt(self, n, bbt_coeffs, bbt_kmax, bbt_hf, bbt_log, eps=1.0e-20):
      "Taylor expansion of unknown Fundamental Generator"
      N_terms = bbt_kmax
      if bbt_hf : N_terms += 1
      if bbt_log: N_terms += 1
      if bbt_mbb: N_terms += 1
      assert(len(bbt_coeffs) + 1 == N_terms), " Number of Coefficients not consistent with k_max"

      # include MBB term
      #if bbt_mbb:
      A = 1.0 - sum(bbt_coeffs)   # MBB term            
      Ds = bbt_coeffs[:bbt_kmax]      # k=1,2,.. terms
      if bbt_log and bbt_hf:
         C = bbt_coeffs[bbt_kmax]       # k=0 term (HF)
         B = bbt_coeffs[bbt_kmax+1]     # log term
      elif not bbt_log and bbt_hf:
         C = bbt_coeffs[bbt_kmax]       # k=0 term (HF)
      elif bbt_log and not bbt_hf:
         B = bbt_coeffs[bbt_kmax]       # log term
      else: pass
                                                        
      # MBB term
      f  = A * self._fij_mbb(n)
      # log term
      if bbt_log: f+= B * self._fij_bb1(n, eps)
      # k=0 term (HF)
      if bbt_hf: f += C * self._fij_hf(n)
      # k=1,2,... terms
      for k in range(1, bbt_kmax+1):
          f += Ds[k-1] * self._fij_bbk(n, k)

      ## exclude MBB term (pure expansion)
      #else:
      #   A = 1.0 - sum(bbt_coeffs)   # k=1 term
      #   Ds = bbt_coeffs[:bbt_kmax-1]      # k=2,.. terms
      #   if bbt_log and bbt_hf:
      #      C = bbt_coeffs[bbt_kmax-1]       # k=0 term (HF)
      #      B = bbt_coeffs[bbt_kmax]         # log term
      #   elif not bbt_log and bbt_hf:
      #      C = bbt_coeffs[bbt_kmax-1]       # k=0 term (HF)
      #   elif bbt_log and not bbt_hf:
      #      B = bbt_coeffs[bbt_kmax-1]       # log term
      #   else: pass
      #                                                     
      #   # MBB term
      #   f  = A * self._fij_bbk(n, 1.0)
      #   # log term
      #   if bbt_log: f+= B * self._fij_bb1(n, eps)
      #   # k=0 term (HF)
      #   if bbt_hf: f += C * self._fij_hf(n)
      #   # k=1,2,... terms
      #   for k in range(2, bbt_kmax+1):
      #       f += Ds[k-2] * self._fij_bbk(n, k)
      return f
  def _fij_bbk(self, n, k):
      f = 2.0**(1.0/k) * self._fij_hf(n)
      f/= ( n[:,numpy.newaxis]**k + n[numpy.newaxis,:]**k )**(1.0/k)
      return f
  def _fij_bbb(self, n, kmax, scale):
      N = n.sum()
      g = N / math.sqrt(self._fij_mbb(n).sum())
      g = g*scale
      A = 1./math.pi**2
      B = 0.5
      #A = 0.103393; B = 0.452121
      a_0 = 1.0 + (g - 1.0) * (B + N*A)
      a_1 = 1.0 - a_0
      s = a_1 / float(kmax)
      f = a_0 * self._fij_mbb(n)
      for k in range(1,kmax+1):
          f += s * self._fij_bbbk(n, k)
      return f
  def _fij_bby(self, n, coeffs, kmax):
      print(" Coefficients: ", coeffs)
      N = n.sum()
      g = N / math.sqrt(self._fij_mbb(n).sum())
      A = 1./math.pi**2
      B = 0.5
      #A = 0.103393; B = 0.452121
      a_0 = 1.0 + (g - 1.0) * (B + N*A)
      f = a_0 * self._fij_mbb(n)
      # coefficients
      a = 0.0 # coeffs[0]
      b = coeffs[0]
      c = 2.0 * (1.0 - a_0)
      c/= -1.0 + math.sinh(b) / (math.cosh(b) - math.cos(a))
      #arg = math.cosh(b) - math.sinh(b) / (2.0 * (1.0 - a_0) - c)
      #print(arg)
      a_sum = a_0
      for k in range(1,kmax+1):
          a_k = c * math.cos(a*k) * math.exp(-b*k)
          f += a_k * self._fij_bbbk(n, k)
          a_sum += a_k
      print( " Sum of a: %13.4f" % a_sum)
      return f
  def _fij_bbz(self, n, coeffs, kmax, pc):
      #print(" Coefficients: ", coeffs[0])
      a = coeffs[0]
      a_sum = a
      f = a * self._fij_mbb(n)
      for k in range(1,kmax+1):
          a_k = a * math.exp(k*math.log(1.0 - a))
          f += a_k * self._fij_bbbk(n, k)
          a_sum += a_k
      print( " Sum of a: %13.4f" % a_sum)
      if pc: f *= self._pc(n)
      return f
  def _fij_bbopt(self, n, kmax, coeffs):
      "well-optimized dmft functional"
      def calc_a0_well_2(I_D, I_N, S, A, B, C, D, E, F, G, H):
          x = math.log(I_D/S + 1.0)
          y =-math.log(2.0 * I_N/S)
          a_0 = A + B*y + C*y*y + D*x*y + E*x*y*y + F*x + G*x*x + H*x*x*y
          a_0 = math.exp(a_0)
          return a_0

      I_n = (n*(1.0 - n)).sum()
      I_d = numpy.sqrt(n*(1.0 - n)).sum() / 2.0 - I_n
      S   = numpy.sqrt(n).sum()
      N   = n.sum()

      A, B, C, D, E, F, G, H = coeffs
      a = calc_a0_well_2(I_d, I_n, S, A, B, C, D, E, F, G, H)
      #print(" A0(BBI) = %13.6f" % a)

      a_sum = a
      f = a * self._fij_mbb(n)
      if a<1:
         for k in range(1,kmax+1):                   
             a_k = a * math.exp(k*math.log(1.0 - a))
             f += a_k * self._fij_bbbk(n, k)
             a_sum += a_k
      #else: f = self._fij_mbb(n)
      #print( " Sum of a: %13.4f" % a_sum)
      return f
  def _fij_bbpade(self, n, kmax, coeffs):
      "well-optimized dmft functional"
      def calc_a0_well_2(I_D, I_N, S, A0, A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6):
          x = math.log(I_D/S + 1.0)
          y =-math.log(2.0 * I_N/S)
          a_0 = (A0+A1*x+A2*y+A3*x*y+A5*y*y+A6*x*y*y)/(1.0+B1*x+B2*y+B3*x*y+B5*y*y+B6*x*y*y)
          return a_0

      I_n = (n*(1.0 - n)).sum()
      I_d = numpy.sqrt(n*(1.0 - n)).sum() / 2.0 - I_n
      S   = numpy.sqrt(n).sum()
      N   = n.sum()

      A0, A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6 = coeffs
      a = calc_a0_well_2(I_d, I_n, S, A0, A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6)
      #print(" A0(BBI) = %13.6f" % a)

      a_sum = a
      f = a * self._fij_mbb(n)
      if a<1:
         for k in range(1,kmax+1):                   
             a_k = a * math.exp(k*math.log(1.0 - a))
             f += a_k * self._fij_bbbk(n, k)
             a_sum += a_k
      #else: f = self._fij_mbb(n)
      #print( " Sum of a: %13.4f" % a_sum)
      return f

  def _fij_bbopt_2(self, n, kmax, coeffs):
      "well-optimized dmft functional"
      def calc_a0_well_2(I_D, I_N, S, A, B, C, D, W):
          x = math.log(I_D/S + 1.0)
          y =-math.log(2.0 * I_N/S)
          a_0 = W/( 1.0 + math.exp(A + B*y + D*y*y + C*x*y)) + math.exp(x) + 1.0 - W
          print(" X = %14.6f  Y = %14.6f  A = %14.6f" % (x, y, a_0))
          return a_0

      I_n = (n*(1.0 - n)).sum()
      I_d = numpy.sqrt(n*(1.0 - n)).sum() / 2.0 - I_n
      S   = numpy.sqrt(n).sum()
      N   = n.sum()

      A, B, C, D, W = coeffs
      a = calc_a0_well_2(I_d, I_n, S, A, B, C, D, W)
      #print(" A0(BBI) = %13.6f" % a)

      a_sum = a
      f = a * self._fij_mbb(n)
      if a<1:
         for k in range(1,kmax+1):                   
             a_k = a * math.exp(k*math.log(1.0 - a))
             f += a_k * self._fij_bbbk(n, k)
             a_sum += a_k
      #else: f = self._fij_mbb(n)
      #print( " Sum of a: %13.4f" % a_sum)
      return f


  def _fij_bbZ(self, n, kmax, pc, ao=None):
      #print(" Coefficients: ", coeffs[0])
      if ao is None:
         c_10 = -0.217033                                                                     
         c_01 = -1.10574
         c_11 = -0.743481
         c_20 =  0.335929
         c_02 =  2.00444
         I_n = (n*(1.0 - n)).sum()
         I_d = numpy.sqrt(n*(1.0 - n)).sum() / 2.0 - I_n
         S   = numpy.sqrt(n).sum()
         if pc:
            c_10 = -0.207746
            c_01 = -0.254876
            c_11 = -0.839247
            c_20 =  0.412964
            c_02 =  0.217702
         a = 1.0 + c_10 * I_d + c_01 * I_n + c_11 * I_d * I_n + c_20 * I_d**2 + c_02 * I_n**2
                                                                                              
         A = 0.245404
         B = 7.9045
         C = -0.407284
         D = 0.114264
         E = 0.258704
         F = 962.558
         G = 0.828886
         H = -0.264658
         #a = G + A*math.exp(-B*I_d) + C*I_n + D*I_d + E*I_n*math.exp(-F*I_d) + H*I_d*I_n

         # new model CCSD/aug-cc-pVTZ for 2-el molecules
         def f1(x, a, b):
             H = 1./(1 + math.exp(-b/a))
             G = 1./(1 + math.exp((1.0-b)/a))
             arg = abs(1./(x*(H-G) + G) - 1.0)
             if arg < 1.e-40: arg = 0.00000000000001
             a0 = b + a * math.log(arg)
             #print(abs(1./(x*(H-G) + G) - 1.0))
             #print(" AO = %14.5f" % a0)
             #exit()
             return a0
         def f2(x, c):
             return math.exp(-c*x)
         def calc_a0(x, y, a, b, c):
             a0 = f2(x, c) * f1(y, a, b)
             return a0
         def calc_a0_lin(I_N, N, S):
             x = math.log(S/N - 1.0)
             y = math.log(2.0 * I_N/N)
             # FCI/STO-3G
             A = 0.739242
             B =-0.349756
             C = 1.56737
             a_0 = A*x + B*y + C
             return a_0
         def calc_a0_well(I_D, I_N, N, A, B, C, D):
             x = math.log(I_D/N + 1.0)
             y =-math.log(2.0 * I_N/N)
             a_0 = math.exp(x)*(A*y*y+B*y+C*x*y+D)
             return a_0
         def calc_a0_well_2(I_D, I_N, S, A, B, C, D, E, F, G, H):
             x = math.log(I_D/S + 1.0)
             y =-math.log(2.0 * I_N/S)
             a_0 = A + B*y + C*y*y + D*x*y + E*x*y*y + F*x + G*x*x + H*x*x*y
             a_0 = math.exp(a_0)
             return a_0


         N = n.sum()

         u = 0.74
         D =-0.04508491
         x = I_d / S**(1.0)
         y = I_n * 2.0 / N
         # with exp(-Cx^2)
         A = 0.04437981 # 0.04450706 # 0.0573076
         B = 0.92331502 # 0.92454078 # 1.05693425
         C = 3.94820784 # 4.00406427 # 4.54670322 

         # with exp(-Cx)
         A = 0.058713  # 0.05871278
         B = 1.5727958 # 1.24226418
         C = 0.88809613 # 0.87847818

         # new fitting (FCI/aug-cc-pvTZ)
         A = 0.03614488 #  0.0389451  #  
         B = 0.90184505 #  0.91592284 #  
         C = 0.70143966 #  0.63247296 #  
         D =-0.04508491
         a = calc_a0(x, y, A, B, C)

         a = calc_a0_lin(I_n, N, S)

         # new fitting (FCI/aug-cc-pVTZ)
         A = 0.0830253
         B =-0.543444
         C = 3.12916
         D = 0.961993
         a = calc_a0_well(I_d, I_n, N, A, B, C, D)

         # new fitting (FCI/aug-cc-pVTZ)
         A =  0.00806763
         B = -0.177793
         C = -0.0139493
         D =-16.8878
         E =  4.10119
         F =  3.00987
         G = 90.0324
         H = 47.6213
         a = calc_a0_well_2(I_d, I_n, S, A, B, C, D, E, F, G, H)
         print(" A0(BBI) = %13.6f" % a)


      else: 
         print(" AO set externally")
         a = ao
      a_sum = a
      f = a * self._fij_mbb(n)
      for k in range(1,kmax+1):
          a_k = a * math.exp(k*math.log(1.0 - a))
          f += a_k * self._fij_bbbk(n, k)
          a_sum += a_k
      print( " Sum of a: %13.4f" % a_sum)
      if pc: f *= self._pc(n)
      return f

  def _fij_bbp(self, n):
      N = n.sum()
      g = N / math.sqrt(self._fij_mbb(n).sum())
      A = 1./math.pi**2
      B = 0.5
      #A = 0.103393; B = 0.452121
      a_0 = 1.0 + (g - 1.0) * (B + N*A)
      f = a_0 * self._fij_mbb(n)
      # coefficients
      a = 0.0
      b = 0.613879 * g * N -5.46879 * g -0.443713 * N + 4.51873
      c = 2.0 * (1.0 - a_0)
      c/= -1.0 + math.sinh(b) / (math.cosh(b) - math.cos(a))
      a_sum = a_0
      for k in range(1,30+1):
          a_k = c * math.cos(a*k) * math.exp(-b*k)
          f += a_k * self._fij_bbbk(n, k)
          a_sum += a_k
      print( " Sum of a: %13.4f" % a_sum)
      return f
  def _fij_BBX(self, n, coeffs, kmax):
      #N_parameters = 2.0
      #assert (len(coeffs) == int(N_parameters)) , "Number of parameters must be 3"
      gamma = n.sum() / math.sqrt(self._fij_mbb(n).sum()) 
      #F = lambda x, A, B, C: A * math.atan(B * x + C)
      F = lambda x, A, B: A * x + B
      f = gamma * self._fij_mbb(n)
      #delta = 1.0 - gamma
      B = coeffs[0]
      kf = float(kmax)
      A = 2.0 * (1.0 - gamma - B * kf) / ((kf + 1.0) * kf)
      AS = numpy.array([F(k, A, B) for k in range(1,kmax+1)])
      #if (delta > 0.0 and AS.sum() > 0.0):
      #   AS/= AS.sum()/delta
      #else: AS = numpy.zeros(kmax, numpy.float64)

      #print(gamma, AS)

      for k in range(1,kmax+1):
          f += AS[k-1] * self._fij_bbbk(n, k)
      return f



  def _fij_1_hf(self, n, m):
      "First analytical derivatives of fij_HF wrt n_m"
      nn = len(n)
      fij_m = numpy.zeros((nn,nn), numpy.float64)
      for i in range(nn):
          for j in range(nn):
              v = 0.0
              if   (i==m) and (j!=m): v = n[j]
              elif (j==m) and (i!=m): v = n[i]
              elif (i==j==m)        : v = n[m]*2.0
              fij_m[i, j] = v
      return fij_m

  def _fij_1_numerical(self, n, m, step, func, **kwargs):
      "First numerical derivatives of fij_func wrt n_m. Uses Forward second-order finite difference with O(step^2)"
      n_p1 = self._n_m(n, m, +1.0*step)
      n_p2 = self._n_m(n, m, +2.0*step)
      fij_m = (-3.0*func(n, **kwargs) + 4.0*func(n_p1, **kwargs) - func(n_p2, **kwargs))/(2.0 * step)
      return fij_m

  def _n_m(self, n, m, step):
      n_new = n.copy()
      n_new[m] += step
      return n_new


  # --- auxiliary --- # 

  def _generalized_density_matrix(self, n, c):
      "Compute occupation-weighted 1-electron density matrix in AO basis"
      ns = n.real.copy(); N = c.shape[0]
      cs = c.real.copy()
      D  = numpy.zeros((N, N), numpy.float64)
      for i in range(n.size):
          D += ns[i] * numpy.outer(cs[:,i], cs[:,i]) 
      return D

  #def _orthogonalizer(self, S):
  #    "Form orthogonalizer"
  #    L, U = numpy.linalg.eigh(S)
  #    L    = numpy.diag(1./numpy.sqrt(L))
  #    X    = numpy.dot(U, numpy.dot(L, U.T))
  #    return X


  # --- COPIED FROM DensityDecomposition -> remove or readapt to reduce code duplication

  def natural_orbitals(self, D, orthogonalize_first=None, order='descending', original_ao_mo=True, renormalize=False, no_cutoff=False):
      "Compute the Natural Orbitals from a given ODPM"
      if orthogonalize_first is not None:
         S = orthogonalize_first
         D_ = self._orthogonalize_OPDM(D, S)
      else:
         D_ = D
      n, U = numpy.linalg.eigh(D_)
      if n.max() > 1.0 or n.min() < 0.0:
         print(" Warning! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))
      #if ((n.max() - 1.0) > 0.00001 or (n.min() < -0.00001)):
      #   raise ValueError("Unphysical NO populations detected! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))
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
         if ( abs(n.sum() - numpy.round(n.sum())) > 1.e-7):
            print(" Warning: nsum=%14.4E delta=%14.4E" % (n.sum(), n.sum() - numpy.round(n.sum())))
         d = numpy.round(n.sum()) - n.sum()
         d/= numpy.float64(n.size)
         n+= d
         n[numpy.where(n<0.0)] = 0.0
         n[numpy.where(n>1.0)] = 1.0
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


  def compute_1el_energy(self, D, Hcore):
      "Compute generalized 1-electron energy"
      energy = numpy.dot(D, Hcore).trace()
      return energy

  def compute_2el_energy(self, D_left, D_right, type='j'):
      "Compute generalized 2-electron energy"
      assert self._jk is not None
      self._jk.C_clear()                                           
      self._jk.C_left_add(psi4.core.Matrix.from_array(D_left, ""))
      I = numpy.identity(D_left.shape[0], numpy.float64)
      self._jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
      self._jk.compute()
      if   type.lower() == 'j': JorK = self._jk.J()[0].to_array(dense=True)
      elif type.lower() == 'k': JorK = self._jk.K()[0].to_array(dense=True)
      else: raise ValueError("Incorrect type of JK matrix. Only J or K allowed.")
      energy = numpy.dot(JorK, D_right).trace()
      return energy

  def matrix_power(self, M, x, eps=1.0e-6):
      "Computes the well-behaved matrix power. All eigenvalues below or equal eps are ignored"
      E, U = numpy.linalg.eigh(M)
      Ex = E.copy(); Ex.fill(0.0)
      for i in range(len(E)):
          if (E[i]>0.0+eps) : Ex[i] = numpy.power(E[i], x)
          else: Ex[i] = 0.0
      Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
      return Mx




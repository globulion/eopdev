#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DFT module.
 Bartosz Błasiak, Gundelfingen, Feb 2019
"""

import os
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
                                              The DFT method
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
                                                                       Last Revision: Gundelfingen, Feb 7th 2019
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
      # Natural Orbitals from current OPDM
      D = self._wfn.Da().to_array(dense=True)
      S = self._wfn.S ().to_array(dense=True)
      self.n, self.c = self.natural_orbitals(D, orthogonalize_first=S, order='descending', no_cutoff=0.0, renormalize=False)
      #print(self.n, self.n.sum())
      return

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
      E_XC = self._compute_XC_energy(self.n, self.c, dmft, **kwargs)
      # Total energy
      E = E_N + E_1 + E_H + E_XC
      return E



  def _compute_hcore_energy(self):
      H = self.H.to_array(dense=True)
      D = self.D.to_array(dense=True)
      E = 2.0 * self.compute_1el_energy(D, H)
      return E
  def _compute_hartree_energy(self):
      D = self.D.to_array(dense=True)
      E = 2.0 * self.compute_2el_energy(D, D, type='j')
      return E
  def _compute_XC_energy(self, n, c, dmft, **kwargs):
      E = None
      MO_list = ('bb1', 'bb2', 'bb3', 'pow', 'gu', 'cga', 'chf', 'mbb_pc', 'bbc1', 'bbc2', 'mbb0',
              'bbh', 'bbh2', 'bbt', 'bbg', 'bbb', 'bbx', 'bbq2', 'bby', 'bbp', 'bbz', 'bbi')

      # Summation in AO space using JK object
      if dmft.lower() not in MO_list:
         if   dmft.lower() == 'mbb':  D = self._generalized_density_matrix(numpy.sqrt(n), c) 
         elif dmft.lower() == 'hf':   D = self._generalized_density_matrix(n, c)
         E = -1.0 * self.compute_2el_energy(D, D, type='k')
      # summation in MO space
      else:
         fij = self.fij(n, dmft, **kwargs)
         Kij = self.Kij(c)
         E = -1.0 * numpy.dot(fij, Kij).trace()
      return E

  def fij(self, n, dmft, **kwargs):
      if   dmft.lower() == 'mbb': f = self._fij_mbb(n)
      elif dmft.lower() == 'mbb_pc': f = self._fij_mbb_PC(n)  # same as bbc1
      elif dmft.lower() == 'mbb0': f = self._fij_mbb_0(n)
      elif dmft.lower() == 'hf' : f = self._fij_hf (n)
      elif dmft.lower() == 'chf': f = self._fij_chf(n)
      elif dmft.lower() == 'cga': f = self._fij_cga(n)
      elif dmft.lower() == 'gu' : f = self._fij_gu (n)
      elif dmft.lower() == 'bbc1':f = self._fij_bbc1(n)
      elif dmft.lower() == 'bbc2':f = self._fij_bbc2(n)
      elif dmft.lower() == 'bb1': f = self._fij_bb1(n)
      elif dmft.lower() == 'bb2': f = self._fij_bb2(n, kwargs["gamma"])
      elif dmft.lower() == 'bb3': f = self._fij_bb3(n)
      elif dmft.lower() == 'bbp': f = self._fij_bbp(n)
      elif dmft.lower() == 'bbb': f = self._fij_bbb(n, kwargs["kmax"], kwargs["scale"])
      elif dmft.lower() == 'bby': f = self._fij_bby(n, kwargs["coeffs"], kwargs["kmax"])
      elif dmft.lower() == 'bbz': f = self._fij_bbz(n, kwargs["coeffs"], kwargs["kmax"], kwargs["pc"])
      elif dmft.lower() == 'bbi': f = self._fij_bbZ(n, kwargs["kmax"], kwargs["pc"], kwargs["ao"])
      elif dmft.lower() == 'bbq2': f = self._fij_bbq2(n)
      elif dmft.lower() == 'bbh': f = self._fij_bbh(n)
      elif dmft.lower() == 'bbh2':f = self._fij_bbh2(n)
      elif dmft.lower() == 'pow': f = self._fij_power(n, kwargs["pow_exp"])
      elif dmft.lower() == 'bbt': f = self._fij_bbt(n, kwargs["bbt_coeffs"], kwargs["bbt_kmax"], kwargs["bbt_hf"], kwargs["bbt_log"])
      elif dmft.lower() == 'bbg': f = self._fij_bbg(n, kwargs["coeffs"], kwargs["kmax"])
      elif dmft.lower() == 'bbx': f = self._fij_BBX(n, kwargs["coeffs"], kwargs["kmax"])
      else: raise NotImplementedError
      return f

  def Kij(self, c):
      "Exchange (ij|ji) integral matrix"
      C = psi4.core.Matrix.from_array(c, "Natural Orbitals LCAO matrix")
      #eri_K_ij = oepdev.calculate_Kij(self._wfn.c1_deep_copy(self._wfn.basisset()), C).to_array(dense=True)
      eri_K_ij = oepdev.calculate_Kij(self._wfn, C).to_array(dense=True)
      psi4.core.clean()

      #bfs = self._wfn.basisset()
      #mints = psi4.core.MintsHelper(bfs)
      #eri = numpy.asarray(mints.ao_eri(bfs, bfs,
      #                                 bfs, bfs))
      #eri_K_ij = numpy.einsum("ijkl,ia,jb,kb,la->ab", eri, c, c, c, c) 
      return eri_K_ij

  # --- exchange-correlation density functionals --- #

  def _pc(self, n):
      "phase correction according to BBC1 functional"
      m = n.copy(); m.fill(0.0)
      m[numpy.where(n>=0.5)] =  1.0  # strong
      m[numpy.where(n< 0.5)] = -1.0  # weak
      pc = numpy.ones((len(n), len(n)))
      for i in range(len(n)):
          for j in range(len(n)):
              if i!=j:
                if m[i] < 0 and m[j] < 0:
                  pc[i,j] = -1.0
      #print(pc)
      return pc

  def _fij_hf(self, n):
      return numpy.outer(n, n)
  def _fij_mbb(self, n):
      ns = numpy.sqrt(n)
      return numpy.outer(ns, ns)
  def _fij_mbb_PC(self, n):
      "Equivalent to BBC1"
      f = self._fij_mbb(n) * self._pc(n)
      return f
  def _fij_mbb_0(self, n):
      "MBB without exchange at all"
      f = numpy.diag(n)
      return f
  def _fij_bbc1(self, n):
      f = self._fij_mbb_PC(n)
      return f
  def _fij_bbc2(self, n):
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
      nn = numpy.outer(n, n)
      return nn**p
  def _fij_chf(self, n):
      f = self._fij_hf(n)
      ns = n * (1.0 - n)
      f += numpy.sqrt(numpy.outer(ns, ns))
      return f
  def _fij_cga(self, n):
      f = self._fij_hf(n)
      ns = n * (2.0 - n)
      f += numpy.sqrt(numpy.outer(ns, ns))
      f/= 2.0
      return f
  def _fij_gu(self, n):
      f = self._fij_mbb(n)
      ns = n - n*n
      f -= numpy.diag(ns)
      return f
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


  def _fij_bbbk(self, n, k):
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



  def _generalized_density_matrix(self, n, c):
      "Compute occupation-weighted 1-electron density matrix in AO basis"
      ns = n.real.copy(); N = c.shape[0]
      cs = c.real.copy()
      D  = numpy.zeros((N, N), numpy.float64)
      for i in range(n.size):
          D += ns[i] * numpy.outer(cs[:,i], cs[:,i]) 
      return D

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
      raise NotImplementedError
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
      if ((n.max() - 1.0) > 0.00001 or (n.min() < -0.00001)):
         raise ValueError("Unphysical NO populations detected! nmax=%14.4E nmin=%14.4E" % (n.max(), n.min()))
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




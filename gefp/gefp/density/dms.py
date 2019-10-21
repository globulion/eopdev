#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 One-Particle Density Matrix Susceptibility Module.
 Bartosz Błasiak, Gundelfingen, May 2019

 Implements the charge-transfer DMS tensors.
"""

import numpy
import psi4, oepdev
import scipy.optimize
import scipy.spatial.transform
from abc import ABC, abstractmethod
from .partitioning import DensityDecomposition
from ..math.matrix import move_atom_rotate_molecule, rotate_ao_matrix, matrix_power

__all__ = ["DMSFit", "DMS", "Computer"]

class Computer: #TODO
  def __init__(self, *dms):
      self._set_work(*dms)

  def _set_work(self, *dms):
      self._data = {}
      self._data["dms"] = dms
      #self._data["nbf"] = numpy.array([x._bfs.nbf() for x in dms])
      #self._data["d"] = numpy.array([x._D for x in dms])

      self._dim = sum(self._data["nbf"])
      self._work = numpy.zeros((self._dim, self._dim))

  def deformation_density(self, return_all=False):
      dD_pau = self._deformation_density_pauli()
      dD_pol = self._deformation_density_polar()
      dD = dD_pau + dD_pol
      if return_all: return dD_pau, dD_pol, dD
      return dD

  # protected
  @abstractmethod
  def _deformation_density_pauli(self): pass

  @abstractmethod
  def _deformation_density_polar(self): pass
#      for i in range(self._nfrag):
#          for j in range(self._nfrag):
#              computer_ij = Computer_Dimer(*self._dms_set[i,j])

class UnperturbedField_Computer(Computer): pass

#class DMS_Computer(Computer):
#
#  def __init__(self, *dms):
#      Computer.__init__(self, *dms)
#
#  def _deformation_density_pauli(self):
#      raise NotImplementedError
#      S_mo_t    = numpy.zeros((self._nmo_t, self._nmo_t), numpy.float64)
#      W_mo_t    = numpy.zeros((self._nmo_t             ), numpy.float64)
#      C_ao_mo_t = numpy.zeros((self._nbf_t, self._nmo_t), numpy.float64)
#
#      for i in range(self._nfrag):
#          bfs_i = self._nbf_set[i]
#          non_i = self._data["non"][i]
#          noc_i = self._data["noc"][i]
#          nmo_i = self._data["nmo"][i]
#          nbf_i = self._data["nbf"][i]
#          ofm_i = self._data["ofm"][i]
#          ofb_i = self._data["ofb"][i]
#          S_ao_ii = self._mints.ao_overlap(
#          S_mo_ii = noc_i.T @ S_ao_ii @ noc_i
#          S_mo_t[ofm_i:ofm_i+nmo_i, ofm_i:ofm_i+nmo_i] = S_mo_ii
#          W_mo_t[ofm_i:ofm_i+nmo_i] = numpy.sqrt(non_i)
#          C_ao_mo_t[ofb_i:ofb_i+nbf_i,ofm_i:ofm_i+nmo_i] = noc_i
#
#          for j in range(i):
#              wfn_j = self.data["wfn"][j]
#              noc_j = self.data["noc"][j]
#              nmo_j = self.data["nmo"][j]
#              nbf_j = self.data["nbf"][j]
#              ofm_j = self.data["ofm"][j]
#              ofb_j = self.data["ofb"][j]
#              S_ao_ij = self._mints.ao_overlap(bfs_i, bfs_j)
#              S_mo_ij = noc_i.T @ S_ao_ij @ noc_j
#              S_mo_t[ofm_i:ofm_i+nmo_i, ofm_j:ofm_j+nmo_j] = S_mo_ij 
#              S_mo_t[ofm_j:ofm_j+nmo_j, ofm_i:ofm_i+nmo_i] = S_mo_ij.T
#
#      n_mo_t = W_mo_t * W_mo_t
#
#      WSW = numpy.diag(W_mo_t) @ S_mo_t @ numpy.diag(W_mo_t)                  # mo::mo
#      WSWm12 = matrix_power(WSW, -0.5)                                        # mo::mo
#
#      K = WSWm12 @ numpy.diag(n_mo_t) @ WSWm12                                # mo::mo
#      CW = C_ao_mo_t @ numpy.diag(W_mo_t))                                    # ao::mo
#      Doo_ao_t = CW @ K @ CW.T                                                # ao::ao
#      D_ao_t = C_ao_mo_t @ numpy.diag(n_mo_t) @ C_ao_mo_t.T                   # ao::ao
#
#      # NO analysis of Antisymmetrized wavefunction
#      if False:
#         noo, coo = Density.natural_orbitals(Doo_ao_t.copy(), None, None, 
#                                          orthogonalize_mo = True,
#                                          order='descending', no_cutoff=0.0,
#                                          return_ao_orthogonal = False,
#                                          ignore_large_n = False,
#                                          renormalize = False,
#                                          n_eps = 0.0)
#         print("Sum of natural orbital occupations for noo= %13.6f" % noo.sum())
#
#      # save
#      return Doo_ao_t - D_ao_t




class DMS(ABC):
  def __init__(self, bfs, D):
      self._bfs = bfs            # BasisSet object
      self._mol = bfs.molecule() # Molecule object
      self._D   = D.copy()       # OPDM associated with molecule
      self._B_ind = None         # DMS (induction)
      self._B_ct  = None         # DMS (polarization)

  def rotate(self, rot): 
      # rotate OPDM
      self._D, R = rotate_ao_matrix(self._D, rot, self._bfs, return_rot=True, aomo=False)
      # rotate molecule
      angles = scipy.spatial.transform.Rotation.from_dcm(rot).as_euler('zxy', degrees=True)
      move_atom_rotate_molecule(self._mol, angles, t='zxy')
      # rotate basis set
      raise NotImplementedError
      # rotate DMS tensor
      self._rotate_dms_ind(rot, R)
      self._rotate_dms_ct(R)

  def translate(self, t): 
      # translate molecule
      xyz = self._mol.geometry().to_array(dense=True)
      self._mol.set_geometry(psi4.core.Matrix.from_array(xyz))
      # translate basis set
      for i in range(self._mol.natom()):
          self._bfs.move_atom(i, t)

  def superimpose(self, xyz, suplist=None): 
      raise NotImplementedError

  def _rotate_dms_ind(self, r, R):
      self._B_ind1 = numpy.einsum("ab,cd,uw,bdw->acu", R, R, r, self._B_ind1)
      self._B_ind2 = numpy.einsum("ab,cd,uw,xz,bdwz->acux", R, R, r, r, self._B_ind2)

  def _rotate_dms_ct(self, R):
      self._B_ct = R.T @ self._B_ct @ R

class Ind_DMS(DMS):
  def __init__(self):
      DMS.__init__(self)
  
class DMSFit(ABC):
  def __init__(self, mol_A, mol_B, method, nsampl_pol, nsampl_ct):
      ABC.__init__(self)
      self._mol_A = mol_A
      self._mol_B = mol_B
      self._method = method
      self._nsampl_pol= nsampl_pol
      self._nsampl_ct = nsampl_ct
      psi4.set_options({"DMATPOL_TRAINING_MODE":"CHARGES",
                        "DMATPOL_FIELD_RANK"   : 2,
                        "DMATPOL_GRADIENT_RANK": 0,
                        "DMATPOL_NSAMPLES"     : nsampl_pol,
                        "DMATPOL_NTEST_CHARGE" : 10,
                        "DMATPOL_TEST_CHARGE"  : 0.001})
      self._wfn_A, self._wfn_B = None, None
      self._dms_pol1_A, self._dms_pol1_B = None, None
      self._dms_pol2_A, self._dms_pol2_B = None, None
      self._dms_ct_A, self._dms_ct_B = None, None
      self._nbf_A, self._nbf_B = None, None

  @classmethod
  def create(cls, mol_A, mol_B = None, type="transl", nsampl_pol=100, nsampl_ct=100, method='scf'):
      if type.lower().startswith("tran"): return Asymmetric_Translation_SameMoleculeDMSFit(mol_A, method, nsampl_pol, nsampl_ct)
      else: raise NotImplementedError("This type of DMS fitting is not implemented yet")

  def run(self):
      self._compute_wfn()
     #self._compute_dms_induction()
      self._prepare_for_ct()
      self._compute_dms_ct()

  # --- protected --- #

  def _compute_wfn(self):
      print(" * Computing Unperturbed Wavefunction")
      _g_A, self._wfn_A = psi4.gradient(self._method, molecule=self._mol_A, return_wfn=True)
      self._nbf_A = self._wfn_A.basisset().nbf()
      if self._mol_B is not None:
         _g_B, self._wfn_B = psi4.gradient(self._method, molecule=self._mol_B, return_wfn=True)
         self._nbf_B = self._wfn_B.basisset().nbf()

  def _compute_dms_induction(self): # this is to be reimplemented in child classes
      solver_A = oepdev.GenEffParFactory.build("POLARIZATION", self._wfn_A, self._wfn_A.options())
      self._dms_pol_A = solver_A.compute()
      if self._mol_A is not None:
         solver_B = oepdev.GenEffParFactory.build("POLARIZATION", self._wfn_B, self._wfn_B.options())
         self._dms_pol_B = solver_B.compute()

  @abstractmethod
  def _prepare_for_ct(self): pass

  @abstractmethod
  def _compute_dms_ct(self): pass

  @abstractmethod
  def _construct_aggregate(self): pass

  def _invert_hessian(self, h):
     det= numpy.linalg.det(h)
     print(" * Hessian Determinant= %14.6f" % det)
     hi = numpy.linalg.inv(h)
     I = numpy.dot(hi, h).diagonal().sum()
     d = I - len(hi)
     if abs(d) > 0.0001: raise ValueError("Hessian is problemmatic! I = %f" % I)
     return hi      

class SameMoleculeDMSFit(DMSFit):
  def __init__(self, mol_A, method, nsampl_pol, nsampl_ct):
      DMSFit.__init__(self, mol_A, None, method, nsampl_pol, nsampl_ct)
      self._g_ind = None
      self._h_ind = None
      self._g_ct  = None
      self._h_ct  = None
      self._dD_indA_set_ref = []
      self._dD_indB_set_ref = []
      self._dD_ctAB_set_ref = []
      self._S_AB_set = []
      self._T_set = []
      self._DIP_set = []


 #@abstractmethod
 #def _compute_gradient(self): pass
 #@abstractmethod
 #def _compute_hessian(self): pass

class Translation_SameMoleculeDMSFit(SameMoleculeDMSFit):
  def __init__(self, mol_A, method, nsampl_pol, nsampl_ct):
      SameMoleculeDMSFit.__init__(self, mol_A, method, nsampl_pol, nsampl_ct)
      self._i = -1
      self._start = 1.5
      self._delta = 0.01
      self._t = numpy.array([3.0, 4.0, 5.0]); self._t/= numpy.linalg.norm(self._t)
      self._natoms = self._mol_A.natom()

  def _compute_dms_induction(self):
      self._dms_ind1 = numpy.zeros((self._nbf_A, self._nbf_A, self._natoms, 3))
      self._dms_ind2 = numpy.zeros((self._nbf_A, self._nbf_A, self._natoms, 3, 3))

      for i in range(self._nbf_A):
          for j in range(self._nbf_A):
              b_1, b_2 = self._compute_dms_induction_ij(i, j)
              self._dms_ind1[i,j] = b1
              self._dms_ind2[i,j] = b2

  def _compute_dms_induction_ij(self, i, j):
      #TODO build up gradient and hessian
      self._g_ind.fill(0.0)
      self._H_ind.fill(0.0)
      # fit
      Hi = self._invert_hessian(self._h_ind)
      par = - self._g_ind @ Hi
      # extract
      b_1 = numpy.zeros((self._natoms, 3))
      b_2 = numpy.zeros((self._natoms, 3, 3))
      #TODO fill up b_1 and b_2
      return b_1, b_2

  def _compute_efield(self):
      D = self._wfn_A.Da().to_array(dense=True)
      mints_A = psi4.core.MintsHelper(self._bfs_A)
      mints_B = psi4.core.MintsHelper(self._bfs_B)
      mol_A = self._bfs_A.molecule()      
      mol_B = self._bfs_B.molecule()

      # field on A atoms
      F_A_set = []
      for i in range(mol_A.natom()):
          xyz_i = [mol_A.x(i), mol_A.y(i), mol_A.z(i)]
          ints = mints_B.electric_field(origin=xyz_i)

          #F_nuc_x = ...
          #F_nuc_y = ...
          #F_nuc_z = ...

          Fi_x = (D @ ints[0].to_array(dense=True)).trace() + F_nuc_x
          Fi_y = (D @ ints[1].to_array(dense=True)).trace() + F_nux_y
          Fi_z = (D @ ints[2].to_array(dense=True)).trace() + F_nuc_z

          F_A_set.append([Fi_x, Fi_y, Fi_z])

      F_A_set = numpy.array(F_A_set)

      # field on B atoms
      F_B_set = []
      for i in range(mol_B.natom()):
          xyz_i = [mol_B.x(i), mol_B.y(i), mol_B.z(i)]
          ints = mints_A.electric_field(origin=xyz_i)

          #F_nuc_x = ...
          #F_nuc_y = ...
          #F_nuc_z = ...

          Fi_x = (D @ ints[0].to_array(dense=True)).trace() + F_nuc_x
          Fi_y = (D @ ints[1].to_array(dense=True)).trace() + F_nux_y
          Fi_z = (D @ ints[2].to_array(dense=True)).trace() + F_nuc_z

          F_B_set.append([Fi_x, Fi_y, Fi_z])

      F_B_set = numpy.array(F_B_set)
      return F_A_set, F_B_set
   
  def _compute_T(self):
      D = self._wfn_A.Da().to_array(dense=True)
      return D @ self._S

  def _rms(self, a, b): return numpy.sqrt(((a-b)**2).sum()/a.size)


  def _construct_aggregate(self): 
      self._i += 1
      log = "\n0 1\n"
      for i in range(self._natoms):
          log += "%s" % self._mol_A.symbol(i)
          log += "%16.6f" % (self._mol_A.x(i) * psi4.constants.bohr2angstroms)
          log += "%16.6f" % (self._mol_A.y(i) * psi4.constants.bohr2angstroms)
          log += "%16.6f" % (self._mol_A.z(i) * psi4.constants.bohr2angstroms)
          log += "\n"
      log += "units angstrom\n"
      log += "symmetry c1\n"
      log += "no_reorient\n"
      log += "no_com\n"

      log += "--\n"
      log += "0 1\n"

      t =  (self._start + self._i * self._delta) * self._t
      for i in range(self._natoms):
          log += "%s" % self._mol_A.symbol(i)
          log += "%16.6f" % (self._mol_A.x(i) * psi4.constants.bohr2angstroms + t[0])
          log += "%16.6f" % (self._mol_A.y(i) * psi4.constants.bohr2angstroms + t[1])
          log += "%16.6f" % (self._mol_A.z(i) * psi4.constants.bohr2angstroms + t[2])
          log += "\n"
      log += "units angstrom\n"
      log += "symmetry c1\n"
      log += "no_reorient\n"
      log += "no_com\n"
      print(log)

      mol = psi4.geometry(log)
      mol.update_geometry()
      return mol
 #def _compute_gradient(self): raise NotImplementedError
 #    #TODO
 #def _compute_hessian(self): raise NotImplementedError
 #    #TODO



class Symmetric_Translation_SameMoleculeDMSFit(Translation_SameMoleculeDMSFit):
  def __init__(self, mol_A, method, nsampl_pol, nsampl_ct):
      Translation_SameMoleculeDMSFit.__init__(self, mol_A, method, nsampl_pol, nsampl_ct)

 
  def _prepare_for_ct(self):
      "Compute electric fields, CT kernels and reference deformation density matrices"
      print(" * Computing DDS for Each Sample")
      mints = psi4.core.MintsHelper(self._wfn_A.basisset())
      for n in range(self._nsampl_ct):
          print(" * - Sample %d" % self._i)
          aggr= self._construct_aggregate()
          dds = DensityDecomposition(aggr, method=self._method, acbs=False, jk_type='direct', 
                                           no_cutoff=0.000, xc_scale=1.0, l_dds=False, n_eps=5.0E-5, cc_relax=True,
                                           verbose=False) 
          dds.compute(polar_approx=False)
          self._S = dds.matrix["sqm"][:self._nbf_A,self._nbf_A:] # S_AB
          self._S_AB_set.append(self._S.copy())
          self._bfs_A = dds.data["wfn"][0].basisset()
          self._bfs_B = dds.data["wfn"][1].basisset()
                                                                                                           
          dD_pol = dds.deformation_density('pol')
          dD     = dds.deformation_density('fqm')

          self._dD_indA_set_ref.append(dD_pol[:self._nbf_A,:self._nbf_A])
          self._dD_indB_set_ref.append(dD_pol[self._nbf_A:,self._nbf_A:])
          self._dD_ctAB_set_ref.append(dD_pol[:self._nbf_A,self._nbf_A:])

          #F_A, F_B = self._compute_efield()
          #self._F_A_set.append(F_A)
          #self._F_B_set.append(F_B)

          self._T_set.append(self._compute_T())

          #DIP = mints.ao_dipole()
          #self._DIP_set.append(DIP)

      self._T_set = numpy.array(self._T_set)
      #self._F_A_set = numpy.array(self._F_A_set)
      #self._F_B_set = numpy.array(self._F_B_set)


  def _compute_dms_ct(self):
      print(" * Computing DMS CT")
      DIM = int(self._nbf_A * (self._nbf_A + 1) / 2)
      self._h_ct = numpy.zeros((DIM, DIM))

      DT = numpy.zeros((self._nbf_A, self._nbf_A))
      TT = numpy.zeros((self._nbf_A, self._nbf_A))
      for n in range(self._nsampl_ct):
          D = self._dD_ctAB_set_ref[n]
          T = self._T_set[n]
          TT += T.T @ T
          DT += D.T @ T
      DT *= -8.0
      TT *= 16.0

      self._g_ct = numpy.zeros(DIM)

      ij = -1
      for i in range(self._nbf_A):
          for j in range(i+1):
              ij += 1
              kl = -1
              self._g_ct[ij] = DT[i,j]
              for k in range(self._nbf_A):
                  for l in range(k+1):
                      kl += 1
                      if j==l: self._h_ct[ij,kl] = TT[i,k]

      Hi = self._invert_hessian(self._h_ct)
      B  = -numpy.dot(Hi, self._g_ct)
      print(((self._h_ct-self._h_ct.T)**2).sum())
      print(B)

      self._dms_ct = numpy.zeros((self._nbf_A, self._nbf_A))
      ij = -1
      for i in range(self._nbf_A):
          for j in range(i+1):
              ij +=1
              v = B[ij]
              self._dms_ct[i,j] = v
              self._dms_ct[j,i] = v

  def check_dms_ct(self):
      print(" * Checking DMS CT")
      error = 0.0
      D = self._wfn_A.Da().to_array(dense=True)

      mints = psi4.core.MintsHelper(self._wfn_A.basisset())
      DIP = mints.ao_dipole()
      for n in range(self._nsampl_ct):
          S      = self._S_AB_set[n]
          dD_ref = self._dD_ctAB_set_ref[n]
          dD_com = 2.0 * D @ S @ self._dms_ct

          rms = self._rms(dD_ref, dD_com)
          error += rms
          av_ref = abs(dD_ref).mean()
          av_com = abs(dD_com).mean()

          # dipole moment due to CT
          mu_ref = numpy.array([ 2.0*(dD_ref @ x.to_array(dense=True)).trace() for x in DIP])
          mu_com = numpy.array([ 2.0*(dD_com @ x.to_array(dense=True)).trace() for x in DIP])
          a_mu_ref = numpy.linalg.norm(mu_ref)
          a_mu_com = numpy.linalg.norm(mu_com)
          print(" * - Sample %d RMS= %14.5f  AV_C= %14.5f AV_R= %14.5f  M_C= %14.4f M_R=%14.4f" \
                      % (n+1,rms, av_com, av_ref, a_mu_com, a_mu_ref))
      error/=self._nsampl_ct
      print(" * - Average RMS= %14.5f" % error)
      return error

class Asymmetric_Translation_SameMoleculeDMSFit(Translation_SameMoleculeDMSFit):
  def __init__(self, mol_A, method, nsampl_pol, nsampl_ct):
      Translation_SameMoleculeDMSFit.__init__(self, mol_A, method, nsampl_pol, nsampl_ct)
      self._R_set = []

  def _compute_R(self):
      D = self._W # self._wfn_A.Da().to_array(dense=True)
      return D @ self._S

  def _compute_T(self):
      D = self._W # self._wfn_A.Da().to_array(dense=True)
      return self._S @ D



  def _prepare_for_ct(self):
      "Compute electric fields, CT kernels and reference deformation density matrices"
      print(" * Computing DDS for Each Sample")
      s = self._wfn_A.S().to_array(dense=True)
      self._bfs_A = self._wfn_A.basisset()
      self._W = self._wfn_A.Da().to_array(dense=True)
      self._W = s @ self._W @ s
      #self._W = matrix_power(self._W, 0.5)
      for n in range(self._nsampl_ct):
          psi4.core.clean()
          print(" * - Sample %d" % self._i)
          aggr= self._construct_aggregate()
          dds = DensityDecomposition(aggr, method=self._method, acbs=False, jk_type='direct', 
                                           no_cutoff=0.000, xc_scale=1.0, l_dds=False, n_eps=5.0E-5, cc_relax=True,
                                           verbose=False) 
          dds.compute(polar_approx=False)
          self._S = dds.matrix["sqm"][:self._nbf_A,self._nbf_A:] # S_AB
          self._S_AB_set.append(self._S.copy())
          #self._bfs_A = dds.data["wfn"][0].basisset()
          #self._bfs_B = dds.data["wfn"][1].basisset()
          mol_B = aggr.extract_subsets(2)
          self._bfs_B = psi4.core.BasisSet.build(mol_B, "BASIS", psi4.core.get_global_option("BASIS"), puream=psi4.core.get_global_option("PUREAM"))
                                                                                                           
          dD_pol = dds.deformation_density('pol')
          dD     = dds.deformation_density('fqm')

          self._dD_indA_set_ref.append(dD_pol[:self._nbf_A,:self._nbf_A])
          self._dD_indB_set_ref.append(dD_pol[self._nbf_A:,self._nbf_A:])
          self._dD_ctAB_set_ref.append(dD_pol[:self._nbf_A,self._nbf_A:])

          #F_A, F_B = self._compute_efield()
          #self._F_A_set.append(F_A)
          #self._F_B_set.append(F_B)

          self._T_set.append(self._compute_T())
          self._R_set.append(self._compute_R())

          #DIP = mints.ao_dipole()
          #self._DIP_set.append(DIP)

      self._T_set = numpy.array(self._T_set)
      self._R_set = numpy.array(self._R_set)
      #self._F_A_set = numpy.array(self._F_A_set)
      #self._F_B_set = numpy.array(self._F_B_set)


  def _compute_dms_ct(self):
      print(" * Computing DMS CT")
      DIM = self._nbf_A**2
      self._h_ct = numpy.zeros((DIM, DIM))

      DT = numpy.zeros((self._nbf_A, self._nbf_A))
      TT = numpy.zeros((self._nbf_A, self._nbf_A))
      RD = numpy.zeros((self._nbf_A, self._nbf_A))
      RR = numpy.zeros((self._nbf_A, self._nbf_A))

      for n in range(self._nsampl_ct):
          D = self._dD_ctAB_set_ref[n]
          T = self._T_set[n]
          R = self._R_set[n]
          TT += T @ T.T
          DT += D @ T.T
          RD += R.T @ D
          RR += R.T @ R 

      p = numpy.einsum("nac,nbd->abcd", self._R_set, self._T_set) * 2.0
      q = numpy.einsum("nca,ndb->abcd", self._R_set, self._T_set) * 2.0
      #p = numpy.einsum("nad,nbc->abcd", self._R_set, self._T_set) * 2.0
      #q = numpy.einsum("ncb,nda->abcd", self._R_set, self._T_set) * 2.0

      DT *= -2.0
      RD *= -2.0
      TT *=  2.0
      RR *=  2.0

      self._g_ct = numpy.zeros(DIM)

      for i in range(self._nbf_A):
          for j in range(self._nbf_A):
              ij = i*self._nbf_A + j
              self._g_ct[ij] = DT[i,j] + RD[i,j]
              #self._g_ct[ij] = DT[i,j] + RD[j,i]
              for k in range(self._nbf_A):
                  for l in range(self._nbf_A):
                      kl = k*self._nbf_A + l
                      v = 0.0
                      #for N in range(self._nsampl_ct):
                      #    v += self._R_set[N][i,k] * self._T_set[N][j,l]
                      #    v += self._R_set[N][k,i] * self._T_set[N][l,j]
                      #v *= 2.0
                      v = p[i,j,k,l] + q[i,j,k,l]
                      if i==k: v += TT[j,l] #+ RR[j,l]
                      if j==l: v += RR[k,i] 
                      
                      self._h_ct[ij,kl] = v
                      

      Hi = self._invert_hessian(self._h_ct)
      B  = -numpy.dot(Hi, self._g_ct)

      self._dms_ct = numpy.zeros((self._nbf_A, self._nbf_A))
      for i in range(self._nbf_A):
          for j in range(self._nbf_A):
              ij = i*self._nbf_A + j
              self._dms_ct[i,j] = B[ij]
      print(" * DMS CT: ")
      print(self._dms_ct.round(2))

  def check_dms_ct(self):
      print(" * Checking DMS CT")
      error = 0.0
      #D = self._wfn_A.Da().to_array(dense=True)
      D = self._W

      mints = psi4.core.MintsHelper(self._wfn_A.basisset())
      DIP = mints.ao_dipole()
      for n in range(self._nsampl_ct):
          S      = self._S_AB_set[n]
          dD_ref = self._dD_ctAB_set_ref[n]
          dD_com = self._dms_ct @ S @ D + D @ S @ self._dms_ct #.T

          rms = self._rms(dD_ref, dD_com)
          error += rms
          av_ref = abs(dD_ref).mean()
          av_com = abs(dD_com).mean()
          print(dD_ref[0])
          print(dD_com[0])

          # dipole moment due to CT
          mu_ref = numpy.array([ 2.0*(dD_ref @ x.to_array(dense=True)).trace() for x in DIP])
          mu_com = numpy.array([ 2.0*(dD_com @ x.to_array(dense=True)).trace() for x in DIP])
          a_mu_ref = numpy.linalg.norm(mu_ref)
          a_mu_com = numpy.linalg.norm(mu_com)
          print(" * - Sample %d RMS= %14.5f  AV_C= %14.5f AV_R= %14.5f  M_C= %14.4f M_R=%14.4f" \
                      % (n+1,rms, av_com, av_ref, a_mu_com, a_mu_ref))
      error/=self._nsampl_ct
      print(" * - Average RMS= %14.5f" % error)
      return error



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
from abc import ABC, abstractmethod
from .partitioning import DensityDecomposition

__all__ = ["DMSFit"]


class DMSFit(ABC):
  def __init__(mol_A, mol_B, method, nsampl_pol, nsampl_ct):
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
      self._dms_pol_A, self._dms_pol_B = None, None
      self._dms_ct_A, self._dms_ct_B = None, None
      self._nbf_A, self._nbf_B = None, None

  @classmethod
  def create(cls, mol_A, mol_B = None, type="transl", nsampl_pol=100, nsampl_ct=100, method='scf'):
      if type.lower().startswith() == "tran": return Translation_SameMoleculeDMSFit(mol_A, method, nsampl_pol, nsampl_ct)
      else: raise NotImplementedError("This type of DMS fitting is not implemented yet")

  def run(self):
      self._compute_wfn()
     #self._compute_dms_induction()
      self._prepare_for_ct()
      self._compute_dms_ct()

  # --- protected --- #

  def _compute_wfn(self):
      _g_A, self._wfn_A = psi4.gradient(self._method, molecule=self._mol_A, return_wfn=True)
      self._nbf_A = self._wfn_A.basisset().nbf()
      if self._mol_B is not None:
         _g_B, self._wfn_B = psi4.gradient(self._method, molecule=self._mol_B, return_wfn=True)
         self._nbf_B = self._wfn_B.basisset().nbf()

  def _compute_dms_induction(self):
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
     hi = numpy.linalg.inv(hi)
     I = numpy.dot(hi, h).diagonal().sum()
     d = I - len(hi)
     if abs(d) > 0.0001: raise ValueError("Hessian is problemmatic!")
     return hi      

class SameMoleculeDMSFit(DMSFit):
  def __init__(mol_A, method, nsampl_pol, nsampl_ct):
      DMSFit.__init__(self, mol_A, None, method, nsampl_pol, nsampl_ct)
      self._g = None
      self._h = None
      self._dD_pol_set_ref = []
      self._dD_pol_set_com = []

  def _prepare_for_ct(self):
      for n in range(self._nsampl_ct):
          aggr= self._construct_aggregate()
          dds = DensityDecomposition(aggr, method=self._method, acbs=True, jk_type='direct', 
                                           no_cutoff=0.000, xc_scale=1.0, l_dds=False, n_eps=5.0E-5, cc_relax=True,
                                           verbose=False) 
          dds.compute(polar_approx=False)
                                                                                                            
          dD_pol = solver.deformation_density('pol') #TODO
          dD     = solver.deformation_density('fqm')

          self._dD_pol_set_ref.append(dD_pol[:self._nbf_A,self._nbf_A:])
          self._dD_pol_set_com.append(numpy.zeros((self._nbf_A, self._nbf_B)))

          
  def _compute_dms_ct(self):
      Hi = self._invert_hessian(self._h)
      self._dms_ct_A = -numpy.dot(self._g, Hi)

  @abstractmethod
  def _compute_gradient(self): pass
  @abstractmethod
  def _compute_hessian(self): pass



class Translation_SameMoleculeDMSFit(SameMoleculeDMSFit):
  def __init__(mol_A, method, nsampl_pol, nsampl_ct):
      SameMoleculeDMSFit.__init__(mol_A, method, wfn, nsampl_pol, nsampl_ct)
      self._i = 0

  def _construct_aggregate(self): raise NotImplementedError
      #self._i += 1
      #TODO
  def _compute_gradient(self): raise NotImplementedError
      #TODO
  def _compute_hessian(self): raise NotImplementedError
      #TODO



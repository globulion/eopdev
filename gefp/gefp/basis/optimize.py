#!/usr/bin/python3
"""
 Auxiliary Basis Set Optimization Library.
 
 The auxiliary basis sets for generalized density fitting (GDF)
 are here optimized.
"""
from abc import ABC, abstractmethod
import re
import psi4
import numpy
import scipy.optimize
import oepdev

__all__ = ["DFBasis", "DFBasisOptimizer", "OEP"]

def make_bastempl(templ, param): return templ % tuple(param)

def oepfitbasis(mol, role='ORBITAL'):
    "Requires setting 'oepfitbasis.templ' and 'oepfitbasis.param' inside bas prior to executing bas"
    basstrings = {}
    mol.set_basis_all_atoms("oepfit", role=role)
    basstrings['oepfit'] = make_bastempl(oepfitbasis.templ, oepfitbasis.param)
    return basstrings

def stripComments(code):
    "Remove comments from the string. Delimiter: #"
    code = str(code)
    return re.sub(r'(?m)^ *#.*\n?', '', code)


# -------------------------------------------------------------------------------------------------- #

class DFBasis:
  """
 Basis set object to be optimized
"""
  def __init__(self, mol, templ_file='templ.dat', param_file='param.dat'):
      "Initialize"
      # bound molecule to this basis set
      self.mol= mol
      # initial parameters
      #self.param_0 = self._read_param(param_file)
      self.param = numpy.mafromtxt(param_file).data
      # template for Psi4 input
      self.templ = open(templ_file).read()
      # current basis set
      self.basis = self.basisset()

  def basisset(self, param=None):
      "Construct basis set from optimization parameters"
      if param is None: param = self.param
      oepfitbasis.param = param
      oepfitbasis.templ = self.templ
      basis = psi4.core.BasisSet.build(self.mol, 'BASIS', oepfitbasis)
      self.basis = basis
      return basis

  def print(self, param=None):
      "Print basis set in Psi4 format"
      if param is None: param = self.param
      basis_str = make_bastempl(self.templ, param)
      return basis_str

  def save(self, out='oepfit.gbs', param=None):
      "Save basis set in Psi4 format to a file"
      nm = out[:-4]
      o  = open(out, 'w')
      o.write("[ %s ]\n" % nm)
      o.write(self.print(param))
      o.close()
      return
  
  # - protected 

  def _read_param(self, param_file):
      "Read parameters to optimize from the input file. Comments are allowed."
      t = open(param_file, 'r')
      line = stripComments(t.readline())
      print(line)
      data = []
      while line:
          data += line.split() 
          line = stripComments(t.readline())
      t.close()
      return numpy.array(data, dtype=numpy.float64)     

# ------------------------------------------------------------------------ #
     
class OEP(ABC):
  """
 OEP object that defines the V matrix necessary for GDF
"""
  def __init__(self, wfn, dfbasis):
      "Initialize"
      # wavefunction
      self.wfn = wfn
      # DF basis to optimize
      self.dfbasis = dfbasis
      # Integral calculator
      self.mints = psi4.core.MintsHelper(wfn.basisset())
      # Testing (T) basis (set to be RI basis)
      self.basis_test = wfn.get_basisset("BASIS_DF_INT")
      # Primary (P) basis
      self.basis_prim = wfn.basisset()
      # Auxiliary (A) basis
      self.basis_aux = dfbasis.basis
      # ERI object (PP|PT) 
      self.eri_pppt = self.mints.ao_eri(self.basis_prim, self.basis_prim, self.basis_prim, self.basis_test)
      # V matrix
      self.V   = self._compute_V()
      # Clear memory
      del self.eri_pppt
      #
      super(OEP, self).__init__()

  @classmethod
  def create(cls, name, wfn, dfbasis):
      "Create OEP instance"
      if   name.lower() == 'pauli': oep = OEP_Pauli(wfn, dfbasis)
      elif name.lower() == 'ct'   : oep = OEP_CT   (wfn, dfbasis)
      else: raise ValueError("Not existent OEP chosen")
      return oep
  
  @abstractmethod
  def _compute_V(self):
      "Compute matrix representation of OEP"
      pass

class OEP_FockLike(OEP):
  """
 OEP for S1 term in Murrell et al.'s theory of Pauli repulsion
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      super(OEP_FockLike, self).__init__(wfn, dfbasis)

  @abstractmethod
  def _compute_Ca(self): pass

  # - implementation - #
  def _compute_V(self):
      "Compute Fock-Like OEP Matrix"
      Ca= self._compute_Ca()
      Da= self.wfn.Da().to_array(dense=True)
      Vao= self.mints.ao_potential(self.basis_test, self.basis_prim).to_array(dense=True)
      V  = numpy.dot(Vao, Ca)
      #
      V += 2.0 * numpy.einsum("bi,mn,mnba->ai", Ca, Da, self.eri_pppt)
      V -=       numpy.einsum("mi,nb,mnba->ai", Ca, Da, self.eri_pppt)
      return V


class OEP_Pauli(OEP_FockLike):
  """
 OEP for S1 term in Murrell et al.'s theory of Pauli repulsion
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      super(OEP_Pauli, self).__init__(wfn, dfbasis)

  # - implementation - #
  def _compute_Ca(self): return self.wfn.Ca_subset("AO", "OCC").to_array(dense=True)


class OEP_CT(OEP_FockLike):
  """
 OEP for Group-(i) term of Otto-Ladik's theory of Charge-Transfer Energy
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      super(OEP_CT, self).__init__(wfn, dfbasis)

  # - implementation - #
  def _compute_Ca(self): return self.wfn.Ca_subset("AO", "VIR").to_array(dense=True)


# -------------------------------------------------------------------------------------------- #

class DFBasisOptimizer:
  """
 Method that optimizes DF basis set.
"""
  def __init__(self, oep):
      "Initialize"
      self.oep = oep
      self.basis_fit = oep.dfbasis
      self.param = oep.dfbasis.param
      self.n_param = len(oep.dfbasis.param)
  
  # ---> public interface <--- #

  def fit(self):
      "Perform fitting of basis set exponents"
      success = self._fit()
      return success

  def compute_error(self, basis):
      "Compute error for a given basis"
      S = self.oep.mints.ao_overlap(self.oep.basis_test, basis).to_array(dense=True)
      V = psi4.core.Matrix.from_array(self.oep.V)
      gdf = oepdev.GeneralizedDensityFit.build_double(basis, self.oep.basis_test, V)
      gdf.compute()
      G = gdf.G().to_array(dense=True)
      er= V - numpy.dot(S, G)
      Z = (er*er).sum()
      #Z = numpy.sqrt(Z/er.size)
      return Z


  # ---> protected interface <--- #

  def _fit(self):
      "The actual fitting procedure"
      param_0 = self.basis_fit.param
      bounds = [(0.01, 10000.0) for x in range(self.n_param)]
      options = {"disp": True, "maxiter": 300, "ftol": 1.0e-9, "iprint": 4}
      res = scipy.optimize.minimize(self._objective_function, param_0, args=(), tol=1.0e-9, method='slsqp',
                                    options=options, bounds=bounds)
      param = res.x
      self.basis_fit = self.oep.dfbasis.basisset(param)
      #
      print("          Initial           Fitted")
      for i in range(self.n_param):
          print("%15.3f %15.3f" % (param_0[i], param[i]))
      self.oep.dfbasis.param = param
      return res.success
     
  def _objective_function(self, param):
      "Objective function for auxiliary basis refinement"
      basis_aux = self.oep.dfbasis.basisset(param) 
      Z = self.compute_error(basis_aux)
      return Z

#!/usr/bin/python3
"""
 Auxiliary Basis Set Optimization Library.
 
 The auxiliary basis sets for generalized density fitting (GDF)
 are here optimized.
"""
from abc import ABC, abstractmethod
import re
import os
import time
import psi4
import numpy
import scipy.optimize
import oepdev
from . import edf
from . import _util

__all__ = ["DFBasis", "oep_ao_basis_set_optimizer", "DFBasisOptimizer", "OEP"]

def make_bastempl(templ, param): return templ % tuple(param)

def oepfitbasis(mol, role='ORBITAL'):
    "Requires setting 'oepfitbasis.templ' and 'oepfitbasis.param' inside bas prior to executing bas"
    basstrings = {}
    mol.set_basis_all_atoms("oepfit", role=role)
    basstrings['oepfit'] = make_bastempl(oepfitbasis.templ, oepfitbasis.param)
    return basstrings

def removeComments(string):
    "Remove comments from the string. Delimiter: #"
    string = re.sub(re.compile("#.*?\n" ) ,"" ,string)
    return string


# -------------------------------------------------------------------------------------------------- #

class DFBasis:
  """
 Basis set object to be optimized.

 Notes:

 o Default bounds can be modified by resetting static variables
     DFBasis.exp_lower_bound   
     DFBasis.exp_upper_bound
     DFBasis.ctr_lower_bound
     DFBasis.ctr_upper_bound
   prior to calling DFBasis if not using the driver.gdf_basis_optimizer.
"""

  # defaults
  exp_lower_bound = 0.01
  exp_upper_bound = 10000.0
  ctr_lower_bound =-2.0
  ctr_upper_bound = 2.0

  def __init__(self, mol, 
                     templ_file='templ.dat', param_file='param.dat', bounds_file=None, 
                     constraints=()):
      "Initialize"
      # ensure physical values are set
      assert(self.exp_lower_bound > 0.0), "Orbital exponents cannot be negative or zero!"

      # bound molecule to this basis set
      self.mol= mol

      # template for Psi4 input
      self.templ = open(templ_file).read()

      # initial parameters input
      self.param = self._read_input(param_file, dtype=numpy.float64)
      self.n_param = len(self.param)

      # parameter bounds: if not provided assume only exponents
      if bounds_file is None: 
         self.bounds = [(self.exp_lower_bound, self.exp_upper_bound) for x in range(self.n_param)]
      else: 
         self.bounds= self._make_bounds(self._read_input(bounds_file, dtype=str))

      # parameter constraints
      self.constraints = constraints

      # current basis set
      self.basis = self.basisset()

  def basisset(self, param=None):
      "Construct basis set from optimization parameters"
      if param is None: param = self.param
      oepfitbasis.param = param
      oepfitbasis.templ = self.templ
      basis = psi4.core.BasisSet.build(self.mol, 'BASIS', oepfitbasis, quiet=True)
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

  def __repr__(self):
      "Print the Basis Set in Psi4 Format"
      log = self.print()
      return str(log)
  
  # - protected 

  def _read_input(self, input_file, dtype):
      "Read parameters to optimize from the input file. Comments are allowed."
      t = open(input_file, 'r')
      data = removeComments(t.read()).split()
      t.close()
      return numpy.array(data, dtype=dtype)

  def _make_bounds(self, bounds_codes):
      "Make list of bounds appropriate for optimization"
      bounds = []
      for item in bounds_codes:
          print(item,bounds_codes)
          item = item.lower()
          if item.startswith('e'):
               if len(item)==1:
                   bound = (self.exp_lower_bound, self.exp_upper_bound)
               else:
                   bound = tuple( float(x) for x in item[2:].split(',') )
               #bound = (self.exp_lower_bound, self.exp_upper_bound) if len(item) == 1 \
               #  else tuple( float(x) for x in item[2:].split(',') )
          elif item.startswith('c'):
               if len(item)==1:
                   bound = (self.ctr_lower_bound, self.ctr_upper_bound)
               else:
                   bound = tuple( float(x) for x in item[2:].split(',') )
               #bound = (self.ctr_lower_bound, self.ctr_upper_bound) if len(item) == 1 \
               #  else tuple( float(x) for x in item[2:].split(',') )
          else: raise(ValueError, 
                "Wrong input in bounds `%s`. Use either one syntax: `X` or `X:min,max` where X = E or C." % item)
          bounds.append(bound)
      return bounds
  
# -------------------------------------------------------------------------------------------- #

def oep_ao_basis_set_optimizer(wfn, interm, 
                  test=None, exemplary=None, target="OCC", cpp=False, more_info=False, 
                  templ_file='templ.dat', param_file='param.dat', bound_file=None, constraints=(), outname='basis.gbs',
                  opt_global=False, use_standardized_inputs=True):
    """
 Method that optimizes DF basis set.
 This is currently the state-of-the-art and recommended.
"""

    # set up
    if test is None:
       test = wfn.basisset()
    if exemplary is None:
       exemplary = wfn.basisset()

    do_quambo = bool(psi4.core.get_global_option("EFP2_WITH_VVO"))
    if do_quambo:
       quambo_solver = oepdev.QUAMBO.build(wfn, bool(psi4.core.get_global_option("QUAMBO_ACBS")))
       quambo_solver.compute()
       Ca = quambo_solver.Ca_subset("AO", target).to_array(dense=True)        # Here are the target orbitals
       localize = False
    else:
       localize = bool(psi4.core.get_global_option("OEPDEV_LOCALIZE"))
       if localize:
          localizer = psi4.core.Localizer.build("BOYS", wfn.basisset(), wfn.Ca_subset("AO", target))
          localizer.localize()
          Ca = localizer.L.to_array(dense=True)
       else:
          Ca = wfn.Ca_subset("AO", target).to_array(dense=True)                  # Here are the target orbitals
          

    if opt_global:
       raise NotImplementedError #TODO
    if use_standardized_inputs:
       raise NotImplementedError #TODO
   

    print("\n ===> Auxiliary Basis Set Optimization Routine <===\n")

    # Notation:
    # T - target (MO)
    # p - primary (AO)
    # i - intermediate (AO)
    # m - minimal auxiliary (AO) - to be optimized
    # M - minimal auxiliary (MO)
    # t - test (AO)
    # e - example AO basis


    primary = wfn.basisset()

    no = Ca.shape[1]                                                       # number of target orbitals
    np = primary.nbf()

    Da   = wfn.Da().to_array(dense=True)                                   # (p, p)
    S_pp = wfn.S().to_array(dense=True)                                    # (p, p)

    mints = psi4.core.MintsHelper(interm)
    S_ii = mints.ao_overlap().to_array(dense=True)                         # (i, i)
    V_iT = edf.compute_v(Ca, Da, prim=primary, left_axis=interm)           # (i, T)

    G = numpy.linalg.inv(S_ii) @ V_iT                                      # (i, T)
    T = edf.find_aux_mo_mini(G, S_ii, I=None)                              # (i, M)

    GM= T.T @ S_ii @ G                                                     # (M, T) #OK

    g = T @ GM # = G                                                       # (i, T) #OK

   # these are true:
   # T.T @ S_ii @ T ===> identity matrix
   # g = G

    # optimize aux_mini (AO)
    dfbasis = DFBasis(wfn.molecule(), templ_file, param_file, bounds_file=bound_file, constraints=constraints)
    param_0 = dfbasis.param

    TIME = -time.time()
    dfbasis = edf.optimize_ao_mini(T, interm, dfbasis, cpp)
    TIME += time.time()

    auxiliary = dfbasis.basisset()
    param = dfbasis.param
    dfbasis.save(outname)

    z_opt = -edf.obj_numpy(param, T, interm, dfbasis)
    z_max = T.shape[1]

    mints = psi4.core.MintsHelper(interm)
    s_im  = mints.ao_overlap(interm, auxiliary).to_array(dense=True)
    s_mm  = mints.ao_overlap(auxiliary, auxiliary).to_array(dense=True)
    o_mi  = edf.projection(T, s_im, s_mm)
    Tm    = edf.projected_t(T, s_im, s_mm)

    print(" Optimized Overlaps:")
    for i in range(len(o_mi.diagonal())):
        print(" %3d %14.7f" % (i+1, o_mi.diagonal()[i])) 

   #Gm = Tm @ numpy.linalg.inv(o_mi) @ T.T @ S_ii @ G # approximate Gm -> this is just guess
   #Gm = Tm @ Tm.T @ s_im.T @ G                       # approximate Gm -> this is derived from 1 = |m) Tm Tm.T (m|
    Gm = Tm @ numpy.linalg.inv(o_mi) @ Tm.T @ s_im.T @ G   # approximate Gm -> above but corrected by the overlap in B

    print(" Optimized AO Basis Set in Psi4 Format:")
    print(dfbasis)
    print(" Inital and final parameters in working order:")
    #
    print("          Initial           Fitted")
    for i in range(len(param)):
        print("%15.3f %15.3f" % (param_0[i], param[i]))


    print(" Running Tests")
    # test 1 
    print(" Test 1: Computed from approximate Gm")
    mints = psi4.core.MintsHelper(test)
    V_= edf.compute_v(Ca, Da, prim=primary, left_axis=test)                  # (t, T)
    S_= mints.ao_overlap(test, auxiliary).to_array(dense=True)           # (t, m)
    Vtest_ = S_ @ Gm                                                     # (t, m) @ (m, T) -> (t, T)
    if more_info: _util.COMPARE(V_, Vtest_, 0)
    Err = Vtest_ - V_
    Err = Err**2
    Err = numpy.sqrt(Err.sum())
    Err_out = Err
    S__ = S_
    Err_1 = Err
    print(" Total error 1 = %14.5f\n" % Err)

    # test 2
    print(" Test 2: Computed from exact GB")
    mints = psi4.core.MintsHelper(test)
    S_= mints.ao_overlap(test, interm).to_array(dense=True) @ T          # (t, i)
    Vtest_ = S_ @ GM                                                     # (t, m) @ (m, T) -> (t, T)
    if more_info: _util.COMPARE(V_, Vtest_, 0)
    Err = Vtest_ - V_
    Err = Err**2
    Err = numpy.sqrt(Err.sum())
    Err_2 = Err
    print(" Total error 2 = %14.5f\n" % Err)

    # test 3 - EDF-2 with auxiliary basis
    print(" Test 3: Computed from optimal Gm from EDF-2")
    V__= edf.compute_v(Ca, Da, prim=primary, left_axis=interm)                  
    V = psi4.core.Matrix.from_array(V__)
    gdf = oepdev.GeneralizedDensityFit.build_double(auxiliary, interm, V)
    gdf.compute()
    G3 = gdf.G().to_array(dense=True)                                   # (m, T) from EDF-2
    Vtest2_ = S__ @ G3
    if more_info: _util.COMPARE(V_, Vtest2_, 0)
    Err = Vtest2_ - V_
    Err = Err**2
    Err = numpy.sqrt(Err.sum())
    Err_3 = Err
    print(" Total error 3 = %14.5f\n" % Err)

    # test 4 - EDF-2 with example basis
    print(" Test 4: Computed from exemplary Gm from EDF-2")
    mints = psi4.core.MintsHelper(test)
    S_= mints.ao_overlap(test, exemplary).to_array(dense=True)              # (t, e)
    gdf = oepdev.GeneralizedDensityFit.build_double(exemplary, interm, V)
    gdf.compute()
    G4 = gdf.G().to_array(dense=True)                                   # (e, T) from EDF-2
    Vtest3_ = S_ @ G4
    if more_info: _util.COMPARE(V_, Vtest3_, 0)
    Err = Vtest3_ - V_
    Err = Err**2
    Err = numpy.sqrt(Err.sum())
    Err_4 = Err
    print(" Total error 4 = %14.5f\n" % Err)

    print(" Summary:")
    print(" ========")
    print()
    print(" AO Basis Sets                               Name    Nbf.")
    print(" -------------")
    print(" Primary           %30s  %4d" % (primary  .name(), primary  .nbf()))
    print(" Intermediate      %30s  %4d" % (interm   .name(), interm   .nbf()))
    print(" Test              %30s  %4d" % (test     .name(), test     .nbf()))
    print(" Example           %30s  %4d" % (exemplary.name(), exemplary.nbf()))
    print(" Optimized         %30s  %4d" % (auxiliary.name(), auxiliary.nbf()))
    print()
    print(" Settings")
    print(" --------")
    print(" Target              %s" % target)
    print(" Use QUAMBO?         %r" % do_quambo)
    print(" Localized?          %r" % localize)
    print(" Global opt?         %r" % opt_global)
    print()
    print(" Objective Function")
    print(" ------------------")
    print(" Z_max = %14.8f" %   z_max)
    print(" Z_opt = %14.8f" %   z_opt)
    print(" Error = %14.8f" % ( z_max - z_opt))
    print(" NrmEr = %14.8f" % ((z_max - z_opt)/z_max))
    print()
    print(" Cumulative Errors")
    print(" -----------------")
    print()
    print(" Fit(EDF-2)   Fit(Approx)  Fit(Exact)   Fit(Example)")
    print(" %8.5f     %8.5f     %8.5f     %8.5f" % (Err_3, Err_1, Err_2, Err_4))
    print()
    print(" Overall AO Basis Set Optimization Time = %f sec" % TIME)
    print()
    print(" Basis set saved to file < %s >" % os.path.abspath(outname))
    print(" AO Basis Optimization Finished Successfully.")
    print()

    return dfbasis, Err_out, Gm

# ------------------------------------------------------------------------ #
     
class OEP(ABC):
  """
 OEP object that defines the V matrix necessary for GDF.
"""
  read_vints = True
  def __init__(self, wfn, dfbasis):
      "Initialize"
      # wavefunction
      self.wfn = wfn
      # DF basis to optimize
      self.dfbasis = dfbasis
      # Integral calculator
      #self.mints = psi4.core.MintsHelper(wfn.basisset())
      # Testing (T) basis (set to be RI basis)
      self.basis_test = wfn.get_basisset("BASIS_DF_INT")
      self._nt = self.basis_test.nbf()
      # Primary (P) basis
      self.basis_prim = wfn.basisset()
      # Auxiliary (A) basis
      self.basis_aux = dfbasis.basis
      # ERI object (PP|PT) 
      #self.eri_pppt = numpy.asarray(self.mints.ao_eri(self.basis_prim, self.basis_prim, self.basis_prim, self.basis_test))
      # V matrix
      self._vints_file_name = 'vints.dat'
      self.V = None
      # Clear memory
      #del self.eri_pppt
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

  def compute(self):
      if self.V is None:
         self.V = self._compute_V()
      else:
         self.V   = self._extract_V()

  def compute_and_save_V(self, name='vints.dat'):
      self._vints_file_name = name
      V = self._compute_V()
      V.tofile(name)
      self.V = V

  def _extract_V(self):
      if OEP.read_vints: return self._read_V()
      else: return self._compute_V()

  def _read_V(self):
      return numpy.fromfile(self._vints_file_name).reshape(self._nt, self._ni)
      

class OEP_FockLike(OEP):
  """
 OEP for S1 term in Murrell et al.'s theory of Pauli repulsion
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      self._Da= wfn.Da().to_array(dense=True)
      super(OEP_FockLike, self).__init__(wfn, dfbasis)


  # - implementation - #
  def _compute_V(self):
      "Compute Fock-Like OEP Matrix"
      Ca= self._Ca_xxx
      Da= self._Da

      # ---> Vao <--- #
      # one-electron part
      mints = psi4.core.MintsHelper(self.wfn.basisset())
      V_nuc= mints.ao_potential(self.basis_test, self.basis_prim).to_array(dense=True)
      V  = numpy.dot(V_nuc, Ca)

      # two-electron part
      f_pppt = psi4.core.IntegralFactory(self.basis_prim, self.basis_prim, self.basis_prim, self.basis_test)
      Ca_ = psi4.core.Matrix.from_array(Ca.copy(), "")
      Da_ = psi4.core.Matrix.from_array(Da.copy(), "")
      V_2el = oepdev.calculate_OEP_basisopt_V(self._nt, f_pppt, Ca_, Da_)
      V += V_2el.to_array(dense=True)
      psi4.core.clean()
      #
      #V += 2.0 * numpy.einsum("bi,mn,mnba->ai", Ca, Da, self.eri_pppt)
      #V -=       numpy.einsum("mi,nb,mnba->ai", Ca, Da, self.eri_pppt)
      return V


class OEP_Pauli(OEP_FockLike):
  """
 OEP for S1 term in Murrell et al.'s theory of Pauli repulsion
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      # - implementation - #
      self._Ca_xxx = wfn.Ca_subset("AO", "OCC").to_array(dense=True)
      self._ni = self._Ca_xxx.shape[1]


      super(OEP_Pauli, self).__init__(wfn, dfbasis)


class OEP_CT(OEP_FockLike):
  """
 OEP for Group-(i) term of Otto-Ladik's theory of Charge-Transfer Energy
 """
  def __init__(self, wfn, dfbasis):
      "Initialize"
      # - implementation - #
      self._Ca_xxx = wfn.Ca_subset("AO", "VIR").to_array(dense=True)
      self._ni = self._Ca_xxx.shape[1]


      super(OEP_CT, self).__init__(wfn, dfbasis)


# -------------------------------------------------------------------------------------------- #

class DFBasisOptimizer:
  """
 Method that optimizes DF basis set.
 This is currently not recommended.
"""
  def __init__(self, oep):
      "Initialize"
      self.oep = oep
      self.basis_fit = oep.dfbasis
      self.param = oep.dfbasis.param
     #self.mints = psi4.core.MintsHelper(oep.basis_test)
  
  # ---> public interface <--- #

  def fit(self, maxiter=1000, tolerance=1e-9, method='slsqp', opt_global=False, temperature=500, stepsize=500,
                take_step=None, accept_test=None):
      "Perform fitting of basis set exponents"
      success = self._fit(maxiter, tolerance, method, opt_global, temperature, stepsize, take_step, accept_test)
      return success

  def compute_error(self, basis, rms=False):
      "Compute error for a given basis"
      mints = psi4.core.MintsHelper(basis)
      S = mints.ao_overlap(self.oep.basis_test, basis).to_array(dense=True)
      V = psi4.core.Matrix.from_array(self.oep.V)
      gdf = oepdev.GeneralizedDensityFit.build_double(basis, self.oep.basis_test, V)
      gdf.compute()
      G = gdf.G().to_array(dense=True)
      er= V - numpy.dot(S, G)
      Z = (er*er).sum()
      if rms: Z = numpy.sqrt(Z/er.size)
      return Z


  # ---> protected interface <--- #

  def _fit(self, maxiter, tolerance, method='slsqp', opt_global=False, temperature=500, stepsize=500, 
                 take_step=None, accept_test=None):
      "The actual fitting procedure"
      param_0 = self.basis_fit.param
      options = {"disp": True, "maxiter": maxiter, "ftol": tolerance, "iprint": 4}
      if not opt_global:
         res = scipy.optimize.minimize(self._objective_function, param_0, args=(), tol=tolerance, method=method,
                     options=options, bounds=self.basis_fit.bounds, constraints=self.basis_fit.constraints)
         success = res.success
      else:
         res = scipy.optimize.basinhopping(self._objective_function, param_0, niter=30, 
                                    T=temperature, stepsize=stepsize, 
                                    minimizer_kwargs={"method": method, "options": options, 
                                        "bounds": self.basis_fit.bounds}, 
                                    callback=None, interval=50, disp=True, niter_success=None,
                                    take_step=take_step, accept_test=accept_test)
         success = True
      param = res.x
      self.oep.dfbasis.basisset(param)
      self.oep.dfbasis.param = param
      #
      print("          Initial           Fitted")
      for i in range(self.basis_fit.n_param):
          print("%15.3f %15.3f" % (param_0[i], param[i]))
      return success # for some reson res.success is not created in global
     
  def _objective_function(self, param):
      "Objective function for auxiliary basis refinement"
      basis_aux = self.oep.dfbasis.basisset(param) 
      Z = self.compute_error(basis_aux, rms=False)
      # setting rms=False results in more sensitive optimization and better final Z
      return Z

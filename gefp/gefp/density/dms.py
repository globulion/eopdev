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
#from ..math.matrix import move_atom_rotate_molecule, rotate_ao_matrix, matrix_power
from gefp.math.matrix import move_atom_rotate_molecule, rotate_ao_matrix, matrix_power


__all__ = ["DMSFit", "DMS"]

PSI4_DRIVER = psi4.energy
numpy.random.seed(0)

class DMS(ABC):
  """
 Density Matrix Susceptibility Tensor. 
 Abstract Base.

 This describes the DMS of Density or Fock matrix as EFP.

 EFPs:
  - nuclear structure
  - unperturbed Density or Fock matrix (alpha or beta)
  - DMS tensors

 Functionalities:
  - rotation, translation and superimposition
  - read/write capabilities based on FRG format (Solvshift)
"""
  def __init__(self, typ):#OK
      ABC.__init__(self)

      self._bfs      = None           # BasisSet
      self._mol      = None           # Molecule
      self._s1       = None           # DMS parameter space
      self._s2       = None           # DMS parameter space
      self._B        = {}             # DMS tensors
      self._M        = None           # Zeroth-order matrix
      self._type     = None           # Type of DMS (Density or Fock)
      self._type_long= None           # ..as above (long version)
      self._n        = None           # Sizing: number of AOs
      self._N        = None           # Sizing: number of distributed sites (atoms)

      self._init(typ)
      self._available_orders = None   # Implemented orders of DMS tensors

  @classmethod
  def create(cls, dms_type='da', order_type='basic'):#OK
      if   order_type.lower() == 'basic': return       Basic_DMS(dms_type)
      elif order_type.lower() == 'ct-1' : return    LinearCT_DMS(dms_type)
      elif order_type.lower() == 'ct-2' : return QuadraticCT_DMS(dms_type)
      else:
         raise NotImplementedError
      

  def available_orders(self): return self._available_orders #OK

  @property
  def N(self): return self._N
  @property
  def n(self): return self._n
  def B(self, m, n): 
      if (m,n) in self._available_orders:
          return self._B[(m,n)]
      else:
          raise ValueError("DMS Error: B(%i,%i) is not available" % (m,n))

  def set_bfs(self, bfs):#OK
      self._bfs = bfs
      self._mol = bfs.molecule()
      self._n   = bfs.nbf()
      self._N   = self._mol.natom()
  def set_M(self, M): self._M = M.copy()
  def set_s1(self, s): self._s1 = s.copy()
  def set_s2(self, s): self._s2 = s.copy()

  def _init(self, t):#OK
      "Initialize information"
      u = {"d": "Density", "f": "Fock", "a": "Alpha", "b": "Beta"}
      m, s = [x for x in t.lower()]
      self._type      = t.lower()
      self._type_long = u[m] + '-' + u[s]

  def B(self, m, n):#OK
      "Retrieve DMS tensor from parameter vector"
      order = (m, n)
      if order in self._available_orders:
         if order in self._B.keys(): 
            return self._B[order]
         else:
            self._generate_B_from_s(m, n)
      else:
         raise ValueError("This DMS object does not include order (%i,%i)" %(m,n))

  def rotate(self, rot):#OK
      "Rotate DMS tensor. Warning: bfs is not rotated. Only DMS and molecule."
      # rotate OPDM
      self._M, R = rotate_ao_matrix(self._M, rot, self._bfs, return_rot=True, aomo=False)
      # rotate molecule
      angles = scipy.spatial.transform.Rotation.from_dcm(rot).as_euler('zxy', degrees=True)
      move_atom_rotate_molecule(self._mol, angles, t='zxy')
      # rotate basis set
      None # no need for it as for now
      # rotate DMS tensor
      for order in self._available_orders:
          self._rotate(rot, R, order)

  def translate(self, t):#OK
      "Translate DMS object"
      # translate molecule
      xyz = self._mol.geometry().to_array(dense=True)
      self._mol.set_geometry(psi4.core.Matrix.from_array(xyz))
      # translate basis set
      for i in range(self._mol.natom()):
          self._bfs.move_atom(i, t)

  def superimpose(self, xyz, suplist=None):#TODO
      "Superimpose DMS object"
      raise NotImplementedError("Please implement superimposition since it is easy")


  def read(self, out):#TODO
      "Create DMS object from file"
      raise NotImplementedError

  def write(self, out):#TODO
      "Write DMS object to file"
      raise NotImplementedError

  def _generate_B_from_s(self, order):#TODO
      "Retrieve DMS tensor from its parameter vector"
      n  = self._n
      n2 = self._n**2 
      N  = self._N

      # ---> Group Z(1) <--- #
      # (1,0)
      # (2,0)
      # (0,2)
      #TODO

      # ---> Group Z(2) <--- #
      # (0,1)
      if (0,1) in self._available_orders:
          START = 0
          END   = START + n2
          self._B[(0,1)] = self._s2[START:END].reshape(n,n)

      # (1,1)
      if (1,1) in self._available_orders:
          START = n2
          END   = START + n2*N*3
          self._B[(1,1)] = self._s2[START:END].reshape(n,n,N,3)

      # (2,1)
      if (2,1) in self._available_orders:
          START = END
          END   = START + n2*N*6
          b     = self._s2[START:END].reshape(n,n,N,6)
          B     = numpy.zeros((n,n,N,3,3))
          xx = b[:,:,:,0]
          xy = b[:,:,:,1]
          xz = b[:,:,:,2]
          yy = b[:,:,:,3]
          yz = b[:,:,:,4]
          zz = b[:,:,:,5]
          B[:,:,:,0,0] = xx
          B[:,:,:,0,1] = xy
          B[:,:,:,0,2] = xz
          B[:,:,:,1,1] = yy
          B[:,:,:,1,2] = yz
          B[:,:,:,2,2] = zz
          B[:,:,:,1,0] = xy.copy()
          B[:,:,:,2,0] = xz.copy()
          B[:,:,:,2,1] = yz.copy()
          self._B[(2,1)] = B



  def _rotate(self, r, R, order):#TODO
      "Rotate DMS tensors"
      B = self.B(*order)
      #
      if   order == (1,0): 
           for i in range(self._N):
               Bi= numpy.einsum("abu,ac,bd,ux->cdx", B[:,:,i,:], R, R, r)
               B[:,:,i,:] = Bi
      #
      elif order == (2,0): 
           for i in range(self._N):
               Bi= numpy.einsum("abuw,ac,bd,xy->cdxy", B[:,:,i,:], R, R, r, r)
               B[:,:,i,:,:] = Bi
      #
      elif order == (0,1): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      elif order == (1,1): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      elif order == (0,2): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      elif order == (1,2): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      elif order == (2,1): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      elif order == (2,2): raise NotImplementedError(" DMS Error: Rotation (%i,%i)" % order)
      else: 
          raise ValueError(" DMS Error: Wrong order for rotation chosen")


class Basic_DMS(DMS):
  """
 Basic model of DMS that handles induction up to second-order and
 first-order pure Pauli effects.

 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,1)          Pauli          Asymmetric    (n,n)         Z(2)
"""
  def __init__(self, type='da'):#OK
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1)]


class LinearCT_DMS(DMS):
  """
 Basic model of DMS that handles induction up to second-order,
 first-order pure Pauli effects and first-order pure CT effects.

 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,1)          Pauli          Asymmetric    (n,n)         Z(2)
  4.            (1,1)          CT             Asymmetric    (n,n,N,3)     Z(2)
"""
  def __init__(self, type='da'):#OK
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1), (1,1)]


  # --> Implementation <-- #

class QuadraticCT_DMS(DMS):
  """
 Basic model of DMS that handles induction up to second-order,
 first-order pure Pauli effects and second-order pure CT effects.

 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,1)          Pauli          Asymmetric    (n,n)         Z(2)
  4.            (1,1)          CT             Asymmetric    (n,n,N,3)     Z(2)
  5.            (2,1)          CT             Asymmetric    (n,n,N,3,3)   Z(2)
"""
  def __init__(self, type='da'):#OK
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1), (1,1), (2,1)]


  # --> Implementation <-- #




class DMSFit(ABC):
  """
 Method to fit DMS tensors.
"""
  minimum_atom_atom_distance = 1.2 / psi4.constants.bohr2angstroms

  def __init__(self, mol, method, 
                     nsamples, dms_type, order_type, use_non_iterative_model):#OK
      ABC.__init__(self)

      self._mol       = mol
      self._method    = method
      self._natoms    = mol.natom()
      self._nsamples  = nsamples
      self._dms_type  = dms_type
      self._order_type= order_type
      psi4.set_options({"DMATPOL_TRAINING_MODE":"CHARGES",
                        "DMATPOL_FIELD_RANK"   : 2,
                        "DMATPOL_GRADIENT_RANK": 0,
                        "DMATPOL_NSAMPLES"     : nsamples,
                        "DMATPOL_NTEST_CHARGE" : 10,
                        "DMATPOL_TEST_CHARGE"  : 0.001})
      self._use_non_iterative_model = use_non_iterative_model
      self._wfn_0= None
      self._bfs_0= None
      self._nbf = None
      self._dms = None
      self._D0  = None

  @classmethod
  def create(cls, mol, fit_type="transl", dms_type="da", order_type='basic',
                  nsamples=100, method='scf', 
                  use_non_iterative_model=True):#OK
      if fit_type.lower().startswith("tran"): 
         return Translation_DMSFit(mol, method, nsamples, dms_type, order_type, use_non_iterative_model)
      else: 
         raise NotImplementedError("This type of DMS fitting is not implemented yet")

  def run(self):#TODO
      "Run the fitting procedure"
      self._compute_wfn()
     #self._compute_dms_external_field()
      self._compute_samples()
     #self._compute_group_1()
      self._compute_group_2()
      self._check()

  def B(self, m, n):
      "Get the susceptibilities"
      return self._dms.B(m,n)

  # --- protected --- #

  def _compute_wfn(self):#OK
      print(" * Computing Unperturbed Wavefunction")
      g, self._wfn_0 = PSI4_DRIVER(self._method, molecule=self._mol, return_wfn=True)
      self._nbf = self._wfn_0.basisset().nbf()
      self._dms = DMS.create(self._dms_type, self._order_type)
      self._dms.set_bfs(self._wfn_0.basisset())

      if   self._dms_type.lower() == 'da': M = self._wfn_0.Da().to_array(dense=True)
      elif self._dms_type.lower() == 'db': M = self._wfn_0.Db().to_array(dense=True)
      elif self._dms_type.lower() == 'fa': M = self._wfn_0.Fa().to_array(dense=True)
      elif self._dms_type.lower() == 'fb': M = self._wfn_0.Fb().to_array(dense=True)
      else:
          raise valueerror(" DMSFit Error: Incorrect type of dms tensors. Available: da, db, fa, fb.")

      self._dms.set_M(M)
      self._D0 = self._wfn_0.Da().to_array(dense=True)
      self._bfs_0 = self._wfn_0.basisset()

  def _compute_dms_external_field(self):#OK
      "DMS for external electric field only (without other molecules)"
      solver = oepdev.GenEffParFactory.build("POLARIZATION", self._wfn_0, self._wfn.options())
      self._dms_ind_0 = solver.compute()

  def _invert_hessian(self, h):#OK
     "Invert Hessian matrix"
     det= numpy.linalg.det(h)
     print(" * Hessian Determinant= %14.6E" % det)
     hi = numpy.linalg.inv(h)
     I = numpy.dot(hi, h).diagonal().sum()
     d = I - len(hi)
     if abs(d) > 0.0001: raise ValueError("Hessian is problemmatic! d=%f I= %f DIM=%f" % (d,I,h.shape[0]))
     return hi      

  def _compute_efield_due_to_fragment(self, mol_j, D_j, ints_j, ri):#OK
      "Compute electric field at ri due to mol_j with D_j and field integrals ints_j"
      fi= numpy.zeros(3)
      for j in range(mol_j.natom()):
          xj= mol_j.x(j)
          yj= mol_j.y(j)
          zj= mol_j.z(j)
          Zj= numpy.float(mol_j.Z(j))

          rij = numpy.array([ri[0]-xj, ri[1]-yj, ri[2]-zj])
          rij_norm = numpy.linalg.norm(rij)
          fi += Zj * rij / rij_norm**3

      fi[0] += (D_j @ ints_j[0]).trace() 
      fi[1] += (D_j @ ints_j[1]).trace() 
      fi[2] += (D_j @ ints_j[2]).trace() 
      return fi

  def _clash(self, geom_1, geom_2):#OK
      "Determine if two geometries clash or not"
      clash = False
      for a in geom_1:
          for b in geom_2:
              r_ab = numpy.sqrt(sum((a-b)**2))
              if r_ab < DMSFit.minimum_atom_atom_distance:
                 clash = True
                 break
      return clash

  # ---> Abstract methods <--- #

  @abstractmethod
  def _construct_aggregate(self): pass

  @abstractmethod
  def _compute_samples(self): pass

  @abstractmethod
  def _compute_group_1(self): pass

  @abstractmethod
  def _compute_group_2(self): pass

  @abstractmethod
  def _check(self): pass

class Translation_DMSFit(DMSFit):
  """
 Translation method to fit DMS tensors.
"""
  def __init__(self, mol, method, nsamples, dms_type, order_type, use_non_iterative_model,
                     start=1.5, srange=2.0):
      DMSFit.__init__(self, mol, method, nsamples, dms_type, order_type, use_non_iterative_model)

      # translation parameters
      self._i = 0
      self._start = start
      self._range = srange

      self._dimer_wfn = None
      self._dimer_mol = None
      self._mints_dimer = None

      self._g_ind = None
      self._h_ind = None
      self._g_ct  = None
      self._h_ct  = None

      self._dD_AA_set_ref = []
      self._dD_BB_set_ref = []
      self._dD_AB_set_ref = []
      self._S_AB_set = []
      self._T_set = []
      self._F_A_set = []
      self._F_B_set = []
      self._W_AB_set= []
      self._W_BA_set= []
      self._w_AB_set= []
      self._w_BA_set= []
     #self._Wf_AB_set= []
     #self._Wf_BA_set= []
      self._A_AB_set= []
      self._A_BA_set= []
      self._a_AB_set= []
      self._a_BA_set= []
      self._k_set   = []
      self._l_set   = []
     #self._Af_AB_set= []
     #self._Af_BA_set= []
     #self._Aff_AB_set= []
     #self._Aff_BA_set= []
      self._K_set   = []
      self._L_set   = []
      self._DIP_set = []      
      self._H_AB_set_ref = []


  # ---> Implementation <--- #

  def _construct_aggregate(self):#OK
      "Create next dimer by translating a molecule"
      psi4.core.clean()

      self._i += 1
      log = "\n0 1\n"
      for i in range(self._natoms):
          log += "%s" % self._mol.symbol(i)
          log += "%16.6f" % (self._mol.x(i) * psi4.constants.bohr2angstroms)
          log += "%16.6f" % (self._mol.y(i) * psi4.constants.bohr2angstroms)
          log += "%16.6f" % (self._mol.z(i) * psi4.constants.bohr2angstroms)
          log += "\n"
      log += "units angstrom\n"
      log += "symmetry c1\n"
      log += "no_reorient\n"
      log += "no_com\n"

      log += "--\n"
      log += "0 1\n"

      t = self._new_translation()

      for i in range(self._natoms):
          log += "%s" % self._mol.symbol(i)
          log += "%16.6f" % (self._mol.x(i) * psi4.constants.bohr2angstroms + t[0])
          log += "%16.6f" % (self._mol.y(i) * psi4.constants.bohr2angstroms + t[1])
          log += "%16.6f" % (self._mol.z(i) * psi4.constants.bohr2angstroms + t[2])
          log += "\n"
      log += "units angstrom\n"
      log += "symmetry c1\n"
      log += "no_reorient\n"
      log += "no_com\n"
      print(log)

      mol = psi4.geometry(log)
      mol.update_geometry()

      e_dimer, wfn_dimer = PSI4_DRIVER(self._method, molecule=mol, return_wfn=True)
      self._dimer_wfn = wfn_dimer
      self._dimer_mol = mol
      self._dimer_bfs = wfn_dimer.basisset()
      self._dimer_mints = psi4.core.MintsHelper(self._dimer_bfs)

      self._mol_A = self._mol
      self._mol_B = mol.extract_subsets(2)
      e_monomer, wfn_B = PSI4_DRIVER(self._method, molecule=self._mol_B, return_wfn=True)
      self._wfn_A = self._wfn_0
      self._wfn_B = wfn_B
      self._bfs_A = self._bfs_0
      self._bfs_B = wfn_B.basisset()

      e_int = e_dimer - 2.0 * e_monomer # MCBS
      print(" Interaction energy= %13.6f [a.u.]" % e_int)
      #
      return mol

  def _compute_samples(self):#OK
      "Compute electric fields, CT kernels and reference deformation density matrices"

      print(" * Computing WFN for Each Sample")
      for n in range(self._nsamples):
          print(" * - Sample %d" % (self._i+1))
          aggr= self._construct_aggregate()

          self._determine_perturbing_densities()

          S_AB = self._dimer_wfn.S().to_array(dense=True)[:self._nbf,self._nbf:]
          H_AB = self._dimer_wfn.H().to_array(dense=True)[:self._nbf,self._nbf:]
                                                                                                        
          dD     = self._dimer_wfn.Da().to_array(dense=True)
          dD[:self._nbf,:self._nbf]-= self._D0
          dD[self._nbf:,self._nbf:]-= self._D0
          self._save_dD(dD)

          F_A, F_B, F_A_mat, F_B_mat = self._compute_efield()


          dD_AA = dD[:self._nbf,:self._nbf].copy()
          dD_BB = dD[self._nbf:,self._nbf:].copy()
          dD_AB = dD[:self._nbf,self._nbf:].copy()

          W_AB = self._DA @ S_AB
          W_BA = self._DB @ S_AB.T

         #w_AB = self._DA @ S_AB    @ S_AB.T @ S_AB
         #w_BA = self._DB @ S_AB.T  @ S_AB   @ S_AB.T

         #Wf_AB= numpy.einsum("ixab,bc->acix",F_A_mat, W_AB) # i.e. F_A_mat[:,:,i,x] @ W_AB 
         #Wf_BA= numpy.einsum("ixab,bc->acix",F_B_mat, W_BA) # i.e. F_B_mat[:,:,i,x] @ W_BA 

          A_AB = W_AB.T @ W_AB
          A_BA = W_BA.T @ W_BA

         #a_AB = w_AB.T @ w_AB
         #a_BA = w_BA.T @ w_BA

         #Af_AB= numpy.einsum("acix,ad->cdix",Wf_AB, W_AB)
         #Af_BA= numpy.einsum("acix,ad->cdix",Wf_BA, W_BA)

         #Aff_AB= numpy.einsum("acix,adjy->cdixjy",Wf_AB, Wf_AB)
         #Aff_BA= numpy.einsum("acix,adjy->cdixjy",Wf_BA, Wf_BA)

          K    = W_AB.T @ dD_AB
          L    = W_BA.T @ dD_AB.T

         #k    = w_AB.T @ dD_AB
         #l    = w_BA.T @ dD_AB.T

          self._S_AB_set.append(S_AB.copy())
          self._H_AB_set_ref.append(H_AB.copy())
          self._W_AB_set.append(W_AB)
          self._W_BA_set.append(W_BA)
          self._A_AB_set.append(A_AB)
          self._A_BA_set.append(A_BA)

         #self._Wf_AB_set.append(Wf_AB)
         #self._Wf_BA_set.append(Wf_BA)

         #self._Af_AB_set.append(Af_AB)
         #self._Af_BA_set.append(Af_BA)

         #self._Aff_AB_set.append(Aff_AB)
         #self._Aff_BA_set.append(Aff_BA)


          self._K_set.append(K)
          self._L_set.append(L)

         #self._k_set.append(k)
         #self._l_set.append(l)
         #self._w_AB_set.append(w_AB)
         #self._w_BA_set.append(w_BA)
         #self._a_AB_set.append(a_AB)
         #self._a_BA_set.append(a_BA)


          self._dD_AA_set_ref.append(dD_AA)
          self._dD_BB_set_ref.append(dD_BB)
          self._dD_AB_set_ref.append(dD_AB)

          self._F_A_set.append(F_A)
          self._F_B_set.append(F_B)

          #DIP = mints.ao_dipole()
          #self._DIP_set.append(DIP)

      self._F_A_set = numpy.array(self._F_A_set)
      self._F_B_set = numpy.array(self._F_B_set)
      self._S_AB_set= numpy.array(self._S_AB_set)
      self._W_AB_set= numpy.array(self._W_AB_set)
      self._W_BA_set= numpy.array(self._W_BA_set)
      self._A_AB_set= numpy.array(self._A_AB_set)
      self._A_BA_set= numpy.array(self._A_BA_set)
     #self._Wf_AB_set= numpy.array(self._Wf_AB_set)
     #self._Wf_BA_set= numpy.array(self._Wf_BA_set)
     #self._Af_AB_set= numpy.array(self._Af_AB_set)
     #self._Af_BA_set= numpy.array(self._Af_BA_set)
     #self._Aff_AB_set= numpy.array(self._Aff_AB_set)
     #self._Aff_BA_set= numpy.array(self._Aff_BA_set)
      self._K_set   = numpy.array(self._K_set)
      self._L_set   = numpy.array(self._L_set)


     #self._w_AB_set= numpy.array(self._w_AB_set)
     #self._w_BA_set= numpy.array(self._w_BA_set)
     #self._a_AB_set= numpy.array(self._a_AB_set)
     #self._a_BA_set= numpy.array(self._a_BA_set)
     #self._k_set   = numpy.array(self._k_set)
     #self._l_set   = numpy.array(self._l_set)

      self._H_AB_set_ref = numpy.array(self._H_AB_set_ref)
      pass

  def _compute_group_1(self):#OK
      "Compute 1st group of parameters"

      # Only (1,0) and (2,0) susceptibilities: fitting for each target matrix element separately
      if (0,2) not in self._dms.available_orders():

          dms_ind1 = numpy.zeros((self._nbf, self._nbf, self._natoms, 3))
          dms_ind2 = numpy.zeros((self._nbf, self._nbf, self._natoms, 3, 3))

          for i in range(self._nbf):                              
              for j in range(self._nbf):
                  b_1, b_2 = self._compute_dms_induction_ij(i, j)
                  dms_ind1[i,j] = b1
                  dms_ind2[i,j] = b2

          par = [ dms_ind1.ravel(), dms_ind2.ravel() ]
          s = numpy.hstack( par )
 
      # higher-order susceptibilities included: need to evaluate full Hessian in S1 subspace
      else:
          raise NotImplementedError

      self._dms.set_s1(s)

  def _compute_group_2(self):#TODO
      "Compute 2nd group of parameters"

      OFFs= [0,]
      # (0,1)
      DIM = self._dms.n**2
      OFFs.append(DIM)
      # (1,1)
      if (1,1) in self._dms.available_orders():
          DIM +=  3*self._dms.N*self._dms.n**2
          OFFs.append(DIM)
      # (2,1)
      if (2,1) in self._dms.available_orders():
          DIM +=  6*self._dms.N*self._dms.n**2
          OFFs.append(DIM)
    
      # allocate 
      g = numpy.zeros(DIM)
      H = numpy.zeros((DIM, DIM))

      # gradient
      # (01) block
      print(" * Computing Gradient (0,1)...")
      START = OFFs[0]
      END   = OFFs[1]
      KL = self._K_set + self._L_set
      g[START:END] = KL.sum(axis=0).ravel()

      # (11) block
      if (1,1) in self._dms.available_orders():
          print(" * Computing Gradient (1,1)...")
          START = OFFs[1]
          END   = OFFs[2]
          #
          g[START:END] = numpy.einsum("nab,niu->abiu", self._K_set, self._F_A_set).ravel()
          g[START:END]+= numpy.einsum("nab,niu->abiu", self._L_set, self._F_B_set).ravel()

      # (21) block
      if (2,1) in self._dms.available_orders():
          print(" * Computing Gradient (2,1)...")
          START = OFFs[2]
          END   = OFFs[3]
          #
          gw = numpy.einsum("nab,niu,niw->abiuw",self._K_set, self._F_A_set, self._F_A_set)
          gw+= numpy.einsum("nab,niu,niw->abiuw",self._L_set, self._F_B_set, self._F_B_set)
          #
          u = numpy.zeros((self._dms.n,self._dms.n,self._dms.N,6))
          #
          u[:,:,:,0] = gw[:,:,:,0,0].copy()       # xx
          u[:,:,:,1] = gw[:,:,:,0,1].copy() * 2.0 # xy
          u[:,:,:,2] = gw[:,:,:,0,2].copy() * 2.0 # xz
          u[:,:,:,3] = gw[:,:,:,1,1].copy()       # yy
          u[:,:,:,4] = gw[:,:,:,1,2].copy() * 2.0 # yz
          u[:,:,:,5] = gw[:,:,:,2,2].copy()       # zz
          #
          g[START:END] = u.ravel().copy()
          del gw, u

          #g_xx = gw[:,:,:,0,0].ravel().copy()
          #g_xy = gw[:,:,:,0,1].ravel().copy() * 2.0 
          #g_xz = gw[:,:,:,0,2].ravel().copy() * 2.0
          #g_yy = gw[:,:,:,1,1].ravel().copy()
          #g_yz = gw[:,:,:,1,2].ravel().copy() * 2.0
          #g_zz = gw[:,:,:,2,2].ravel().copy()

          #g[START:END] = numpy.hstack([g_xx, g_xy, g_xz, g_yy, g_yz, g_zz])
          #g[START:END][:,:,:,0] = g_xx.ravel()  # xx
          #g[START:END][:,:,:,1] = g_xy.ravel()  # xy
          #g[START:END][:,:,:,2] = g_xz.ravel()  # xz
          #g[START:END][:,:,:,3] = g_yy.ravel()  # yy
          #g[START:END][:,:,:,4] = g_yz.ravel()  # yz
          #g[START:END][:,:,:,5] = g_zz.ravel()  # zz

      g *= -2.0

      # Hessian
      I = numpy.identity(self._nbf)
      A = self._A_AB_set + self._A_BA_set
     #Af= self._Af_AB_set + self._Af_BA_set
     #Aff= self._Aff_AB_set + self._Aff_BA_set

      # (01) susceptibility
      print(" * Computing Hessian (0,1) (0,1)...")
      START = OFFs[0]
      END   = OFFs[1]
      L = END - START
      H_01_01 = numpy.einsum("bd,ac->abcd",I,A.sum(axis=0))
      H_01_01+= numpy.einsum("nda,nbc->abcd",self._W_AB_set,self._W_BA_set)
      H_01_01+= numpy.einsum("nda,nbc->abcd",self._W_BA_set,self._W_AB_set)
      H[START:END,START:END] = H_01_01.reshape(L,L).copy() ; del H_01_01

      # (11) susceptibility
      if (1,1) in self._dms.available_orders():
          PREV  = OFFs[0]
          START = OFFs[1]
          END   = OFFs[2]
          L = END - START
          M = START - PREV
          #
          print(" * Computing Hessian (1,1) (1,1)...")
          H_11_11 = numpy.einsum("bd,nac,niu,njw->abiucdjw",I,self._A_AB_set,self._F_A_set,self._F_A_set)               
          H_11_11+= numpy.einsum("bd,nac,niu,njw->abiucdjw",I,self._A_BA_set,self._F_B_set,self._F_B_set)
          H_11_11+= numpy.einsum("nda,nbc,niu,njw->abiucdjw",self._W_AB_set,self._W_BA_set,self._F_A_set,self._F_B_set)
          H_11_11+= numpy.einsum("nda,nbc,niu,njw->abiucdjw",self._W_BA_set,self._W_AB_set,self._F_B_set,self._F_A_set)
          #
          print(" * Computing Hessian (0,1) (1,1)...")
          H_01_11 = numpy.einsum("bd,nac,njw->abcdjw",I,self._A_AB_set,self._F_A_set)
          H_01_11+= numpy.einsum("bd,nac,njw->abcdjw",I,self._A_BA_set,self._F_B_set)
          H_01_11+= numpy.einsum("nda,nbc,njw->abcdjw",self._W_AB_set,self._W_BA_set,self._F_B_set)
          H_01_11+= numpy.einsum("nda,nbc,njw->abcdjw",self._W_BA_set,self._W_AB_set,self._F_A_set)
          #                                                                                                              
          H[START:  END,START:  END] = H_11_11.reshape(L,L).copy() ; del H_11_11
          H[PREV :START,START:  END] = H_01_11.reshape(M,L).copy()
          H[START:  END, PREV:START] = H_01_11.reshape(M,L).T.copy() ; del H_01_11

      # (21) susceptibility
      if (2,1) in self._dms.available_orders():
          DPREV = OFFs[0]  
          PREV  = OFFs[1]
          START = OFFs[2]
          END   = OFFs[3]
          L = END - START
          M = START - PREV
          N = PREV - DPREV
          #
          print(" * Computing Hessian (2,1) (2,1)...")
          H_21_21 = numpy.einsum("bd,nac,niu,nix,njw,njy->abiuxcdjwy", I,self._A_AB_set,self._F_A_set,self._F_A_set,
                                                                                      self._F_A_set,self._F_A_set)

          H_21_21+= numpy.einsum("bd,nac,niu,nix,njw,njy->abiuxcdjwy", I,self._A_BA_set,self._F_B_set,self._F_B_set,
                                                                                      self._F_B_set,self._F_B_set)

          H_21_21+= numpy.einsum("nda,nbc,niu,nix,njw,njy->abiuxcdjwy", self._W_AB_set,self._W_BA_set,
                                                                                      self._F_A_set,self._F_A_set,
                                                                                      self._F_B_set,self._F_B_set)

          H_21_21+= numpy.einsum("nda,nbc,niu,nix,njw,njy->abiuxcdjwy", self._W_BA_set,self._W_AB_set,
                                                                                      self._F_B_set,self._F_B_set,
                                                                                      self._F_A_set,self._F_A_set)

          U=numpy.zeros((self._dms.n,self._dms.n,self._dms.N,6,self._dms.n,self._dms.n,self._dms.N,6))
          #
          U[:,:,:,0,:,:,:,0] = H_21_21[:,:,:,0,0,:,:,:,0,0].copy()       # xx,xx
          U[:,:,:,0,:,:,:,1] = H_21_21[:,:,:,0,0,:,:,:,0,1].copy() * 2.0 # xx,xy
          U[:,:,:,0,:,:,:,2] = H_21_21[:,:,:,0,0,:,:,:,0,2].copy() * 2.0 # xx,xz
          U[:,:,:,0,:,:,:,3] = H_21_21[:,:,:,0,0,:,:,:,1,1].copy()       # xx,yy
          U[:,:,:,0,:,:,:,4] = H_21_21[:,:,:,0,0,:,:,:,1,2].copy() * 2.0 # xx,yz
          U[:,:,:,0,:,:,:,5] = H_21_21[:,:,:,0,0,:,:,:,2,2].copy()       # xx,zz
          #
          U[:,:,:,1,:,:,:,0] = H_21_21[:,:,:,0,1,:,:,:,0,0].copy() * 2.0 # xy,xx
          U[:,:,:,1,:,:,:,1] = H_21_21[:,:,:,0,1,:,:,:,0,1].copy() * 4.0 # xy,xy
          U[:,:,:,1,:,:,:,2] = H_21_21[:,:,:,0,1,:,:,:,0,2].copy() * 4.0 # xy,xz
          U[:,:,:,1,:,:,:,3] = H_21_21[:,:,:,0,1,:,:,:,1,1].copy() * 2.0 # xy,yy
          U[:,:,:,1,:,:,:,4] = H_21_21[:,:,:,0,1,:,:,:,1,2].copy() * 4.0 # xy,yz
          U[:,:,:,1,:,:,:,5] = H_21_21[:,:,:,0,1,:,:,:,2,2].copy() * 2.0 # xy,zz
          #
          U[:,:,:,2,:,:,:,0] = H_21_21[:,:,:,0,2,:,:,:,0,0].copy() * 2.0 # xz,xx
          U[:,:,:,2,:,:,:,1] = H_21_21[:,:,:,0,2,:,:,:,0,1].copy() * 4.0 # xz,xy
          U[:,:,:,2,:,:,:,2] = H_21_21[:,:,:,0,2,:,:,:,0,2].copy() * 4.0 # xz,xz
          U[:,:,:,2,:,:,:,3] = H_21_21[:,:,:,0,2,:,:,:,1,1].copy() * 2.0 # xz,yy
          U[:,:,:,2,:,:,:,4] = H_21_21[:,:,:,0,2,:,:,:,1,2].copy() * 4.0 # xz,yz
          U[:,:,:,2,:,:,:,5] = H_21_21[:,:,:,0,2,:,:,:,2,2].copy() * 2.0 # xz,zz
          #
          U[:,:,:,3,:,:,:,0] = H_21_21[:,:,:,1,1,:,:,:,0,0].copy()       # yy,xx
          U[:,:,:,3,:,:,:,1] = H_21_21[:,:,:,1,1,:,:,:,0,1].copy() * 2.0 # yy,xy
          U[:,:,:,3,:,:,:,2] = H_21_21[:,:,:,1,1,:,:,:,0,2].copy() * 2.0 # yy,xz
          U[:,:,:,3,:,:,:,3] = H_21_21[:,:,:,1,1,:,:,:,1,1].copy()       # yy,yy
          U[:,:,:,3,:,:,:,4] = H_21_21[:,:,:,1,1,:,:,:,1,2].copy() * 2.0 # yy,yz
          U[:,:,:,3,:,:,:,5] = H_21_21[:,:,:,1,1,:,:,:,2,2].copy()       # yy,zz
          #
          U[:,:,:,4,:,:,:,0] = H_21_21[:,:,:,1,2,:,:,:,0,0].copy() * 2.0 # yz,xx
          U[:,:,:,4,:,:,:,1] = H_21_21[:,:,:,1,2,:,:,:,0,1].copy() * 4.0 # yz,xy
          U[:,:,:,4,:,:,:,2] = H_21_21[:,:,:,1,2,:,:,:,0,2].copy() * 4.0 # yz,xz
          U[:,:,:,4,:,:,:,3] = H_21_21[:,:,:,1,2,:,:,:,1,1].copy() * 2.0 # yz,yy
          U[:,:,:,4,:,:,:,4] = H_21_21[:,:,:,1,2,:,:,:,1,2].copy() * 4.0 # yz,yz
          U[:,:,:,4,:,:,:,5] = H_21_21[:,:,:,1,2,:,:,:,2,2].copy() * 2.0 # yz,zz
          #
          U[:,:,:,5,:,:,:,0] = H_21_21[:,:,:,2,2,:,:,:,0,0].copy()       # zz,xx
          U[:,:,:,5,:,:,:,1] = H_21_21[:,:,:,2,2,:,:,:,0,1].copy() * 2.0 # zz,xy
          U[:,:,:,5,:,:,:,2] = H_21_21[:,:,:,2,2,:,:,:,0,2].copy() * 2.0 # zz,xz
          U[:,:,:,5,:,:,:,3] = H_21_21[:,:,:,2,2,:,:,:,1,1].copy()       # zz,yy
          U[:,:,:,5,:,:,:,4] = H_21_21[:,:,:,2,2,:,:,:,1,2].copy() * 2.0 # zz,yz
          U[:,:,:,5,:,:,:,5] = H_21_21[:,:,:,2,2,:,:,:,2,2].copy()       # zz,zz
          # 
          H[START:END,START:END] = U.reshape(L,L).copy()
          del H_21_21, U

          # (01,21)
          print(" * Computing Hessian (0,1) (2,1)...")
          H_01_21 = numpy.einsum("bd,nac,njw,njx->abcdjwx", I,self._A_AB_set,self._F_A_set,self._F_A_set)
          H_01_21+= numpy.einsum("bd,nac,njw,njx->abcdjwx", I,self._A_BA_set,self._F_B_set,self._F_B_set)
          H_01_21+= numpy.einsum("nda,nbc,njw,njx->abcdjwx", self._W_AB_set,self._W_BA_set,self._F_B_set,self._F_B_set)
          H_01_21+= numpy.einsum("nda,nbc,njw,njx->abcdjwx", self._W_BA_set,self._W_AB_set,self._F_A_set,self._F_A_set)
          #
          U = numpy.zeros((self._dms.n,self._dms.n,self._dms.n,self._dms.n,self._dms.N,6))
          #
          U[:,:,:,:,:,0] = H_01_21[:,:,:,:,:,0,0].copy()       # xx
          U[:,:,:,:,:,1] = H_01_21[:,:,:,:,:,0,1].copy() * 2.0 # xy
          U[:,:,:,:,:,2] = H_01_21[:,:,:,:,:,0,2].copy() * 2.0 # xz
          U[:,:,:,:,:,3] = H_01_21[:,:,:,:,:,1,1].copy()       # yy
          U[:,:,:,:,:,4] = H_01_21[:,:,:,:,:,1,2].copy() * 2.0 # yz
          U[:,:,:,:,:,5] = H_01_21[:,:,:,:,:,2,2].copy()       # zz
          #
          H[DPREV:PREV,START:END] = U.reshape(N,L).copy()
          del H_01_21, U
          H[START:END,DPREV:PREV] = H[DPREV:PREV,START:END].T.copy()

          # (11,21)
          print(" * Computing Hessian (1,1) (2,1)...")
          H_11_21 = numpy.einsum("bd,nac,nix,nju,njw->abixcdjuw", I,self._A_AB_set,self._F_A_set,self._F_A_set,self._F_A_set)
          H_11_21+= numpy.einsum("bd,nac,nix,nju,njw->abixcdjuw", I,self._A_BA_set,self._F_B_set,self._F_B_set,self._F_B_set)
          H_11_21+= numpy.einsum("nda,nbc,nix,nju,njw->abixcdjuw",self._W_AB_set,self._W_BA_set,
                                                                  self._F_A_set,self._F_B_set,self._F_B_set)
          H_11_21+= numpy.einsum("nda,nbc,nix,nju,njw->abixcdjuw",self._W_BA_set,self._W_AB_set,
                                                                  self._F_B_set,self._F_A_set,self._F_A_set)
          #
          U = numpy.zeros((self._dms.n,self._dms.n,self._dms.N,3,self._dms.n,self._dms.n,self._dms.N,6))
          #
          U[:,:,:,:,:,:,:,0] = H_11_21[:,:,:,:,:,:,:,0,0].copy()       # xx
          U[:,:,:,:,:,:,:,1] = H_11_21[:,:,:,:,:,:,:,0,1].copy() * 2.0 # xy
          U[:,:,:,:,:,:,:,2] = H_11_21[:,:,:,:,:,:,:,0,2].copy() * 2.0 # xz
          U[:,:,:,:,:,:,:,3] = H_11_21[:,:,:,:,:,:,:,1,1].copy()       # yy
          U[:,:,:,:,:,:,:,4] = H_11_21[:,:,:,:,:,:,:,1,2].copy() * 2.0 # yz
          U[:,:,:,:,:,:,:,5] = H_11_21[:,:,:,:,:,:,:,2,2].copy()       # zz
          # 
          H[PREV:START,START:END] = U.reshape(M,L).copy()
          del H_11_21, U
          H[START:END,PREV:START] = H[PREV:START,START:END].T.copy()


      H *= 2.0

      # fit
      print(" * Computing Hessian Inverse...")
      Hi = self._invert_hessian(H)

      print(" * Computing the DMS tensors...")
      s = - g @ Hi

      # extract susceptibilities
      print(" * Setting the DMS tensors...")
      self._dms.set_s2(s)
      for order in [(0,1),(1,1),(2,1)]:
          if order in self._dms.available_orders():
             self._dms._generate_B_from_s(order)

  def _compute_dms_induction_ij(self, i, j):#TODO
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

  def _check(self):
      "Check the quality of fitting on training set"

      B_01 = self.B(0,1)
      if (1,1) in self._dms.available_orders(): 
          B_11 = self.B(1,1)
      if (2,1) in self._dms.available_orders(): 
          B_21 = self.B(2,1)

      print(" DMSFit: Check - training set")
      for s in range(self._nsamples):
          
          W_AB = self._W_AB_set[s]
          W_BA = self._W_BA_set[s]
         #w_AB = self._w_AB_set[s]
         #w_BA = self._w_BA_set[s]
         #Wf_AB= self._Wf_AB_set[s]
         #Wf_BA= self._Wf_BA_set[s]
          F_A  = self._F_A_set[s]
          F_B  = self._F_B_set[s]

          H_AB_ref = self._H_AB_set_ref[s]

          dD_AB_ref = self._dD_AB_set_ref[s]
          dD_AB_com = W_AB @ B_01 + B_01.T @ W_BA.T
         #dD_AB_com+= numpy.einsum("acix,cbix->ab",Wf_AB, B_11)
         #dD_AB_com+= numpy.einsum("acix,cbix->ab",Wf_BA, B_11).T
          if (1,1) in self._dms.available_orders():
              dD_AB_com+= numpy.einsum("ac,ix,cbix->ab",W_AB, F_A, B_11)
              dD_AB_com+= numpy.einsum("ac,ix,cbix->ab",W_BA, F_B, B_11).T
          if (2,1) in self._dms.available_orders():
              dD_AB_com+= numpy.einsum("ac,ix,iy,cbixy->ab",W_AB, F_A, F_A, B_21)
              dD_AB_com+= numpy.einsum("ac,ix,iy,cbixy->ab",W_BA, F_B, F_B, B_21).T

          rms = self._rms(dD_AB_ref, dD_AB_com)

          r = (dD_AB_ref @ H_AB_ref).trace()
          c = (dD_AB_com @ H_AB_ref).trace()

          print(" Sample=%03d  RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f" \
                   % (s+1, rms, dD_AB_ref[0,0], dD_AB_com[0,0], r, c))

          self._i = s
          self._save_dD(dD_AB_ref, prefix='ref')
          self._save_dD(dD_AB_com, prefix='com')





  # -----> Utility methods <----- #

  def _compute_efield(self):#OK
      "Compute instantaneous electric field on fragment due to other fragment"

      # field due to B evaluated on A atoms
      F_B_set = []
      F_B_mat_set = []
      for i in range(self._mol_A.natom()):
          xi= self._mol_A.x(i)
          yi= self._mol_A.y(i)
          zi= self._mol_A.z(i)

          ri = numpy.array([xi,yi,zi])
          ints = [x.to_array(dense=True)[self._nbf:,self._nbf:] for x in self._dimer_mints.electric_field(origin=ri)]

          f = self._compute_efield_due_to_fragment(self._mol_B, self._DB, ints, ri)
          F_B_set.append(f)
          F_B_mat_set.append(numpy.array(ints))

      F_B_set = numpy.array(F_B_set)
      F_B_mat_set = numpy.array(F_B_mat_set)

      # field due to A evaluated on B atoms
      F_A_set = []
      F_A_mat_set = []
      for i in range(self._mol_B.natom()):
          xi= self._mol_B.x(i)
          yi= self._mol_B.y(i)
          zi= self._mol_B.z(i)

          ri = numpy.array([xi,yi,zi])
          ints = [x.to_array(dense=True)[:self._nbf,:self._nbf] for x in self._dimer_mints.electric_field(origin=ri)]

          f = self._compute_efield_due_to_fragment(self._mol_A, self._DA, ints, ri)
          F_A_set.append(f)
          F_A_mat_set.append(numpy.array(ints))

      F_A_set = numpy.array(F_A_set)
      F_A_mat_set = numpy.array(F_A_mat_set)

      return F_A_set, F_B_set, F_A_mat_set, F_B_mat_set
   


  def _rms(self, a, b): return numpy.sqrt(((a-b)**2).sum()/a.size)

  def _draw_translation(self):#OK
      theta = numpy.arccos(numpy.random.random())
      phi   = 2.0 * numpy.pi * numpy.random.random()
      r     = self._start + self._range * numpy.sqrt(numpy.random.random())
      x     = r * numpy.sin(theta) * numpy.cos(phi)
      y     = r * numpy.sin(theta) * numpy.sin(phi)
      z     = r * numpy.cos(theta)                 

      t     = numpy.array([x,y,z])
      return t

  def _new_translation(self):#OK
      "Select new random translation"
      done = False
      geom = self._mol.geometry().to_array(dense=True) 
      while not done:
         t = self._draw_translation()
         if not self._clash(geom, geom+t):
            done = True 
      return t

  def _save_dD(self, D, prefix='dd'):#OK
      name= prefix + '_%03d.dat' 
      out = open(name % self._i, 'w')
      log = ''
      n = len(D)
      for row in D:
          log+= n*"%13.5E" % tuple(row) + "\n"
      out.write(log)
      out.close() 

  def _determine_perturbing_densities(self):
      if self._use_non_iterative_model:
         DA = self._D0.copy()
         DB = DA.copy()
      else:
         D = self._dimer_wfn.Da().to_array(dense=True)
         DA= D[:self._nbf,:self._nbf]
         DB= D[self._nbf:,self._nbf:]

      self._DA = DA
      self._DB = DB




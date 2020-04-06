#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 One-Particle Density Matrix Susceptibility Module.
 Bartosz Błasiak, Gundelfingen, May 2019

 Implements the charge-transfer DMS tensors.
"""

import numpy, math
import psi4, oepdev
import scipy.optimize
import scipy.spatial.transform
from abc import ABC, abstractmethod
from ..math.matrix import move_atom_rotate_molecule, rotate_ao_matrix, matrix_power
from ..math import composite
from ..core.utilities import psi_molecule_from_file
from ..density.opdm import Density
from ..solvshift.mdinput import MDInput
from ..solvshift.fragment import Fragment


__all__ = ["DMSFit", "DMS", "Computer", "Aggregate"]

PSI4_DRIVER = psi4.energy
numpy.random.seed(0)

class Aggregate:
  def __init__(self, psi4_molecule, qm_frags=None):
      self.all = psi4_molecule
      self.qm_frags = qm_frags
      self.nfrags = psi4_molecule.nfragments()
      self.qm = None
      if qm_frags is not None: self.qm = psi4_molecule.extract_subsets(qm)
      self.bath = [] #if self.nfrags == 1 else [psi4_molecule.extract_subsets(2+i) for i in range(self.nfrags-1)]
      self._mrec = 1./(numpy.array([self.all.mass(i) for i in range(self.all.natom())])) * psi4.constants.au2amu
  def update(self, xyz):
      self.all.set_geometry(xyz)
      if self.qm_frags is not None: self.qm = self.all.extract_subsets(self.qm_frags)
      self.bath = [] 
  def tofile(self, out, center_mode=None, format='xyz', misc=None):
      geom = self.all.geometry().clone()
      geom.scale(psi4.constants.bohr2angstroms)
      if   center_mode is None         : com = [0.0,0.0,0.0]
      elif center_mode.lower() == 'qm' : com = self.qm.center_of_mass()
      elif center_mode.lower() == 'all': com = self.all.center_of_mass()
      elif isinstance(center_mode, int): com = self.bath[center_mode].center_of_mass()
      else: raise ValueError("Centering mode - %s - is not supported" % center_mode)
      if format.lower() =='xyz':
         out.write("%d\n%s\n" % (self.all.natom(),misc))                                                                          
         for i in range(self.all.natom()):                                                              
             sym = self.all.label(i)
             out.write("%s %10.6f %10.6f %10.6f\n"%(sym,geom.get(i,0)-com[0],geom.get(i,1)-com[1],geom.get(i,2)-com[2]))
      else:
         raise ValueError("Incorrect format for structure file provided.")

class Computer(ABC): 
   """
 ---------------------------------------------------------------------------------------------------
 Method to solve DMS-SCF problems.

 Usage:
  computer = Computer(aggregate, fragments)
  energy = computer.run(field=[0,0,0])
  mu     = computer.dipole_moment()
  mu     = computer.ff_dipole_moment()
  alpha  = computer.ff_polarizability(energy=False)
  beta   = computer.ff_hyperpolarizability(energy=False)
  gamma  = computer.ff_hyper2polarizability(energy=Fanse)
  camm   = computer.camm()
  charges= computer.atomic_charges(kappa=0.0)

  computer.update(psi_molecule)
 ---------------------------------------------------------------------------------------------------
                                                                           Gundelfingen, 31 Mar 2020
"""
   # Global defaults
   verbose                             = True
   e_convergence                       = 1.0e-6
   dmsscf_maximum_number_of_iterations = 1000
   raise_error_when_unconverged        = True

   def __init__(self, aggregate, fragments): # use solvshift.MDInput class in the future
       ABC.__init__(self)
       #
       self._aggregate = aggregate
       if isinstance(fragments, list):
          assert(len(fragments) == aggregate.nfrags)
          self._fragments = fragments
       elif isinstance(fragments, str):
          args, idx = MDInput(fragments).get()
          ind, nmol, bsm, supl, reord = args
          assert len(idx)==aggregate.all.natom()
          self._fragments = [bsm[x].copy() for x in ind]
       else: raise ValueError("Wrong type of data passed into fragments. Only list or str are possible.")
       #
       self._number_of_fragments = aggregate.nfrags
       #
       self._bfs = None
       self._nbf = None
       self._mints = None 
       self._jk = None
       self._total_energy = None
       self._S = None
       self._H = None
       self._electric_field_ao_integrals = None
       self._Ca= None
       self._Da= None; self._Db = None
       self._Da_0 = None; self._Db_0 = None
       self._Fa= None; self._Fb = None
       self._ao_offsets_by_fragment = None
       self._natom_offsets_by_fragment = None
       self._n = None # --> deprecate (no need of NO occupancies anymore)
       #
       self._build_bfs()
       self._allocate_memory()

   def run(self, field=None, verbose=None, converge_by='energy'):
       "Perform single-point energy or density calculation"
       verbose_stash = Computer.verbose
       if verbose is not None: Computer.verbose = verbose
       E = self._run_dmsscf(field, converge_by)
       Computer.verbose = verbose_stash
       return E

   def update(self, psi_molecule):
       "Updates for the next calculation"
       psi4.core.clean()
       del self._mints
       self._jk.finalize()
       self._aggregate.update(psi_molecule.geometry())
       self._build_bfs()

   def Da(self): 
       "Return OPDM alpha matrix"
       return self._Da.copy()

   def Da_0(self):
       "Return zeroth-order OPDM alpha matrix"
       return self._Da_0.copy()

   def Fa(self): 
       "Return Fock alpha matrix"
       return self._Fa.copy()

   def total_energy(self): 
       "Return total energy of the system"
       return self._total_energy

   def E(self): return self._total_energy

   # --> Protected Interface <-- #

   @abstractmethod
   def _orthogonalize_density_matrix(self): pass

   @abstractmethod
   def _run_dmsscf(self, field): pass

   def _allocate_memory(self):
       n = self._nbf
       self._Da = numpy.zeros((n,n))
       self._Fa = numpy.zeros((n,n))

       t = numpy.cumsum([self._fragments[x].get_nbf() for x in range(self._number_of_fragments)])
       self._ao_offsets_by_fragment = numpy.hstack(([0],t[:-1]))
       t = numpy.cumsum([self._fragments[x].get_natoms() for x in range(self._number_of_fragments)])
       self._natom_offsets_by_fragment = numpy.hstack(([0],t[:-1]))

   def _build_bfs(self):#OK
       self._bfs       = psi4.core.BasisSet.build(self._aggregate.all, "BASIS", psi4.core.get_global_option("BASIS"), puream=False)
       self._mints     = psi4.core.MintsHelper(self._bfs)
       self._jk        = psi4.core.JK.build(self._bfs, jk_type="direct")
       self._jk.set_memory(int(5e8))
       self._jk.initialize()

       self._nbf = self._bfs.nbf()

class _DMS_SCF_Procedure(Computer):
   "The DMS-SCF Procedure"
   def __init__(self, aggregate, fragments):
       Computer.__init__(self, aggregate, fragments)

   # --> Implementation <-- #

   def _run_dmsscf(self, field=None, converge_by='energy'):
       "Perform single-point energy or density calculation"
       self.__superimpose_fragments()
       self.__compute_one_electron_integrals(field)
       self.__initialize_opdm()
       self.__compute_total_energy(field, converge_by)
       return self._total_energy

   # --> Private Interface <-- #

   def __initialize_opdm(self):
       "Initializes the OPDM of the entire system. As for now, it takes the last OPDM as guess, otherwise gas-phase guess."
       if 1: #self._Da is not None:
          self._Da = numpy.zeros((self._nbf, self._nbf))    
          for I in range(self._number_of_fragments):
              off = self._ao_offsets_by_fragment[I]
              nbf = self._fragments[I].get_nbf()
              Da  = self._fragments[I].get()["opdm"]
              self._Da[off:off+nbf,off:off+nbf] = Da.copy()
          self._Da_0 = self._Da.copy()

   def __superimpose_fragments(self):
       "Perform superposition of all fragments"
       for I in range(self._number_of_fragments):
           mol = self._aggregate.all.extract_subsets(I+1)
           xyz = mol.geometry().to_array(dense=True)
           rms = self._fragments[I].superimpose(xyz)
 
   def __compute_one_electron_integrals(self, field):
       "Compute necessary one-electron integrals"
       T = self._mints.ao_kinetic().to_array(dense=True)
       V = self._mints.ao_potential().to_array(dense=True)         
       H = T + V
       if field is not None:
          DIP = self._mints.ao_dipole()
          for x in range(3): H -= DIP[x].to_array(dense=True) * field[x]
       self._H = H

       self._S = self._mints.ao_overlap().to_array(dense=True)

       self._electric_field_ao_integrals = []
       for i in range(self._aggregate.all.natom()):
           xi = self._aggregate.all.x(i)
           yi = self._aggregate.all.y(i)
           zi = self._aggregate.all.z(i)
           ri = numpy.array([xi,yi,zi])
           ints = [x.to_array(dense=True) for x in self._mints.electric_field(origin=ri)] # ints at ri
           self._electric_field_ao_integrals.append(ints)
       self._electric_field_ao_integrals = numpy.array(self._electric_field_ao_integrals)

       self._Ca = numpy.zeros((self._nbf, self._nbf)) # --> can probably be deprecated
       self._n  = None #numpy.zeros(self._nbf)
       for I in range(self._number_of_fragments):
           par = self._fragments[I].get()
           nbf = self._fragments[I].get_nbf()
           off = self._ao_offsets_by_fragment[I]
           Ca = numpy.hstack([par['caocc'], par['cavir']])
           self._Ca[off:off+nbf,off:off+nbf] = Ca.copy()

   def __compute_electric_field_due_to_fragment(self, J_list, D_J, ints_J, ri):
       "Compute electric field at ri due to mol of dimer with D_J and field integrals ints_J"
       fi= numpy.zeros(3)
       for J in J_list:
           off_natom_J = self._natom_offsets_by_fragment[J]
           for j in range(self._fragments[J].get_natoms()):          
               xj=             self._aggregate.all.x(off_natom_J+j) 
               yj=             self._aggregate.all.y(off_natom_J+j)
               zj=             self._aggregate.all.z(off_natom_J+j)
               Zj= numpy.float(self._aggregate.all.Z(off_natom_J+j))
                                                                     
               rij = numpy.array([ri[0]-xj, ri[1]-yj, ri[2]-zj])
               rij_norm = numpy.linalg.norm(rij)
               fi += Zj * rij / rij_norm**3

       fi[0] += 2.0 * (D_J @ ints_J[0]).trace() 
       fi[1] += 2.0 * (D_J @ ints_J[1]).trace() 
       fi[2] += 2.0 * (D_J @ ints_J[2]).trace() 
       return fi

   def __rms(self, a, b): return ((a-b)**2).sum()

   def __compute_total_energy(self, field, converge_by):
       "DMS SCF Procedure for N-Body System"
       # Initialize
       Iter = 0
       current_dmsscf_energy = 1.0e+10
       error = 1.0e+10
       success = True

       # Start
       if Computer.verbose:
          print(" @DMS-SCF: Starting iterative process for aggregate with %d fragments" % self._number_of_fragments)

       while error > Computer.e_convergence:

           # Check if maximum number of iterations exceeded
           if Iter == Computer.dmsscf_maximum_number_of_iterations: 
              success = False 
              break

           # Run microiterations for each pair of fragments: Older Code
           A = 2 # the best many-body model so far
           if A==0:
             for i in range(self._number_of_fragments):
               energy = self.__dmsscf_iteration_fragment(i, field)
               if Computer.verbose and psi4.core.get_global_option("PRINT")>1:
                  print(" @DMS-SCF: MicroIter I=%4d E=%14.6f" % (i+1, energy))

           elif A==1:
             for i in range(self._number_of_fragments):
               energy = self.__dmsscf_iteration_fragment(i, field)
               for j in range(i):
                   energy = self.__dmsscf_for_pair_of_fragments(i,j, field)
                   if Computer.verbose and psi4.core.get_global_option("PRINT")>1:
                      print(" @DMS-SCF: MicroIter I=%4d J=%4d E=%14.6f" % (i+1, j+1, energy))

           # Run globally including all possible many-body interactions
           else:
              energy = self.__dmsscf_iteration(field, converge_by)

           # Include Pauli deformation (experimental: use OtherComputer class) --> deprecate
           self._orthogonalize_density_matrix()

           # Compute total energy and error
           if converge_by == 'energy':
              error = abs(energy - current_dmsscf_energy)
              current_dmsscf_energy = energy
           elif converge_by == 'dipole':
              error = self.__rms(energy, current_dmsscf_energy)
              current_dmsscf_energy = energy
              energy = numpy.linalg.norm(energy)
           else:
              error = self.__rms(energy, current_dmsscf_energy)
              current_dmsscf_energy = energy
              energy = energy.sum()

           if Computer.verbose:
              print(" @DMS-SCF: Iter.%4d E=%14.6f" % (Iter, energy))

           Iter += 1

       # Finalize
       self._total_energy = energy

       if success:
          if Computer.verbose:
             print(" @DMS-SCF: Converged. Final Energy = %14.8f" % (energy))
       else:
          if Computer.verbose:
             print(" @DMS-SCF: Did not converged after %i iterations" % Computer.dmsscf_maximum_number_of_iterations)
          if Computer.raise_error_when_unconverged: 
             raise ValueError(" DMS-SCF did not converge!")

   def __dmsscf_iteration(self, field, converge_by):
       "1 Global Iteration of DMS SCF Procedure. So far the best version of the DMS SCF."

       do_23_body = False
       do_45_body = False
       do_67_body = False

       # Extract OPDM and DMS tensors
       nat = [x.get_natoms() for x in self._fragments]
       par = [x.get() for x in self._fragments]
       D0  = [x['opdm'  ] for x in par]
       B_10= [x['dms_10'] for x in par] 
       B_20= [x['dms_20'] for x in par] 
       try: 
            B_01= [x['dms_01'] for x in par] 
            do_23_body = True
       except KeyError: 
            B_01= [None for x in par] 
            do_23_body = False
       try: B_02= [x['dms_02'] for x in par] 
       except KeyError: B_02= [None for x in par] 

       try: 
            B_03= [x['dms_03'] for x in par] 
            do_45_body = True
       except KeyError: 
            B_03= [None for x in par] 
            do_45_body = False
       try: B_04= [x['dms_04'] for x in par] 
       except KeyError: B_04= [None for x in par] 

       try: 
            B_05= [x['dms_05'] for x in par] 
            do_67_body = True
       except KeyError: 
            B_05= [None for x in par] 
            do_67_body = False
       try: B_06= [x['dms_06'] for x in par] 
       except KeyError: B_06= [None for x in par] 

       # Compute perturbations: Electric field
       F = []
       for I in range(self._number_of_fragments):
           nbf_I = self._fragments[I].get_nbf()
           natom_I = nat[I]
           off_ao_I = self._ao_offsets_by_fragment[I]
           off_natom_I = self._natom_offsets_by_fragment[I]

           D_JJ = self._Da.copy()
           D_JJ[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].fill(0.0)

           JJ_list = []
           for J in range(self._number_of_fragments):
               if I!=J: 
                  JJ_list.append(J)
                  # local field approximation:
                  off_ao_J = self._ao_offsets_by_fragment[J]
                  nbf_J = self._fragments[J].get_nbf()
                  D_JJ[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J].fill(0.0)
                  D_JJ[off_ao_J:off_ao_J+nbf_J,off_ao_I:off_ao_I+nbf_I].fill(0.0)

           F_I = []

           for i in range(natom_I):
               xi = self._aggregate.all.x(off_natom_I+i)
               yi = self._aggregate.all.y(off_natom_I+i)
               zi = self._aggregate.all.z(off_natom_I+i)
               ri = numpy.array([xi,yi,zi])

               ints = self._electric_field_ao_integrals[off_natom_I+i] # ints_ALL (at I)

               f = self.__compute_electric_field_due_to_fragment(JJ_list, D_JJ, ints, ri)
               if field is not None: f+= field

               F_I.append(f)

           F_I = numpy.array(F_I)
           F.append(F_I)

           # print mean field values
           ffx = F_I[:,0]
           ffy = F_I[:,1]
           ffz = F_I[:,2]
           ff  = numpy.sqrt(ffx*ffx + ffy*ffy + ffz*ffz).mean()
           psi4.core.print_out(" Average E-Field on fragment %d: %14.6f [A.U.]\n" % (I+1,ff))



       # Compute perturbations: CT force
       W_2 = {}
       W_3 = {}
       W_4 = {}
       W_5 = {}
       W_6 = {}
       W_7 = {}

       # 2-body
       if do_23_body:
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              off_ao_I = self._ao_offsets_by_fragment[I]
                                                                                             
              D_I = self._Da[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].copy()
                                                                                             
              W_2[(I,I)] = numpy.zeros((nbf_I, nbf_I))
                                                                                             
              for J in range(self._number_of_fragments):
                if I!=J:
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                                                                                             
                  S_IJ = self._S[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                             
                  W_2[(I,J)] = D_I @ S_IJ
                                                                                             
          # 3-body
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              for J in range(self._number_of_fragments):
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                                                                                             
                  w_IJ = numpy.zeros((nbf_I, nbf_J))
                  #print("Considering IJ = %d%d" % (I,J))
                  for K in range(self._number_of_fragments):
                    if K!=I and K!=J : #and I==J:
                      #print("Adding K=%d to %dK%d" % (K,I,J))
                                                                                             
                      nbf_K = self._fragments[K].get_nbf()
                      off_ao_K = self._ao_offsets_by_fragment[K]
                                                                                             
                      S_KJ = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                             
                      w_IJ+= W_2[(I,K)] @ S_KJ
                                                                                             
                  W_3[(I,J)] = w_IJ

       if do_45_body:
          # 4-body                                                                                                  
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              for J in range(self._number_of_fragments):
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                  w_IJ = numpy.zeros((nbf_I, nbf_J))
                                                                                                                    
                  for K in range(self._number_of_fragments):
                   if I!=K:
                      nbf_K = self._fragments[K].get_nbf()
                      off_ao_K = self._ao_offsets_by_fragment[K]
                                                                                                                    
                      W_IK = W_2[(I,K)]
                                                                                                                    
                      for L in range(self._number_of_fragments):
                        if K!=L and L!=J: #and I==L and K==J:
                         #if I!=J and I==L and K==J:
                          nbf_L = self._fragments[L].get_nbf()
                          off_ao_L = self._ao_offsets_by_fragment[L]
                                                                                                                    
                          S_KL = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_L:off_ao_L+nbf_L].copy()
                          S_LJ = self._S[off_ao_L:off_ao_L+nbf_L,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                                                    
                          w_IJ+= W_IK @ S_KL @ S_LJ
                                                                                                                    
                      #if K!=I and K!=J and I!=J:
                      #   S_KJ = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_J:off_ao_J+nbf_J].copy()
                      #   w_IJ+= W_IK @ S_KJ
                                                                                                                    
                  W_4[(I,J)] = w_IJ
                                                                                                                    
          # 5-body
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              for J in range(self._number_of_fragments):
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                                                                                                                    
                  w_IJ = numpy.zeros((nbf_I, nbf_J))
                                                                                                                    
                  for K in range(self._number_of_fragments):
                   if I!=K:
                      nbf_K = self._fragments[K].get_nbf()
                      off_ao_K = self._ao_offsets_by_fragment[K]
                                                                                                                    
                      W_IK = W_2[(I,K)]
                                                                                                                    
                      for L in range(self._number_of_fragments):
                       if L!=K:
                          nbf_L = self._fragments[L].get_nbf()                                   
                          off_ao_L = self._ao_offsets_by_fragment[L]
                                                                                                                    
                          S_KL = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_L:off_ao_L+nbf_L].copy()
                                                                                                                    
                          for M in range(self._number_of_fragments):
                            if L!=M and M!=J: #  and I==L==J  and K==M: #and J==I==L:#  and I==L and K==M and L==J:
                             #if I==J==L and K==M:
                                                                                                                    
                              nbf_M = self._fragments[M].get_nbf()                                   
                              off_ao_M = self._ao_offsets_by_fragment[M]
                                                                                                     
                              S_LM = self._S[off_ao_L:off_ao_L+nbf_L,off_ao_M:off_ao_M+nbf_M].copy()
                              S_MJ = self._S[off_ao_M:off_ao_M+nbf_M,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                                    
                              w_IJ+= W_IK @ S_KL @ S_LM @ S_MJ
                                                                                                                    
                  W_5[(I,J)] = w_IJ



       if do_67_body:
          # 6-body                                                                                       
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              for J in range(self._number_of_fragments):
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                                                                                                         
                  w_IJ = numpy.zeros((nbf_I, nbf_J))
                                                                                                         
                  for K in range(self._number_of_fragments):
                   if I!=K:
                      nbf_K = self._fragments[K].get_nbf()
                      off_ao_K = self._ao_offsets_by_fragment[K]
                                                                                                         
                      W_IK = W_2[(I,K)]
                                                                                                         
                      for L in range(self._number_of_fragments):
                       if L!=K:
                          nbf_L = self._fragments[L].get_nbf()                                   
                          off_ao_L = self._ao_offsets_by_fragment[L]
                                                                                                         
                          S_KL = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_L:off_ao_L+nbf_L].copy()
                                                                                                         
                          for M in range(self._number_of_fragments):
                           if L!=M:
                            nbf_M = self._fragments[M].get_nbf()                                   
                            off_ao_M = self._ao_offsets_by_fragment[M]
                                                                                                         
                            S_LM = self._S[off_ao_L:off_ao_L+nbf_L,off_ao_M:off_ao_M+nbf_M].copy()
                                                                                                         
                            for N in range(self._number_of_fragments):
                             if M!=N and N!=J:
                                                                                                         
                              nbf_N = self._fragments[N].get_nbf()                                   
                              off_ao_N = self._ao_offsets_by_fragment[N]
                                                                                                     
                              S_MN = self._S[off_ao_M:off_ao_M+nbf_M,off_ao_N:off_ao_N+nbf_N].copy()
                              S_NJ = self._S[off_ao_N:off_ao_N+nbf_N,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                                    
                              w_IJ+= W_IK @ S_KL @ S_LM @ S_MN @ S_NJ
                                                                                                         
                  W_6[(I,J)] = w_IJ
                                                                                                         
          # 7-body
          for I in range(self._number_of_fragments):
              nbf_I = self._fragments[I].get_nbf()
              for J in range(self._number_of_fragments):
                  nbf_J = self._fragments[J].get_nbf()
                  off_ao_J = self._ao_offsets_by_fragment[J]
                                                                                                         
                  w_IJ = numpy.zeros((nbf_I, nbf_J))
                                                                                                         
                  for K in range(self._number_of_fragments):
                   if I!=K:
                      nbf_K = self._fragments[K].get_nbf()
                      off_ao_K = self._ao_offsets_by_fragment[K]
                                                                                                         
                      W_IK = W_2[(I,K)]
                                                                                                         
                      for L in range(self._number_of_fragments):
                       if L!=K:
                          nbf_L = self._fragments[L].get_nbf()                                   
                          off_ao_L = self._ao_offsets_by_fragment[L]
                                                                                                         
                          S_KL = self._S[off_ao_K:off_ao_K+nbf_K,off_ao_L:off_ao_L+nbf_L].copy()
                                                                                                         
                          for M in range(self._number_of_fragments):
                           if L!=M:
                            nbf_M = self._fragments[M].get_nbf()                                   
                            off_ao_M = self._ao_offsets_by_fragment[M]
                                                                                                         
                            S_LM = self._S[off_ao_L:off_ao_L+nbf_L,off_ao_M:off_ao_M+nbf_M].copy()
                                                                                                         
                            for N in range(self._number_of_fragments):
                             if N!=M:
                              nbf_N = self._fragments[N].get_nbf()                                   
                              off_ao_N = self._ao_offsets_by_fragment[N]
                                                                                                         
                              S_MN = self._S[off_ao_M:off_ao_M+nbf_M,off_ao_N:off_ao_N+nbf_N].copy()
                              
                              for O in range(self._number_of_fragments):
                               if O!=N and O!=J:
                                                                                                         
                                nbf_O = self._fragments[O].get_nbf()                                    
                                off_ao_O = self._ao_offsets_by_fragment[O]
                                                                                                       
                                S_NO = self._S[off_ao_N:off_ao_N+nbf_N,off_ao_O:off_ao_O+nbf_O].copy()
                                S_OJ = self._S[off_ao_O:off_ao_O+nbf_O,off_ao_J:off_ao_J+nbf_J].copy()
                                                                                                      
                                w_IJ+= W_IK @ S_KL @ S_LM @ S_MN @ S_NO @ S_OJ
                                                                                                         
                  W_7[(I,J)] = w_IJ



       # Compute blocks of deformation density matrix, dD
       # Reconstruct OPDM for entire system
       D_new = numpy.zeros((self._nbf, self._nbf))

       for I in range(self._number_of_fragments):
           nbf_I = self._fragments[I].get_nbf()
           off_ao_I = self._ao_offsets_by_fragment[I]

           dD_II_extField = numpy.einsum("acix,ix->ac"    , B_10[I], F[I])                     
           dD_II_extField+= numpy.einsum("acixy,ix,iy->ac", B_20[I], F[I], F[I])

           D_new[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I] = D0[I].copy() + dD_II_extField

           for J in range(I+1):
               nbf_J = self._fragments[J].get_nbf()
               off_ao_J = self._ao_offsets_by_fragment[J]

               IJ = (I,J)
               JI = (J,I)

               dD_IJ = D_new[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J].copy()

               # 2-body terms
               if I!=J:
                if B_01[J] is not None:           
                   dD_IJ+= W_2[IJ] @ B_01[J] 
                if B_01[I] is not None:
                   dD_IJ+= B_01[I].T @ W_2[JI].T

               # 3-body terms
               if 1: #I==J:
                if B_02[J] is not None:                                                    
                   dD_IJ+= W_3[IJ] @ B_02[J]
                if B_02[I] is not None:
                   dD_IJ+= B_02[I].T @ W_3[JI].T

               # 4-body terms
               if 1: #I!=J:
                if B_03[J] is not None:                                      
                   dD_IJ+= W_4[IJ] @ B_03[J] 
                if B_03[I] is not None:
                   dD_IJ+= B_03[I].T @ W_4[JI].T                            

               # 5-body terms
               if 1: #I==J:
                if B_04[J] is not None:          
                   dD_IJ+= W_5[IJ] @ B_04[J]
                if B_04[I] is not None:
                   dD_IJ+= B_04[I].T @ W_5[JI].T

               # 6-body terms
               if 1:   
                if B_05[J] is not None:          
                   dD_IJ+= W_6[IJ] @ B_05[J]
                if B_05[I] is not None:
                   dD_IJ+= B_05[I].T @ W_6[JI].T

               # 7-body terms
               if 1:   
                if B_06[J] is not None:          
                   dD_IJ+= W_7[IJ] @ B_06[J]
                if B_06[I] is not None:
                   dD_IJ+= B_06[I].T @ W_7[JI].T


               D_new[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J] = dD_IJ.copy()
               D_new[off_ao_J:off_ao_J+nbf_J,off_ao_I:off_ao_I+nbf_I] = dD_IJ.T.copy()

       # save
       self._Da = D_new.copy()
       if converge_by == 'energy':
          # Construct Fock matrix                                                        
          II = numpy.identity(self._nbf)
          self._jk.C_clear()                                           
          self._jk.C_left_add(psi4.core.Matrix.from_array(D_new, ""))
          self._jk.C_right_add(psi4.core.Matrix.from_array(II, ""))
          self._jk.compute()
          H = self._H
          J = self._jk.J()[0].to_array(dense=True)
          K = self._jk.K()[0].to_array(dense=True)
          G = 2.0 * J - K
          F = H + G
          # save
          self._Fa = F.copy()
          # compute energy
          E = (D_new @ (H + F)).trace() + self._aggregate.all.nuclear_repulsion_energy()
          if field is not None:
             mu_nuc = self._aggregate.all.nuclear_dipole()
             for x in range(3): E-= mu_nuc[x] * field[x]
                                                                                         
          return E
       elif converge_by == 'dipole':
          return self.dipole_moment()
       elif converge_by == 'opdm':
          return D_new
       else:
          raise NotImplementedError("Convergence by only energy, dipole or opdm is implemented")

   def __superslice(self, A, x, y, A_new=None, add=False):
       X, Y = numpy.meshgrid(x, y, indexing='ij')
       if A_new is None: 
          return A[X,Y]
       else: 
          if add: A[X,Y]+= A_new
          else:   A[X,Y] = A_new.copy()

   def __dmsscf_iteration_fragment(self, I, field):
       "1 Global Iteration of DMS SCF Procedure"

       # Extract OPDM and DMS tensors
       nat = [x.get_natoms() for x in self._fragments]
       par = [x.get() for x in self._fragments]
       D0  = [x['opdm'  ] for x in par]
       B_10= [x['dms_10'] for x in par] 
       B_20= [x['dms_20'] for x in par] 
       try: B_01= [x['dms_01'] for x in par] 
       except KeyError: B_01= [None for x in par] 
       try: B_02= [x['dms_02'] for x in par] 
       except KeyError: B_02= [None for x in par] 
       try: 
           B_03= [x['dms_03'] for x in par] 
       except KeyError: B_03= [None for x in par] 
       try: 
           B_04= [x['dms_04'] for x in par] 
       except KeyError: B_04= [None for x in par] 

       # AO spaces
       nbf_I = self._fragments[I].get_nbf()
       off_ao_I = self._ao_offsets_by_fragment[I]
       I_idx = numpy.arange(off_ao_I, off_ao_I+nbf_I)

       idx = numpy.array([], dtype=int)
       for J in range(self._number_of_fragments):
         if I!=J:
           nbf_J = self._fragments[J].get_nbf()
           off_ao_J = self._ao_offsets_by_fragment[J]
           J_idx = numpy.arange(off_ao_J, off_ao_J+nbf_J)
           idx = numpy.hstack([idx, J_idx])

       # B-tensor of entire aggregate
       B_01_all = numpy.zeros((self._nbf, self._nbf))
       B_02_all = numpy.zeros((self._nbf, self._nbf))
       B_03_all = numpy.zeros((self._nbf, self._nbf))
       B_04_all = numpy.zeros((self._nbf, self._nbf))
       for J in range(self._number_of_fragments):
           nbf_J = self._fragments[J].get_nbf()
           off_ao_J = self._ao_offsets_by_fragment[J]
           J_idx = numpy.arange(off_ao_J, off_ao_J+nbf_J)
           if B_01[J] is not None: self.__superslice(B_01_all, J_idx, J_idx, B_01[J]) 
           if B_02[J] is not None: self.__superslice(B_02_all, J_idx, J_idx, B_02[J]) 
           if B_03[J] is not None: self.__superslice(B_03_all, J_idx, J_idx, B_03[J]) 
           if B_04[J] is not None: self.__superslice(B_04_all, J_idx, J_idx, B_04[J]) 

       # Compute perturbations: Electric field
       F = []

       for II in range(self._number_of_fragments):
           nbf_II = self._fragments[II].get_nbf()                                           
           natom_II = nat[II]
           off_ao_II = self._ao_offsets_by_fragment[II]
           off_natom_II = self._natom_offsets_by_fragment[II]
                                                                                          
           D_JJ = self._Da.copy()
           D_JJ[off_ao_II:off_ao_II+nbf_II,off_ao_II:off_ao_II+nbf_II].fill(0.0)
                                                                                          
           JJ_list = []
           for J in range(self._number_of_fragments):
               if II!=J: 
                  JJ_list.append(J)
                  ## local field approximation:
                  #off_ao_J = self._ao_offsets_by_fragment[J]
                  #nbf_J = self._fragments[J].get_nbf()
                  #D_JJ[off_ao_II:off_ao_II+nbf_II,off_ao_J:off_ao_J+nbf_J].fill(0.0)
                  #D_JJ[off_ao_J:off_ao_J+nbf_J,off_ao_II:off_ao_II+nbf_II].fill(0.0)
                                                                                          
           F_II = []
                                                                                          
           for i in range(natom_II):
               xi = self._aggregate.all.x(off_natom_II+i)
               yi = self._aggregate.all.y(off_natom_II+i)
               zi = self._aggregate.all.z(off_natom_II+i)
               ri = numpy.array([xi,yi,zi])
                                                                                          
               ints = self._electric_field_ao_integrals[off_natom_II+i] # ints_ALL (at I)
                                                                                          
               f = self.__compute_electric_field_due_to_fragment(JJ_list, D_JJ, ints, ri)
               if field is not None: f+= field
                                                                                          
               F_II.append(f)

           F.append(numpy.array(F_II))



       # Compute perturbations: PT force
       D_I  = self.__superslice(self._Da, I_idx, I_idx)
       D_J  = self.__superslice(self._Da,   idx, idx)
       S_IJ = self.__superslice(self._S , I_idx, idx)
       S_JI = S_IJ.T.copy()

       W_IJ = D_I @ S_IJ       
       W_JI = D_J @ S_JI 
                               
       W_IJI = W_IJ @ S_JI
       W_JIJ = W_JI @ S_IJ
                               
       W_IJIJ = W_IJI @ S_IJ
       W_JIJI = W_JIJ @ S_JI
                               
       W_IJIJI = W_IJIJ @ S_JI
       W_JIJIJ = W_JIJI @ S_IJ
       #

       # Compute blocks of deformation density matrix, dD
       # Reconstruct OPDM for entire system
       D_new = self._Da.copy() #numpy.zeros((self._nbf, self._nbf))

       # ---> electric fields
       for J in range(self._number_of_fragments):
       # if J!=I:
           nbf_J = self._fragments[J].get_nbf()
           off_ao_J = self._ao_offsets_by_fragment[J]
           J_idx = numpy.arange(off_ao_J, off_ao_J+nbf_J)

           dD_JJ_extField = numpy.einsum("acix,ix->ac"    , B_10[J], F[J])                     
           dD_JJ_extField+= numpy.einsum("acixy,ix,iy->ac", B_20[J], F[J], F[J])

           self.__superslice(D_new, J_idx, J_idx, D0[J].copy() + dD_JJ_extField, add=False)
        #else: self.__superslice(D_new, J_idx, J_idx, self.__superslice(self._Da, J_idx, J_idx), add=False)

       # ---> wavefunction overlaps
       dD_II = numpy.zeros((len(I_idx),len(I_idx)))
       dD_JJ = numpy.zeros((len(  idx),len(  idx)))
       dD_IJ = numpy.zeros((len(I_idx),len(  idx)))

       # B_I <-- W_others
       if B_01[I] is not None:
          dD_IJ+= B_01[I].T @ W_JI.T
       if B_03[I] is not None:
          dD_IJ+= B_03[I].T @ W_JIJI.T

       if B_02[I] is not None:
          dD_II+= W_IJI @ B_02[I] + B_02[I].T @ W_IJI.T
       if B_04[I] is not None:
          dD_II+= W_IJIJI @ B_04[I] + B_04[I].T @ W_IJIJI.T

       # W_I --> B_others
       B_01_J = self.__superslice(B_01_all, idx, idx)
       B_02_J = self.__superslice(B_02_all, idx, idx)
       B_03_J = self.__superslice(B_03_all, idx, idx)
       B_04_J = self.__superslice(B_04_all, idx, idx)

       dD_IJ+= W_IJ @ B_01_J
       dD_IJ+= W_IJIJ @ B_03_J

       dD_JJ+= W_JIJ @ B_02_J + B_02_J.T @ W_JIJ.T
       dD_JJ+= W_JIJIJ @ B_04_J + B_04_J.T @ W_JIJIJ.T

       self.__superslice(D_new, I_idx,I_idx, dD_II  , add=True)
       self.__superslice(D_new,   idx,  idx, dD_JJ  , add=True)
       self.__superslice(D_new, I_idx,  idx, dD_IJ  , add=False)
       self.__superslice(D_new,   idx,I_idx, dD_IJ.T, add=False)


       # Construct Fock matrix
       II = numpy.identity(self._nbf)
       self._jk.C_clear()                                           
       self._jk.C_left_add(psi4.core.Matrix.from_array(D_new, ""))
       self._jk.C_right_add(psi4.core.Matrix.from_array(II, ""))
       self._jk.compute()
       H = self._H
       J = self._jk.J()[0].to_array(dense=True)
       K = self._jk.K()[0].to_array(dense=True)
       G = 2.0 * J - K
       F = H + G
       # compute total energy
       E = (D_new @ (H + F)).trace() + self._aggregate.all.nuclear_repulsion_energy()
       if field is not None:
          mu_nuc = self._aggregate.all.nuclear_dipole()
          for x in range(3): E-= mu_nuc[x] * field[x]

       # save
       self._Da = D_new.copy()
       self._Fa = F.copy()
       return E


   def __dmsscf_iteration_older(self):
       "1 Global Iteration of DMS SCF Procedure"

       # Extract OPDM and DMS tensors
       nat = [x.get_natoms() for x in self._fragments]
       par = [x.get() for x in self._fragments]
       D0  = [x['opdm'  ] for x in par]
       B_10= [x['dms_10'] for x in par] 
       B_20= [x['dms_20'] for x in par] 
       try: B_01= [x['dms_01'] for x in par] 
       except KeyError: B_01= [None for x in par] 
       try: B_02= [x['dms_02'] for x in par] 
       except KeyError: B_02= [None for x in par] 
       try: B_03= [x['dms_03'] for x in par] 
       except KeyError: B_03= [None for x in par] 
       try: B_04= [x['dms_04'] for x in par] 
       except KeyError: B_04= [None for x in par] 

       # Compute perturbations: Electric field
       F = []
       for I in range(self._number_of_fragments):
           nbf_I = self._fragments[I].get_nbf()
           natom_I = nat[I]
           off_ao_I = self._ao_offsets_by_fragment[I]
           off_natom_I = self._natom_offsets_by_fragment[I]

           JJ_list = []
           for J in range(self._number_of_fragments):
               if I!=J: JJ_list.append(J)

           F_I = []
           for i in range(natom_I):
               xi = self._aggregate.all.x(off_natom_I+i)
               yi = self._aggregate.all.y(off_natom_I+i)
               zi = self._aggregate.all.z(off_natom_I+i)
               ri = numpy.array([xi,yi,zi])

               D_JJ = self._Da.copy()
               D_JJ[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].fill(0.0)

               ints = self._electric_field_ao_integrals[off_natom_I+i] # ints_ALL (at I)

               f = self.__compute_electric_field_due_to_fragment(JJ_list, D_JJ, ints, ri)
               F_I.append(f)
           F.append(numpy.array(F_I))

       # Compute perturbations: CT force
       W_IJ = {}
       W_JI = {}
       W_IJI = {}
       W_JIJ = {}
       W_IJIJ = {}
       W_JIJI = {}
       W_IJIJI = {}
       W_JIJIJ = {}

       for I in range(self._number_of_fragments):
           nbf_I = self._fragments[I].get_nbf()
           off_ao_I = self._ao_offsets_by_fragment[I]

           D_I = self._Da[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].copy()

           for J in range(self._number_of_fragments):
             if I!=J:
               nbf_J = self._fragments[J].get_nbf()
               off_ao_J = self._ao_offsets_by_fragment[J]

               D_J = self._Da[off_ao_J:off_ao_J+nbf_J,off_ao_J:off_ao_J+nbf_J].copy()

               S_IJ = self._S[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J].copy()
               S_JI = S_IJ.T.copy()

               pair = (I,J)

               w_IJ = D_I @ S_IJ       
               w_JI = D_J @ S_JI 
                                       
               w_IJI = w_IJ @ S_JI
               w_JIJ = w_JI @ S_IJ
                                       
               w_IJIJ = w_IJI @ S_IJ
               w_JIJI = w_JIJ @ S_JI
                                       
               w_IJIJI = w_IJIJ @ S_JI
               w_JIJIJ = w_JIJI @ S_IJ
               #
               W_IJ[pair] = w_IJ
               W_JI[pair] = w_JI
               W_IJI[pair] = w_IJI
               W_JIJ[pair] = w_JIJ
               W_IJIJ[pair] = w_IJIJ
               W_JIJI[pair] = w_JIJI
               W_IJIJI[pair] = w_IJIJI
               W_JIJIJ[pair] = w_JIJIJ

       # Compute blocks of deformation density matrix, dD
       # Reconstruct OPDM for entire system
       D_new = self._Da.copy()
       # dD_II
       for I in range(self._number_of_fragments):
           nbf_I = self._fragments[I].get_nbf()
           off_ao_I = self._ao_offsets_by_fragment[I]

           dD_II = numpy.einsum("acix,ix->ac"    , B_10[I], F[I])                     
           dD_II+= numpy.einsum("acixy,ix,iy->ac", B_20[I], F[I], F[I])
           for J in range(self._number_of_fragments):
             if I!=J:
               nbf_J = self._fragments[J].get_nbf()
               off_ao_J = self._ao_offsets_by_fragment[J]

               IJ = (I,J)
               if B_02[I] is not None:                                                  
                  dD_II+= 1.0 * (W_IJI[IJ] @ B_02[I] + B_02[I].T @ W_IJI[IJ].T)
               if B_04[I] is not None: 
                  dD_II+= 1.0 * (W_IJIJI[IJ] @ B_04[I] + B_04[I].T @ W_IJIJI[IJ].T)

               dD_IJ = numpy.zeros((nbf_I, nbf_J))
               if B_01[J] is not None:
                  dD_IJ+= W_IJ[IJ] @ B_01[J] 
               if B_01[I] is not None:
                  dD_IJ+= B_01[I].T @ W_JI[IJ].T                            
               if B_03[J] is not None:
                  dD_IJ+= W_IJIJ[IJ] @ B_03[J] 
               if B_03[I] is not None:
                  dD_IJ+= B_03[I].T @ W_JIJI[IJ].T                            

               D_new[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J] = dD_IJ

           D_new[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I] = D0[I].copy() + dD_II

       # Construct Fock matrix
       II = numpy.identity(self._nbf)
       self._jk.C_clear()                                           
       self._jk.C_left_add(psi4.core.Matrix.from_array(D_new, ""))
       self._jk.C_right_add(psi4.core.Matrix.from_array(II, ""))
      #self._jk.C_left_add(wfn.Da())
      #self._jk.C_right_add(psi4.core.Matrix.from_array(II, ""))
       self._jk.compute()
       H = self._H
       J = self._jk.J()[0].to_array(dense=True)
       K = self._jk.K()[0].to_array(dense=True)
       G = 2.0 * J - K
       F = H + G
       #J_test = jk.J()[1].to_array(dense=True)
       #K_test = jk.K()[1].to_array(dense=True)
       #G_test = 2.0 * J_test - K_test
       #F_test = H + G_test
       #D_test = wfn.Da().to_array(dense=True)
       # compute energy
       E = (D_new @ (H + F)).trace() + self._aggregate.all.nuclear_repulsion_energy()
       #E_test = (D_test @ (H + F_test)).trace() + wfn.molecule().nuclear_repulsion_energy()

       # save
       self._Da = D_new.copy()
       self._Fa = F.copy()
       return E

   def __dmsscf_for_pair_of_fragments(self, I, J, field):
       "DMS SCF Procedure for Pair: Original code written for dimer case - kept for debugging"
       frg_I = self._fragments[I]; par_I = frg_I.get()
       frg_J = self._fragments[J]; par_J = frg_J.get()

       natom_I = frg_I.get_natoms(); nbf_I = frg_I.get_nbf()
       natom_J = frg_J.get_natoms(); nbf_J = frg_J.get_nbf()

       off_ao_I = self._ao_offsets_by_fragment[I]
       off_ao_J = self._ao_offsets_by_fragment[J]
       off_natom_I = self._natom_offsets_by_fragment[I]
       off_natom_J = self._natom_offsets_by_fragment[J]

       # Extract OPDM and DMS tensors
       D_A0 = par_I['opdm']
       D_B0 = par_J['opdm']
       BA_10 = par_I['dms_10']; BB_10 = par_J['dms_10']
       BA_20 = par_I['dms_20']; BB_20 = par_J['dms_20']
       try: BA_01 = par_I['dms_01']; BB_01 = par_J['dms_01']
       except: BA_01, BB_01 = None, None
       try: BA_02 = par_I['dms_02']; BB_02 = par_J['dms_02']
       except: BA_02, BB_02 = None, None
       try: BA_03 = par_I['dms_03']; BB_03 = par_J['dms_03']
       except: BA_03, BB_03 = None, None
       try: BA_04 = par_I['dms_04']; BB_04 = par_J['dms_04']
       except: BA_04, BB_04 = None, None

       # Compute perturbations
       D_A = self._Da[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].copy()
       D_B = self._Da[off_ao_J:off_ao_J+nbf_J,off_ao_J:off_ao_J+nbf_J].copy()
       S_AB = self._S[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J].copy()
       S_BA = S_AB.T.copy()
       W_AB = D_A @ S_AB
       W_BA = D_B @ S_BA 

       W_ABA = W_AB @ S_BA
       W_BAB = W_BA @ S_AB

       W_ABAB = W_ABA @ S_AB
       W_BABA = W_BAB @ S_BA

       W_ABABA = W_ABAB @ S_BA
       W_BABAB = W_BABA @ S_AB

       F_B = [] # field due to B evaluated on A atoms
       JJ_list = []
       for JJ in range(self._number_of_fragments):
           if JJ!=I: JJ_list.append(JJ)
       for i in range(natom_I):
        xi = self._aggregate.all.x(off_natom_I+i)
        yi = self._aggregate.all.y(off_natom_I+i)
        zi = self._aggregate.all.z(off_natom_I+i)
        ri = numpy.array([xi,yi,zi])


        D_JJ = self._Da.copy()
        D_JJ[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I].fill(0.0)

        ints = self._electric_field_ao_integrals[off_natom_I+i] # ints_B (at A)

        f = self.__compute_electric_field_due_to_fragment(JJ_list, D_JJ, ints, ri)
        if field is not None: f+= field
        F_B.append(f)
           
       F_A = [] # field due to A evaluated on B atoms
       II_list = []
       for II in range(self._number_of_fragments):
           if II!=J: II_list.append(II)
       for j in range(natom_J):
        xj = self._aggregate.all.x(off_natom_J+j)
        yj = self._aggregate.all.y(off_natom_J+j)
        zj = self._aggregate.all.z(off_natom_J+j)
        rj = numpy.array([xj,yj,zj])

        D_II = self._Da.copy()
        D_II[off_ao_J:off_ao_J+nbf_J,off_ao_J:off_ao_J+nbf_J].fill(0.0)

        ints = self._electric_field_ao_integrals[off_natom_J+j,:] # ints_A (at B)

        f = self.__compute_electric_field_due_to_fragment(II_list, D_II, ints, rj)
        if field is not None: f+= field

        F_A.append(f)
       F_A = numpy.array(F_A); F_B = numpy.array(F_B)

       # Compute blocks of deformation density matrix, dD
       # dD_AA
       dD_AA = numpy.einsum("acix,ix->ac", BA_10, F_B)
       dD_AA+= numpy.einsum("acixy,ix,iy->ac", BA_20, F_B, F_B)
       if BA_02 is not None: 
          dD_AA+= W_ABA @ BA_02 + BA_02.T @ W_ABA.T                      
       if BA_04 is not None: 
          dD_AA+= W_ABABA @ BA_04 + BA_04.T @ W_ABABA.T                      
       # dD_BB
       dD_BB = numpy.einsum("acix,ix->ac", BB_10, F_A)
       dD_BB+= numpy.einsum("acixy,ix,iy->ac", BB_20, F_A, F_A)
       if BB_02 is not None:
          dD_BB+= W_BAB @ BB_02 + BB_02.T @ W_BAB.T                      
       if BB_04 is not None:
          dD_BB+= W_BABAB @ BB_04 + BB_04.T @ W_BABAB.T                      
       # dD_AB
       dD_AB = numpy.zeros((nbf_I, nbf_J))
       if BA_01 is not None:
          dD_AB+= W_AB @ BB_01 + BA_01.T @ W_BA.T                            
       if BA_03 is not None:
          dD_AB+= W_ABAB @ BB_03 + BA_03.T @ W_BABA.T                            

       # Reconstruct OPDM for entire system
       D_new = self._Da.copy()
       D_new[off_ao_I:off_ao_I+nbf_I,off_ao_I:off_ao_I+nbf_I] = D_A0+ dD_AA
       D_new[off_ao_J:off_ao_J+nbf_J,off_ao_J:off_ao_J+nbf_J] = D_B0+ dD_BB
       D_new[off_ao_I:off_ao_I+nbf_I,off_ao_J:off_ao_J+nbf_J] =       dD_AB
       D_new[off_ao_J:off_ao_J+nbf_J,off_ao_I:off_ao_I+nbf_I] =       dD_AB.T.copy()

       # Fock matrix
       II = numpy.identity(self._nbf)
       self._jk.C_clear()                                           
       self._jk.C_left_add(psi4.core.Matrix.from_array(D_new, ""))
       self._jk.C_right_add(psi4.core.Matrix.from_array(II, ""))
       self._jk.compute()
       H = self._H
       J = self._jk.J()[0].to_array(dense=True)
       K = self._jk.K()[0].to_array(dense=True)
       G = 2.0 * J - K
       F = H + G
       # compute energy
       E = (D_new @ (H + F)).trace() + self._aggregate.all.nuclear_repulsion_energy()
       if field is not None:
          mu_nuc = self._aggregate.all.nuclear_dipole()
          for x in range(3): E-= mu_nuc[x] * field[x]

       # save
       self._Da = D_new.copy()
       self._Fa = F.copy()
       return E


class Property_Computer(_DMS_SCF_Procedure):
   "Implements properties other than total energy of the system."
   def __init__(self, aggregate, fragments):
       _DMS_SCF_Procedure.__init__(self, aggregate, fragments)

   def atomic_charges(self, kappa=0.0):
       "Compute atomic partial charges"
       raise NotImplementedError

   def camm(self):
       "Compute ground-state cumulative atomic multipole moments"
       wfn = psi4.core.Wavefunction(self._aggregate.all, self._bfs)
       Da  = psi4.core.Matrix.from_array(self._Da, "")
       camm = oepdev.DMTPole.build("CAMM", wfn, 1)
       camm.compute([Da], [False])
       return camm

   def dipole_moment(self):
       "Compute dipole moment from nuclear dipole moment and density matrix"
       mu = self._aggregate.all.nuclear_dipole()
       mu = numpy.array([mu[0], mu[1], mu[2]])
       DIP = self._mints.ao_dipole()
       mu += 2.0 * numpy.array([(DIP[x].to_array(dense=True) @ self._Da).trace() for x in range(3)])
       return mu

   def ff_dipole_moment(self, step=0.001, verbose=None):
       "Finite-field dipole moment: differentiation of total energy"
       self.run(field=numpy.array([0.0,0.0,step]), verbose=verbose)
       Ez1 = self.E()
       self.run(field=numpy.array([0.0,0.0,-step]), verbose=verbose)
       Ezm1= self.E()
       #
       self.run(field=numpy.array([0.0,step,0.0]), verbose=verbose)
       Ey1 = self.E()
       self.run(field=numpy.array([0.0,-step,0.0]), verbose=verbose)
       Eym1= self.E()
       #
       self.run(field=numpy.array([step,0.0,0.0]), verbose=verbose)
       Ex1 = self.E()
       self.run(field=numpy.array([-step,0.0,0.0]), verbose=verbose)
       Exm1= self.E()
       #
       h = 2.0 * step
       mu_x =-(Ex1 - Exm1)/ h
       mu_y =-(Ey1 - Eym1)/ h
       mu_z =-(Ez1 - Ezm1)/ h
       mu = numpy.array([mu_x, mu_y, mu_z])
       return mu 

   def ff_polarizability(self, step=0.001, energy=False, verbose=None):
       "Finite-field dipole-dipole polarizability"
       # differentiate total energy wrt electric field
       if energy:
          self.run(field=numpy.array([0.0,0.0,0.0]), verbose=verbose)
          E_0 = self.E()
          #
          self.run(field=numpy.array([step,0.0,0.0]), verbose=verbose)
          E_1 = self.E()
          self.run(field=numpy.array([-step,0.0,0.0]), verbose=verbose)
          E_2 = self.E()
          #
          self.run(field=numpy.array([0.0,step,0.0]), verbose=verbose)
          E_3 = self.E()
          self.run(field=numpy.array([0.0,-step,0.0]), verbose=verbose)
          E_4 = self.E()
          #
          self.run(field=numpy.array([0.0,0.0,step]), verbose=verbose)  
          E_5 = self.E()
          self.run(field=numpy.array([0.0,0.0,-step]), verbose=verbose)
          E_6 = self.E()
          ##
          self.run(field=numpy.array([ step, step,0.0]), verbose=verbose)  
          E_11= self.E()
          self.run(field=numpy.array([-step, step,0.0]), verbose=verbose)
          E_12= self.E()
          self.run(field=numpy.array([-step,-step,0.0]), verbose=verbose)  
          E_13= self.E()
          self.run(field=numpy.array([ step,-step,0.0]), verbose=verbose)
          E_14= self.E()
          #
          self.run(field=numpy.array([ step,0.0, step]), verbose=verbose)  
          E_21= self.E()                             
          self.run(field=numpy.array([-step,0.0, step]), verbose=verbose)
          E_22= self.E()                             
          self.run(field=numpy.array([-step,0.0,-step]), verbose=verbose)  
          E_23= self.E()                             
          self.run(field=numpy.array([ step,0.0,-step]), verbose=verbose)
          E_24= self.E()
          #
          self.run(field=numpy.array([0.0, step, step]), verbose=verbose)  
          E_31= self.E()                  
          self.run(field=numpy.array([0.0,-step, step]), verbose=verbose)
          E_32= self.E()                  
          self.run(field=numpy.array([0.0,-step,-step]), verbose=verbose)  
          E_33= self.E()                  
          self.run(field=numpy.array([0.0, step,-step]), verbose=verbose)
          E_34= self.E()

          h = step
          a_xx = ((E_1-E_0) + (E_2-E_0)) / h**2
          a_yy = ((E_3-E_0) + (E_4-E_0)) / h**2
          a_zz = ((E_5-E_0) + (E_6-E_0)) / h**2

          a_xy = ((E_11-E_12) + (E_13-E_14)) / (4.0 * h**2)
          a_xz = ((E_21-E_22) + (E_23-E_24)) / (4.0 * h**2)
          a_yz = ((E_31-E_32) + (E_33-E_34)) / (4.0 * h**2)

          alpha = -numpy.array([[a_xx, a_xy, a_xz],[a_xy, a_yy, a_yz],[a_xz, a_yz, a_zz]])

       # differentiate total dipole moment wrt electric field
       else:
          self.run(field=numpy.array([0.0,0.0,step]), verbose=verbose)  
          Ez1 = self.dipole_moment()
          self.run(field=numpy.array([0.0,0.0,-step]), verbose=verbose)
          Ezm1= self.dipole_moment()
          #
          self.run(field=numpy.array([0.0,step,0.0]), verbose=verbose)
          Ey1 = self.dipole_moment()
          self.run(field=numpy.array([0.0,-step,0.0]), verbose=verbose)
          Eym1= self.dipole_moment()
          #
          self.run(field=numpy.array([step,0.0,0.0]), verbose=verbose)
          Ex1 = self.dipole_moment()
          self.run(field=numpy.array([-step,0.0,0.0]), verbose=verbose)
          Exm1= self.dipole_moment()
          #
          h = 2.0 * step
          alpha_Ox = (Ex1 - Exm1)/ h
          alpha_Oy = (Ey1 - Eym1)/ h
          alpha_Oz = (Ez1 - Ezm1)/ h
          alpha = numpy.vstack([alpha_Ox, alpha_Oy, alpha_Oz])
       return alpha

   def ff_hyperpolarizability(self, step=0.001, energy=False, verbose=None):
       "Finite-field dipole-dipole first hyperpolarizability"
       # differentiate total energy wrt electric field
       if energy:
          raise NotImplementedError
          self.run(field=numpy.array([0.0,0.0,0.0]), verbose=verbose)
          E0 = self.E()
       # differentiate total dipole moment wrt electric field
       else:
          raise NotImplementedError
       return beta

   def ff_hyper2polarizability(self, step=0.001, energy=False, verbose=None):
       "Finite-field dipole-dipole second hyperpolarizability"
       # differentiate total energy wrt electric field
       if energy:
          raise NotImplementedError
          self.run(field=numpy.array([0.0,0.0,0.0]), verbose=verbose)
          E0 = self.E()
       # differentiate total dipole moment wrt electric field
       else:
          raise NotImplementedError
       return gamma

          

class SimpleComputer(Property_Computer):
   def __init__(self, *args, **kwargs):
       Property_Computer.__init__(self, *args, **kwargs)

   def _orthogonalize_density_matrix(self):
       # nothing to do here
       pass

class OtherComputer(Property_Computer): # --> deprecate
   def __init__(self, *args, **kwargs):
       Property_Computer.__init__(self, *args, **kwargs)

   def _orthogonalize_density_matrix(self):
       if self._n is not None:
          n, c = Density.natural_orbitals(self._Da, self._S, self._Ca, orthogonalize_mo = True, 
                                                    n_eps = 0.001, ignore_large_n= True,
                                                    no_cutoff=0.000, renormalize=False)
       else:
          n = numpy.zeros(self._nbf); 
          for I in range(self._number_of_fragments):
              naocc = self._fragments[I].get()['naocc']
              off = self._ao_offsets_by_fragment[I]
              for i in range(naocc):
                  n[off+i] = 1.0
          s = self._Ca.T @ self._S @ self._Ca
          x = matrix_power(s,-0.5)
          c = self._Ca @ x
         #self._Ca = c.copy()

       n += 0.000000000001

       w = numpy.sqrt(n); print(n)
       s_mo = c.T @ self._S @ c

       wsw = numpy.diag(w) @ s_mo @ numpy.diag(w)
       wswm12 = matrix_power(wsw, -0.5)
       K = wswm12 @ numpy.diag(n) @ wswm12
       cw = c @ numpy.diag(w)
       Da_oo = cw @ K @ cw.T
      #print(Da_oo); exit()

       self._Da = Da_oo
       self._Ca = c
       self._n = n


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
  def __init__(self, typ):
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
      self._Ca_occ   = None           # LCAO-MO Occ Alpha
      self._Ca_vir   = None           # LCAO-MO Vir Alpha

      self._init(typ)
      self._available_orders = None   # Implemented orders of DMS tensors

  @classmethod
  def create(cls, dms_type='da', order_type='basic'):
      if   order_type.lower() == 'field'   : return          Field_DMS(dms_type)
      elif order_type.lower() == 'basic'   : return          Basic_DMS(dms_type)
      elif order_type.lower() == 'bd'      : return         BasicD_DMS(dms_type)
      elif order_type.lower() == 'new'     : return            New_DMS(dms_type)
      elif order_type.lower() == 'new2'    : return           New2_DMS(dms_type)
      else:
         raise NotImplementedError
      

  def available_orders(self): return self._available_orders

  def N(self): return self._N
  def n(self): return self._n
  def M(self): return self._M.copy()
  def B(self, m, n): 
      if (m,n) in self._available_orders:
          return self._B[(m,n)]
      else:
          raise ValueError("DMS Error: B(%i,%i) is not available" % (m,n))

  def set_bfs(self, bfs):
      self._bfs = bfs
      self._mol = bfs.molecule()
      self._n   = bfs.nbf()
      self._N   = self._mol.natom()
  def set_M(self, M): self._M = M.copy()
  def set_s1(self, s): self._s1 = s.copy()
  def set_s2(self, s): self._s2 = s.copy()
  def set_ca(self, Ca_occ, Ca_vir):
      self._Ca_occ = Ca_occ.copy()
      self._Ca_vir = Ca_vir.copy()
  def Ca_occ(self): return self._Ca_occ.copy()
  def Ca_vir(self): return self._Ca_vir.copy()
  def M(self): return self._M.copy()

  def B(self, m, n):
      "Retrieve DMS tensor from parameter vector"
      order = (m, n)
      if order in self._available_orders:
         if order in self._B.keys(): 
            return self._B[order]
         else:
            self._generate_B_from_s(order)
      else:
         raise ValueError("This DMS object does not include order (%i,%i)" %(m,n))

  def read(self, out):#TODO
      "Create DMS object from file"
      raise NotImplementedError

  def write(self, out):#TODO
      "Write DMS object to file"
      raise NotImplementedError

  def _generate_B_from_s(self, order):
      "Retrieve DMS tensor from its parameter vector"
      n  = self._n
      n2 = self._n**2 
      N  = self._N
      n_ = composite.number_of_elements(n)

      # ---> Group Z(1) <--- #
      # (1,0)
      if (1,0) in self._available_orders and order == (1,0):
          START = 0
          END   = START + n_*N*3
          b = self._s1[START:END].reshape(n_,N,3)
          self._B[(1,0)] = composite.retrieve_tensor(b)

      # (2,0)
      if (2,0) in self._available_orders and order == (2,0):
          START = n_*N*3
          END   = START + n_*N*6
          b     = self._s1[START:END].reshape(n_,N,6)
          b     = composite.retrieve_tensor(b)
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
          self._B[(2,0)] = B

      # (0,2)
      if (0,2) in self._available_orders and order == (0,2):
          START = n_*(N*3  + N*6)
          END   = START + n2
          self._B[(0,2)] = self._s1[START:END].reshape(n,n)

      # (0,4)
      if (0,4) in self._available_orders and order == (0,4):
          START = n_*(N*3  + N*6) + n2
          END   = START + n2
          self._B[(0,4)] = self._s1[START:END].reshape(n,n)

      # (0,6)
      if (0,6) in self._available_orders and order == (0,6):
          START = n_*(N*3  + N*6) + n2 + n2
          END   = START + n2
          self._B[(0,6)] = self._s1[START:END].reshape(n,n)


      # ---> Group Z(2) <--- #
      # (0,1)
      if (0,1) in self._available_orders and order == (0,1):
          START = 0
          END   = START + n2
          self._B[(0,1)] = self._s2[START:END].reshape(n,n)

      # (0,3)
      if (0,3) in self._available_orders and order == (0,3):
          START = n2
          END   = START + n2
          self._B[(0,3)] = self._s2[START:END].reshape(n,n)

      # (0,5)
      if (0,5) in self._available_orders and order == (0,5):
          START = n2 + n2
          END   = START + n2
          self._B[(0,5)] = self._s2[START:END].reshape(n,n)


  def _init(self, t):
      "Initialize information"
      u = {"d": "Density", "g": "Fock", "a": "Alpha", "b": "Beta"}
      m = [x for x in t.lower()]
      self._type      = t.lower()
      self._type_long = u[m[0]] 
      if len(m)>1: self._type_long+= '-' + u[m[1]]


class Field_DMS(DMS):
  """
 Basic model of DMS that handles induction up to second-order
 in the external electric field without any other electronic densities.

 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
"""
  def __init__(self, type='da'):
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0)]


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
  def __init__(self, type='da'):
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1)]


class BasicD_DMS(DMS):
  """
 Basic model of DMS that handles induction up to second-order and
 Pauli effects up to secod-order.

 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,2)          Pauli          Asymmetric    (n,n)         Z(1)
  4.            (0,1)          Pauli          Asymmetric    (n,n)         Z(2)
"""
  def __init__(self, type='da'):
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,2), (0,1)]

class New_DMS(DMS):
  """
 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,1)          Pauli          Asymmteric    (n,n)         Z(1)
  4.            (0,3)          Pauli          Asymmetric    (n,n)         Z(1)
  5.            (0,2)          Pauli          Asymmetric    (n,n)         Z(2)
  6.            (0,4)          Pauli          Asymmetric    (n,n)         Z(2)
"""
  def __init__(self, type='da'):
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1), (0,2), (0,3), (0,4)]


class New2_DMS(DMS):
  """
 The order of DMS blocks:

 Block         DMS order      Interaction    Symmetry       Dimension     Group
 -----         ---------      -----------    ------------  -----------   -------
  1.            (1,0)          Induction      Symmetric     (n,n,N,3)     Z(1)
  2.            (2,0)          Induction      Symmetric     (n,n,N,3,3)   Z(1)
  3.            (0,1)          Pauli          Asymmteric    (n,n)         Z(1)
  4.            (0,3)          Pauli          Asymmetric    (n,n)         Z(1)
  5.            (0,5)          Pauli          Asymmetric    (n,n)         Z(1)
  6.            (0,2)          Pauli          Asymmetric    (n,n)         Z(2)
  7.            (0,4)          Pauli          Asymmetric    (n,n)         Z(2)
  8.            (0,6)          Pauli          Asymmetric    (n,n)         Z(2)
"""
  def __init__(self, type='da'):
      DMS.__init__(self, type)

      self._available_orders = [(1,0), (2,0), (0,1), (0,2), (0,3), (0,4), (0,5), (0,6)]




class DMSFit(ABC):
   """
 Method to fit DMS tensors.
"""

   # ---> Global defaults <--- #
                                                                                                       
   minimum_atom_atom_distance       = 1.2 / psi4.constants.bohr2angstroms
   minimum_atom_charge_distance     = 4.5
   srange                           = 4.0
   stop_when_error_hessian          = True
   compute_dmatpol_susceptibilities = False
   generate_random_samples          = True
   read_samples                     = False
   dmatpol_ntest_charge             = 40
   dmatpol_test_charge              = 0.05
   dmatpol_esp_pad_shpere           = 6.0 # Bohr
   dmatpol_gradient_rank            = 0
   replace_extfield_with_dmatpol    = False
                                                                                                       
   def __init__(self, mol, method, 
                      nsamples, dms_types, order_type, use_iterative_model, use_external_field_model):
       ABC.__init__(self)
                                                                                                       
       if mol.multiplicity() > 1:
          raise NotImplementedError(" DMSFit Error: Open-shell systems are not implemented yet.")
                                                                                                       
       # Molecule
       self._mol       = mol
       # QM Method
       self._method    = method
       # Number of Fitting Samples
       self._nsamples  = nsamples
       # Types of DMS to Fit
       self._dms_types = dms_types.split(',')
       # DMS Model
       self._order_type= order_type
       # Iterative Model
       self._use_iterative_model = use_iterative_model
       # External Field Model
       self._use_external_field_model = use_external_field_model
                                                                                                       
                                                                                                       
   @classmethod
   def create(cls, mol, fit_type="transl", dms_types="da,g", order_type='basic',
                   nsamples=100, method='scf', 
                   use_iterative_model=True, use_external_field_model=False):
       if fit_type.lower().startswith("tran"): 
          return Translation_DMSFit(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)
       elif fit_type.lower().startswith("rot"): 
          return Rotation_DMSFit(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)
       else: 
          raise NotImplementedError("This type of DMS fitting is not implemented yet")


   # ---> Abstract methods <--- #

   @abstractmethod
   def run(self): pass

   @abstractmethod
   def B(self): pass

   @abstractmethod
   def _compute_samples(self): pass

   @abstractmethod
   def _check(self): pass


class _Global_Settings_DMSFit(DMSFit):
   """
 Global settings and utilities to fit DMS tensors.
"""
   def __init__(self, mol, method,                                                         
                      nsamples, dms_types, order_type, use_iterative_model, use_external_field_model):
       super().__init__(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)
                                                                                           
       # Number of atoms under DMS fitting
       self._natoms    = self._mol.natom()
                                                                                           
       # Current ID of statistical sample
       self._i         = 0


   # ---> Protected Utilities <--- #

   def _invert_hessian(self, h):                                                                                
       "Invert Hessian matrix"                                                                               
       det= numpy.linalg.det(h)
       psi4.core.print_out(" ** Hessian Determinant= %14.6E\n" % det)
       hi = numpy.linalg.inv(h)
       I = numpy.dot(hi, h).diagonal().sum()
       d = I - len(hi)
       psi4.core.print_out(" ** Hessian Dimension= (%8d) x (%8d)\n" % h.shape)
       psi4.core.print_out(" ** delta=%14.6E %14.2f %d\n" % (d,I,h.shape[0]))
       if DMSFit.stop_when_error_hessian:
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
                                                                                                                
       fi[0] += 2.0 * (D_j @ ints_j[0]).trace() 
       fi[1] += 2.0 * (D_j @ ints_j[1]).trace() 
       fi[2] += 2.0 * (D_j @ ints_j[2]).trace() 
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
                                                                                                                
   def _save_xyz(self, mol, out_name, misc=None):
       "Save molecule mol to xyz file out_name with additional information misc"
       out=open(out_name,'w')
       log = mol.to_string('xyz', units='Angstrom')
       if misc is not None:
          log = log.split('\n')
          log[1] = "%s" % misc
          log = '\n'.join(log)
       out.write(log+"\n"); out.close()

   def _save_frg(self, dms, out_name):
       "Save Solvshift FRG file with DMS tensors"
       frg = Fragment()
       frg.set(basis=psi4.core.get_global_option("BASIS"), method=self._method, mol=self._mol, dms=dms)
       frg.write(out_name)

   def _save_dD(self, D, prefix='dd'):          
       "Save matrix to a file"
       name= prefix + '_%03d.dat' 
       out = open(name % self._i, 'w')
       log = ''
       n = len(D)
       for row in D:
           log+= n*"%13.5E" % tuple(row) + "\n"
       out.write(log)
       out.close() 

   def _rms(self, a, b): 
       "RMS between arrays a and b"
       return numpy.sqrt(((a-b)**2).sum()/a.size)

   def _random_double(self):
       "Random double between -1 and 1"
       return 2.0*numpy.random.random() - 1.0





class EFP_DMSFit(_Global_Settings_DMSFit):
  """
 DMS Fitting for EFP-like Molecular Fragments
"""

  def __init__(self, mol, method, 
                     nsamples, dms_types, order_type, use_iterative_model, use_external_field_model):
      super().__init__(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)

      psi4.set_options({"DMATPOL_TRAINING_MODE":"CHARGES",
                        "DMATPOL_FIELD_RANK"   : 2,
                        "DMATPOL_GRADIENT_RANK": DMSFit.dmatpol_gradient_rank,
                        "DMATPOL_NSAMPLES"     : nsamples,
                        "DMATPOL_NTEST_CHARGE" : DMSFit.dmatpol_ntest_charge,
                        "ESP_PAD_SPHERE"       : DMSFit.dmatpol_esp_pad_shpere,
    "cphf_diis"                    : True,
    "cphf_diis_dim"                : 8,
    "cphf_maxiter"                 : 200,
    "cphf_conver"                  : 1e-8,
    "cphf_localize"                : True,
    "cphf_localizer"               : "PIPEK_MEZEY",
    "save_jk"                      : True,
                        "DMATPOL_TEST_CHARGE"  : DMSFit.dmatpol_test_charge})

      self._e0 = None
      self._wfn_0= None
      self._bfs_0= None
      self._nbf = None
      self._dms_da = None
      self._dms_g  = None
      self._D0  = None
      self._G0  = None

      self._dimer_wfn = None
      self._dimer_mol = None
      self._dimer_mints = None

      # DMS for external electric field (JCP 2018)
      self._B_ind_10 = None
      self._B_ind_20 = None
      self._dms_extField_da = None
      self._dms_extField_g = None
      self._dms_extField_da_dmatpol = None



  # ---> Implementation <--- #

  def run(self):#TODO
      "Run the fitting procedure"
      # compute quantities for isolated molecule
      self.__compute_isolated_molecule()

      s = self._nsamples
      n = self._dms_da.n()
      N = self._dms_da.N()

      # compute DMS for isolated molecule in external electric field from DMATPOL routine (Da only)
      if DMSFit.compute_dmatpol_susceptibilities:
         self.__compute_dms_external_field_dmatpol()

      # prepare for the DMS fittings
      if not DMSFit.read_samples:
         self._compute_samples()

      # compute DMS for isolated molecule in external electric field
      dms_extField_da = None
      dms_extField_g  = None
      if self._use_external_field_model:
         F_extField_set = numpy.fromfile('temp_F_extField_set.dat').reshape(s,N,3)
         dD_extField_ref_set = numpy.fromfile('temp_dD_extField_ref_set.dat').reshape(s,n,n)
         dG_extField_ref_set = numpy.fromfile('temp_dG_extField_ref_set.dat').reshape(s,n,n)
         self._compute_group_extField(self._dms_extField_da, F_extField_set, dD_extField_ref_set)
         self._compute_group_extField(self._dms_extField_g , F_extField_set, dG_extField_ref_set)
         dms_extField_da = self._dms_extField_da
         dms_extField_g  = self._dms_extField_g 
         self._save_frg(dms_extField_da, 'dms_extField_da.frg') # save before substitution
         self._save_frg(dms_extField_g , 'dms_extField_g.frg' )
         del F_extField_set, dD_extField_ref_set, dG_extField_ref_set
         # replace with DMATPOL susceptibilities (OPDM only)
         if DMSFit.replace_extfield_with_dmatpol and DMSFit.compute_dmatpol_susceptibilities:
            self._dms_extField_da = self._dms_extField_da_dmatpol
            dms_extField_da = self._dms_extField_da

      # fit DMS Z-1: (1,0), (2,0), (0,2), (1,2) and (2,2)
      F_A_set  = numpy.fromfile('temp_F_A_set.dat').reshape(s,N,3)
      F_B_set  = numpy.fromfile('temp_F_B_set.dat').reshape(s,N,3)

      dD_AA_ref_set = numpy.fromfile('temp_dD_AA_ref_set.dat').reshape(s,n,n)
      dD_BB_ref_set = numpy.fromfile('temp_dD_BB_ref_set.dat').reshape(s,n,n)

      W_AB_set = numpy.fromfile('temp_W_AB_set.dat').reshape(s,n,n)
      W_BA_set = numpy.fromfile('temp_W_BA_set.dat').reshape(s,n,n)
      A_AB_set = numpy.fromfile('temp_A_AB_set.dat').reshape(s,n,n)
      A_BA_set = numpy.fromfile('temp_A_BA_set.dat').reshape(s,n,n)
      
      self._compute_group_1(self._dms_da, dms_extField_da,
                            F_A_set, F_B_set, dD_AA_ref_set, dD_BB_ref_set, A_AB_set, A_BA_set, W_AB_set, W_BA_set)
      del dD_AA_ref_set, dD_BB_ref_set

      dG_AA_ref_set = numpy.fromfile('temp_dG_AA_ref_set.dat').reshape(s,n,n)
      dG_BB_ref_set = numpy.fromfile('temp_dG_BB_ref_set.dat').reshape(s,n,n)

      w_AB_set = numpy.fromfile('temp_w_AB_set.dat').reshape(s,n,n)
      w_BA_set = numpy.fromfile('temp_w_BA_set.dat').reshape(s,n,n)
      a_AB_set = numpy.fromfile('temp_a_AB_set.dat').reshape(s,n,n)
      a_BA_set = numpy.fromfile('temp_a_BA_set.dat').reshape(s,n,n)


      self._compute_group_1(self._dms_g , dms_extField_g,
                            F_A_set, F_B_set, dG_AA_ref_set, dG_BB_ref_set, a_AB_set, a_BA_set, w_AB_set, w_BA_set)
      del dG_AA_ref_set, dG_BB_ref_set


      # fit DMS Z-2: (0,1), (1,1) and (2,1)
      K_set    = numpy.fromfile('temp_K_set.dat').reshape(s,n,n)
      L_set    = numpy.fromfile('temp_L_set.dat').reshape(s,n,n)
      dD_AB_ref_set = numpy.fromfile('temp_dD_AB_ref_set.dat').reshape(s,n,n)

      self._compute_group_2(self._dms_da, K_set, L_set, dD_AB_ref_set,
                                          W_AB_set, W_BA_set, A_AB_set, A_BA_set)
      del K_set, L_set,                   W_AB_set, W_BA_set, A_AB_set, A_BA_set

      k_set    = numpy.fromfile('temp_k_set.dat').reshape(s,n,n)
      l_set    = numpy.fromfile('temp_l_set.dat').reshape(s,n,n)
      dG_AB_ref_set = numpy.fromfile('temp_dG_AB_ref_set.dat').reshape(s,n,n)

      self._compute_group_2(self._dms_g , k_set, l_set, dG_AB_ref_set,
                                          w_AB_set, w_BA_set, a_AB_set, a_BA_set)
      del k_set, l_set, F_A_set, F_B_set, w_AB_set, w_BA_set, a_AB_set, a_BA_set

      # save DMS
      self._save_frg(self._dms_da, "dms_full_da.frg")
      self._save_frg(self._dms_g, "dms_full_g.frg")

      # test the DMS fittings
      self._check()

  def B(self, m, n, dtype='da'):
      "Get the susceptibilities"
      if   dtype.lower() == 'da': return self._dms_da.B(m,n)
      elif dtype.lower() == 'g' : return self._dms_g .B(m,n)
      elif 'field' in dtype.lower():
           if dtype.lower().startswith('da'):  return self._dms_extField_da.B(m,n)
           elif dtype.lower().startswith('g'): return self._dms_extField_g .B(m,n)
           else: raise ValueError("DMSFit Error: Wrong value!!")
      elif dtype.lower() == 'dmatpol': 
           if   (m,n) == (1,0): return self._B_ind_10
           elif (m,n) == (2,0): return self._B_ind_20
           else: 
               raise ValueError("DMSFit Error: Only 10 and 20 electric field susceptibilities are available")

      else: raise NotImplementedError("DMSFit Error: Only 'da', 'g' and 'fa' available now.")



  # ---> Protected Interface <--- #
  
  def _determine_perturbing_densities(self):                                                       
      "Determine the perturbations depending on the type of DMS model"
      if not self._use_iterative_model:
         DA = self._dms_da._M.copy()
         DB = DA.copy()
         #
         GA = self._dms_g ._M.copy()
         GB = GA.copy()
      else:
         D = self._dimer_wfn.Da().to_array(dense=True)
         DA= D[:self._nbf,:self._nbf]
         DB= D[self._nbf:,self._nbf:]
         #
         if   'e' in self._dms_types: # G = 2J - K
               G = self._dimer_wfn.Fa().to_array(dense=True) - self._dimer_wfn.H().to_array(dense=True)
         elif 'g' in self._dms_types: # G = 2J - K
              #T = self._dimer_mints.ao_kinetic().to_array(dense=True)
              #G = self._dimer_wfn.Fa().to_array(dense=True) - T
               G = self._dimer_wfn.Fa().to_array(dense=True) - self._dimer_wfn.H().to_array(dense=True)
         #elif 'f' in self._dms_types: # G = T + V + 2J - K = F
         #      G = self._dimer_wfn.Fa().to_array(dense=True).copy()
         #elif 'g1' in self._dms_types: # G = 2V + 2J - K 
         #      T = self._dimer_mints.ao_kinetic().to_array(dense=True)
         #      V = self._dimer_mints.ao_potential().to_array(dense=True)
         #      G = self._dimer_wfn.Fa().to_array(dense=True) - T + V
         #else:
         #      raise ValueError(" Only g or e or f or g1 types for Fock DMS are available.")
            
         GA= G[:self._nbf,:self._nbf]
         GB= G[self._nbf:,self._nbf:]
      #
      self._DA = DA
      self._DB = DB
      #
      self._GA = GA
      self._GB = GB


  def _compute_efield(self):
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

      s = 1.0
      return F_A_set*s, F_B_set*s, F_A_mat_set, F_B_mat_set

   
  # ---> Private Interface <--- #

  def __compute_isolated_molecule(self):
      psi4.core.print_out(" ---> Computing Unperturbed Wavefunction <---\n\n")

      # Compute HF or DFT wavefunction of molecule of interest
      self._e0, self._wfn_0 = PSI4_DRIVER(self._method, molecule=self._mol, return_wfn=True)

      # Set the AO basis set
      self._bfs_0 = self._wfn_0.basisset()
      self._nbf = self._wfn_0.basisset().nbf()

      Da = self._wfn_0.Da().to_array(dense=True)
      Ca_occ = self._wfn_0.Ca_subset('AO','OCC').to_array(dense=True)
      Ca_vir = self._wfn_0.Ca_subset('AO','VIR').to_array(dense=True)
      #
      G_extField = self._wfn_0.Fa().to_array(dense=True) - self._wfn_0.H().to_array(dense=True) # G = 2J - K
      if   'e' in self._dms_types: # G = 2J - K
            G = G_extField.copy()
      elif 'g' in self._dms_types: # G = 2J - K
            G = G_extField.copy()
           #mints = psi4.core.MintsHelper(self._bfs_0)
           #T = mints.ao_kinetic().to_array(dense=True)
           #G = self._wfn_0.Fa().to_array(dense=True) - T
      #elif 'f' in self._dms_types: # G = T + V + 2J - K = F
      #      G = self._wfn_0.Fa().to_array(dense=True).copy()
      #elif 'g1' in self._dms_types: # G = 2V + 2J - K
      #      mints = psi4.core.MintsHelper(self._bfs_0)
      #      T = mints.ao_kinetic().to_array(dense=True)
      #      V = mints.ao_potential().to_array(dense=True)
      #      G = self._wfn_0.Fa().to_array(dense=True) - T + V


      # Initialize DMS tensors
      self._dms_da= DMS.create('da', self._order_type)
      self._dms_g = DMS.create('g' , self._order_type)

      self._dms_da.set_bfs(self._bfs_0)
      self._dms_g .set_bfs(self._bfs_0)

      self._dms_da.set_M(Da)
      self._dms_g .set_M(G )

      self._dms_da.set_ca(Ca_occ, Ca_vir)
      self._dms_g .set_ca(Ca_occ, Ca_vir)

      self._D0 = Da.copy()
      self._G0 = G .copy()

      if self._use_external_field_model:
         self._dms_extField_da = DMS.create('da', "Field") 
         self._dms_extField_g  = DMS.create('g' , "Field")
         self._dms_extField_da.set_bfs(self._bfs_0)
         self._dms_extField_g .set_bfs(self._bfs_0)
         self._dms_extField_da.set_M(Da.copy())
         self._dms_extField_g .set_M(G_extField.copy())
         self._dms_extField_da.set_ca(Ca_occ, Ca_vir)
         self._dms_extField_g .set_ca(Ca_occ, Ca_vir)

      #if   self._dms_type.lower() == 'da': M = self._wfn_0.Da().to_array(dense=True)
      #elif self._dms_type.lower() == 'db': M = self._wfn_0.Db().to_array(dense=True)
      #elif self._dms_type.lower() == 'fa': M = self._wfn_0.Fa().to_array(dense=True)
      #elif self._dms_type.lower() == 'fb': M = self._wfn_0.Fb().to_array(dense=True)
      #else:
      #    raise valueerror(" DMSFit Error: Incorrect type of dms tensors. Available: da, db, fa, fb.")
      psi4.core.clean()


  def __compute_dms_external_field_dmatpol(self):
      "DMS for external electric field only (without other molecules) from DMATPOL routine"
      solver = oepdev.GenEffParFactory.build("POLARIZATION", self._wfn_0, psi4.core.get_options())
      genEffPar = solver.compute()

      n = self._dms_da.n()
      N = self._dms_da.N()

      B_ind_10 = numpy.zeros((n,n,N,3))
      B_ind_20 = numpy.zeros((n,n,N,3,3))

      for i in range(N):
          B_10_ix = genEffPar.susceptibility(1,0,i,0).to_array(dense=True) 
          B_10_iy = genEffPar.susceptibility(1,0,i,1).to_array(dense=True) 
          B_10_iz = genEffPar.susceptibility(1,0,i,2).to_array(dense=True) 

          B_20_ixx= genEffPar.susceptibility(2,0,i,0).to_array(dense=True) # xx 0
          B_20_ixy= genEffPar.susceptibility(2,0,i,1).to_array(dense=True) # xy 1
          B_20_ixz= genEffPar.susceptibility(2,0,i,2).to_array(dense=True) # xz 2
          B_20_iyy= genEffPar.susceptibility(2,0,i,4).to_array(dense=True) # yx 3
          B_20_iyz= genEffPar.susceptibility(2,0,i,5).to_array(dense=True) # yy 4 
          B_20_izz= genEffPar.susceptibility(2,0,i,8).to_array(dense=True) # yz 5 zx 6 zy 7 zz 8

          B_ind_10[:,:,i,0] = B_10_ix.copy()
          B_ind_10[:,:,i,1] = B_10_iy.copy()
          B_ind_10[:,:,i,2] = B_10_iz.copy()

          B_ind_20[:,:,i,0,0] = B_20_ixx.copy()
          B_ind_20[:,:,i,0,1] = B_20_ixy.copy()
          B_ind_20[:,:,i,0,2] = B_20_ixz.copy()
          B_ind_20[:,:,i,1,1] = B_20_iyy.copy()
          B_ind_20[:,:,i,1,2] = B_20_iyz.copy()
          B_ind_20[:,:,i,2,2] = B_20_izz.copy()
          B_ind_20[:,:,i,1,0] = B_20_ixy.copy()
          B_ind_20[:,:,i,2,0] = B_20_ixz.copy()
          B_ind_20[:,:,i,2,1] = B_20_iyz.copy()


      self._B_ind_10 = B_ind_10 
      self._B_ind_20 = B_ind_20 

      # create DMS object
      s1 = composite.form_futf(B_ind_10).ravel()
      B2 = numpy.array([B_ind_20[:,:,:,0,0],B_ind_20[:,:,:,0,1],B_ind_20[:,:,:,0,2],
                        B_ind_20[:,:,:,1,1],B_ind_20[:,:,:,1,2],B_ind_20[:,:,:,2,2]]).transpose(1,2,3,0)
      s2 = composite.form_futf(B2).ravel()
      s = numpy.hstack([s1,s2])
      self._dms_extField_da_dmatpol = DMS.create('da', "Field") 
      self._dms_extField_da_dmatpol.set_bfs(self._bfs_0)
      self._dms_extField_da_dmatpol.set_M(self._D0.copy())
      self._dms_extField_da_dmatpol.set_s1(s)
      self._dms_extField_da_dmatpol._generate_B_from_s((1,0))
      self._dms_extField_da_dmatpol._generate_B_from_s((2,0))
      Ca_occ = self._wfn_0.Ca_subset('AO','OCC').to_array(dense=True)
      Ca_vir = self._wfn_0.Ca_subset('AO','VIR').to_array(dense=True)
      self._dms_extField_da_dmatpol.set_ca(Ca_occ, Ca_vir)
      self._save_frg(self._dms_extField_da_dmatpol, "dms_field_da_dmatpol.frg")

      # clean up
      psi4.core.clean()


  # ---> Abstract methods <--- #

  @abstractmethod
  def _compute_group_extField(self): pass

  @abstractmethod
  def _compute_group_1(self): pass

  @abstractmethod
  def _compute_group_2(self): pass

  @abstractmethod
  def _construct_aggregate(self): pass

  @abstractmethod
  def _extField(self): pass



class ExternalField_EFP_DMSFit(EFP_DMSFit):

  def __init__(self, mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model):
      super().__init__(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)

  # ---> Implementation <--- #

  def _compute_group_extField(self, dms, F_set, dM_ref_set):
      "Compute B(10) and B(20) in external electric field"
      psi4.core.print_out(" ---> Computing DMS for Ext-Field group of type %s <---\n\n" % dms._type_long)

      n = dms.n()
      N = dms.N()

      # First compute Hessian                                                                                                
      dim_1 = N*3
      dim_2 = N*6
      H = numpy.zeros((dim_1 + dim_2, dim_1 + dim_2), numpy.float64)
                                                                                            
      H_10_10 = numpy.einsum("niu,njw->iujw", F_set, F_set) 
                                                                                            
      H[:dim_1,:dim_1] = H_10_10.copy().reshape(dim_1, dim_1)
      del H_10_10
                                                                                            
      H_20_20 = numpy.einsum("niu,nix,njw,njy->iuxjwy", F_set, F_set, F_set, F_set)
                                                                                            
      u = numpy.zeros((N,6,N,6))
      u[:,0,:,0] = H_20_20[:,0,0,:,0,0].copy() * 1.0
      u[:,0,:,1] = H_20_20[:,0,0,:,0,1].copy() * 2.0
      u[:,0,:,2] = H_20_20[:,0,0,:,0,2].copy() * 2.0
      u[:,0,:,3] = H_20_20[:,0,0,:,1,1].copy() * 1.0
      u[:,0,:,4] = H_20_20[:,0,0,:,1,2].copy() * 2.0
      u[:,0,:,5] = H_20_20[:,0,0,:,2,2].copy() * 1.0
                                                                                            
      u[:,1,:,0] = H_20_20[:,0,1,:,0,0].copy() * 2.0
      u[:,1,:,1] = H_20_20[:,0,1,:,0,1].copy() * 4.0
      u[:,1,:,2] = H_20_20[:,0,1,:,0,2].copy() * 4.0
      u[:,1,:,3] = H_20_20[:,0,1,:,1,1].copy() * 2.0
      u[:,1,:,4] = H_20_20[:,0,1,:,1,2].copy() * 4.0
      u[:,1,:,5] = H_20_20[:,0,1,:,2,2].copy() * 2.0
                                                                                            
      u[:,2,:,0] = H_20_20[:,0,2,:,0,0].copy() * 2.0
      u[:,2,:,1] = H_20_20[:,0,2,:,0,1].copy() * 4.0
      u[:,2,:,2] = H_20_20[:,0,2,:,0,2].copy() * 4.0
      u[:,2,:,3] = H_20_20[:,0,2,:,1,1].copy() * 2.0
      u[:,2,:,4] = H_20_20[:,0,2,:,1,2].copy() * 4.0
      u[:,2,:,5] = H_20_20[:,0,2,:,2,2].copy() * 2.0
                                                                                            
      u[:,3,:,0] = H_20_20[:,1,1,:,0,0].copy() * 1.0
      u[:,3,:,1] = H_20_20[:,1,1,:,0,1].copy() * 2.0
      u[:,3,:,2] = H_20_20[:,1,1,:,0,2].copy() * 2.0
      u[:,3,:,3] = H_20_20[:,1,1,:,1,1].copy() * 1.0
      u[:,3,:,4] = H_20_20[:,1,1,:,1,2].copy() * 2.0
      u[:,3,:,5] = H_20_20[:,1,1,:,2,2].copy() * 1.0
                                                                                            
      u[:,4,:,0] = H_20_20[:,1,2,:,0,0].copy() * 2.0
      u[:,4,:,1] = H_20_20[:,1,2,:,0,1].copy() * 4.0
      u[:,4,:,2] = H_20_20[:,1,2,:,0,2].copy() * 4.0
      u[:,4,:,3] = H_20_20[:,1,2,:,1,1].copy() * 2.0
      u[:,4,:,4] = H_20_20[:,1,2,:,1,2].copy() * 4.0
      u[:,4,:,5] = H_20_20[:,1,2,:,2,2].copy() * 2.0
                                                                                            
      u[:,5,:,0] = H_20_20[:,2,2,:,0,0].copy() * 1.0
      u[:,5,:,1] = H_20_20[:,2,2,:,0,1].copy() * 2.0
      u[:,5,:,2] = H_20_20[:,2,2,:,0,2].copy() * 2.0
      u[:,5,:,3] = H_20_20[:,2,2,:,1,1].copy() * 1.0
      u[:,5,:,4] = H_20_20[:,2,2,:,1,2].copy() * 2.0
      u[:,5,:,5] = H_20_20[:,2,2,:,2,2].copy() * 1.0
                                                                                            
      H[dim_1:,dim_1:] = u.copy().reshape(dim_2, dim_2)
      del u, H_20_20
                                                                                            
      H_10_20 = numpy.einsum("niu,njw,njy->iujwy", F_set, F_set, F_set)
                                                                                            
      u = numpy.zeros((N,3,N,6))
      u[:,:,:,0] = H_10_20[:,:,:,0,0].copy() * 1.0
      u[:,:,:,1] = H_10_20[:,:,:,0,1].copy() * 2.0
      u[:,:,:,2] = H_10_20[:,:,:,0,2].copy() * 2.0
      u[:,:,:,3] = H_10_20[:,:,:,1,1].copy() * 1.0
      u[:,:,:,4] = H_10_20[:,:,:,1,2].copy() * 2.0
      u[:,:,:,5] = H_10_20[:,:,:,2,2].copy() * 1.0
                                                                                            
      H[:dim_1,dim_1:] = u.copy().reshape(dim_1, dim_2)
      H[dim_1:,:dim_1] = H[:dim_1,dim_1:].copy().T
      del u, H_10_20
                                                                                            
      H *= 2.0
                                                                                                               
      # Fit for each matrix element separately
      psi4.core.print_out(" * Computing Hessian Inverse...\n")
      Hi = self._invert_hessian(H)
                                                                                            
      S_1 = numpy.zeros((n,n,N,3))
      S_2 = numpy.zeros((n,n,N,6))
                                                                                            
      for i in range(n):
          for j in range(n):
             #psi4.core.print_out(" * Computing the DMS tensors for (%i %i)...\n" % (i,j))
                                                                                            
              # gradient                                                           
              g_10 = numpy.einsum("n,niu->iu", dM_ref_set[:,i,j], F_set)
                                                                                   
              u = numpy.zeros((N,6))
              g_20 = numpy.einsum("n,niu,niw->iuw", dM_ref_set[:,i,j], F_set, F_set)
              u[:,0] = g_20[:,0,0].copy()
              u[:,1] = g_20[:,0,1].copy() * 2.0
              u[:,2] = g_20[:,0,2].copy() * 2.0
              u[:,3] = g_20[:,1,1].copy()
              u[:,4] = g_20[:,1,2].copy() * 2.0
              u[:,5] = g_20[:,2,2].copy()
                                                                                   
              g = numpy.zeros(dim_1 + dim_2)
                                                                                   
              g[:dim_1] = g_10.ravel().copy()
              g[dim_1:] = u.ravel().copy()
              del u, g_10, g_20
              g *= -2.0
                                                                                   
              s = - g @ Hi
                                                                                            
              S_1[i,j] = s[:dim_1].copy().reshape(N,3)
              S_2[i,j] = s[dim_1:].copy().reshape(N,6)
                                                                                                           
      S_1 = composite.form_futf(S_1)
      S_2 = composite.form_futf(S_2)
                                                                                            
      s = numpy.hstack([S_1.ravel(), S_2.ravel()])
                                                                                                               
                                                                                                              
      # Extract susceptibilities
      psi4.core.print_out(" * Setting the DMS tensors...\n")
      dms.set_s1(s)
      for order in [(1,0),(2,0)]:
          if order in dms.available_orders(): dms._generate_B_from_s(order)


  def _extField(self):
      "Compute electric field and wavefunction in the presence of point charges"
      psi4.core.clean()
      psi4.core.print_out(" ===> ExtField Run For Sample %i <===\n" % self._i)

      opt_stash = DMSFit.minimum_atom_atom_distance
      DMSFit.minimum_atom_atom_distance = DMSFit.minimum_atom_charge_distance
      charges = self.__generate_random_charges()
      DMSFit.minimum_atom_atom_distance = opt_stash

      F = self.__generate_field_from_charges(charges, self._mol)
      F_aver = numpy.linalg.norm(F, axis=1)
      for i in range(self._mol.natom()):
          psi4.core.print_out(" Sample %3d Field on Atom %2i = %14.8f [a.u.]\n" % (self._i, i+1, F_aver[i]))

      MM_charges = psi4.QMMM()
      BohrToAngstrom = 0.5291772086
      for charge in charges:
          q = charge[0]
          x = charge[1] * BohrToAngstrom
          y = charge[2] * BohrToAngstrom
          z = charge[3] * BohrToAngstrom
          MM_charges.extern.addCharge(q,x,y,z)
      psi4.core.set_global_option_python('EXTERN', MM_charges.extern)

      e, wfn = psi4.energy(self._method, molecule=self._mol, return_wfn=True)
      psi4.core.print_out(" Sample %3d ExtField Total Energy= %16.6f\n" % (self._i, e))

      V = MM_charges.extern.computePotentialMatrix(self._bfs_0).to_array(dense=True)
      E_nuc_field = self.__compute_nuclear_field_interaction_energy(self._mol, charges)

      MM_charges.extern.clear()

      psi4.core.set_global_option_python('EXTERN', None)
      return F, V, wfn, E_nuc_field


  # ---> Private interface <--- #

  def __generate_field_from_charges(self, charges, mol):
      fields = numpy.zeros((mol.natom(), 3))
      geom = mol.geometry().to_array(dense=True)
      for i in range(len(geom)):
          fields[i] = self.__field_due_to_charges(charges, geom[i])
      return fields

  def __field_due_to_charges(self, charges, r):
      x,y,z = r
      fx = 0.0; fy = 0.0; fz = 0.0
      for charge in charges:
          q,X,Y,Z = charge
          dx = x-X; dy = y-Y; dz = z-Z
          R = math.sqrt(dx*dx+dy*dy+dz*dz)
          R = q/(R*R*R)
          fx+= R * dx; fy+= R * dy; fz+= R * dz
      f = numpy.array([fx,fy,fz])
      return f

  def __draw_random_point(self):
      geom = self._mol.geometry().to_array(dense=True)
      com = self._mol.center_of_mass()
      cx_ = com[0]; cy_=com[1]; cz_ = com[2]
      pad = psi4.core.get_options().get_double("ESP_PAD_SPHERE")
      minval = geom[0,0] - cx_; maxval = minval
      for i in range(self._mol.natom()):
          for j in range(3):
              val = geom[i,j] - com[j]
              if val > maxval: maxval = val
              if val < minval: minval = val
      rad = max(abs(maxval), abs(minval))
      radius_ = rad + pad

      theta = numpy.arccos(self._random_double())
      phi   = 2.0 * numpy.pi * self._random_double()
      r     = radius_ * numpy.cbrt(self._random_double())
      x     = cx_ + r * numpy.sin(theta) * numpy.cos(phi)
      y     = cy_ + r * numpy.sin(theta) * numpy.sin(phi)
      z     = cz_ + r * numpy.cos(theta)                 

      point = numpy.array([x,y,z])
      return point

  def __generate_random_charges(self):
      geom = self._mol.geometry().to_array(dense=True)
      nc = psi4.core.get_options().get_int("DMATPOL_NTEST_CHARGE")
      charges = numpy.zeros((nc,4))

     #print(geom * 0.5291772086)
      for i in range(nc):
          done = False
          while not done:
              point = self.__draw_random_point()
              if not self._clash(geom, [point,]):
                 done = True 
                 q = self.__draw_random_charge()
                 x,y,z = point
                 charges[i] = numpy.array([q,x,y,z])
         #print("X %14.6f %14.6f %14.6f" % tuple(charges[i,1:] * 0.5291772086))
      return charges

  def __draw_random_charge(self):
      scale = psi4.core.get_options().get_double("DMATPOL_TEST_CHARGE")
      q = self._random_double() * scale
      return q

  def __compute_nuclear_field_interaction_energy(self, mol, charges):
      E = 0.0
      for i in range(self._mol.natom()):
          Zi = self._mol.Z(i)
          xi = self._mol.x(i)
          yi = self._mol.y(i)
          zi = self._mol.z(i)
          for j in range(len(charges)):
              qj, xj, yj, zj = charges[j]
              rij = math.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
              E += Zi * qj / rij
      return E


class Translation_DMSFit(ExternalField_EFP_DMSFit):
  """
 Translation method to fit DMS tensors.
"""
  def __init__(self, mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model):
      super().__init__(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)

      # translation parameters
      self._start = 0.0
      self._range = DMSFit.srange


  # ---> Implementation <--- #

  def _construct_aggregate(self):
      "Create next dimer by translating a molecule"
      psi4.core.clean()

      self._i += 1
      log = "\n%d %d\n" % (self._mol.molecular_charge(), self._mol.multiplicity())
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
      log += "%d %d\n" % (self._mol.molecular_charge(), self._mol.multiplicity())

      t = self.__new_translation()

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
     #print(log)

      mol = psi4.geometry(log)
      mol.update_geometry()

      e_dimer, wfn_dimer = PSI4_DRIVER(self._method, molecule=mol, return_wfn=True)
      self._dimer_wfn = wfn_dimer
      self._dimer_mol = mol
      self._dimer_bfs = wfn_dimer.basisset()
      self._dimer_mints = psi4.core.MintsHelper(self._dimer_bfs)

      DIP_nuc = self._dimer_mol.nuclear_dipole()
      DIP_nuc = [DIP_nuc[0], DIP_nuc[1], DIP_nuc[2],]
      parcel = (e_dimer, self._dimer_mol.nuclear_repulsion_energy(), DIP_nuc)

      self._mol_A = self._mol
      self._mol_B = mol.extract_subsets(2); 
      e_monomer, wfn_B = PSI4_DRIVER(self._method, molecule=self._mol_B, return_wfn=True)
      self._wfn_A = self._wfn_0
      self._wfn_B = wfn_B
      self._bfs_A = self._bfs_0
      self._bfs_B = wfn_B.basisset()

      e_int = e_dimer - 2.0 * e_monomer # MCBS
      message = " Sample %3d. Interaction energy= %13.6f [a.u.]   %13.6f [kcal/mol]" % (self._i, e_int, e_int*627.5 )
      psi4.core.print_out(message+'\n'); print(message)
      # 
      out = "geom_%03d.xyz" % self._i
      self._save_xyz(mol, out, misc="Eint(MCBS)=%13.3f [kcal/mol]" % (e_int*627.5))
      #
      psi4.core.clean()
      return mol, parcel



  def _compute_samples(self):
      "Compute electric fields, CT channels and reference deformation matrices"

      e_set = []
      E_nuc_set = []
      DIP_nuc_set = []

      dD_AA_set_ref = []
      dD_BB_set_ref = []
      dD_AB_set_ref = []
      dD_extField_set_ref = []
      
      dG_AA_set_ref = []
      dG_BB_set_ref = []
      dG_AB_set_ref = []
      dG_extField_set_ref = []
      
      S_AB_set = []
      H_set = []
      T_set = []
      V_set = []
      DIP_set = []      
      
      F_A_set = []
      F_B_set = []
      F_extField_set = []
      V_extField_set = []
      E_nuc_field_set= []
      
      W_AB_set= []
      W_BA_set= []
      w_AB_set= []
      w_BA_set= []
      
      A_AB_set= []
      A_BA_set= []
      a_AB_set= []
      a_BA_set= []
      
      K_set   = []
      L_set   = []
      k_set   = []
      l_set   = []

      #self._bfs_dimer_set = []

      psi4.core.print_out(" ---> Computing wavefunction for each sample <---\n\n")
      for s in range(self._nsamples):

          # Dimer system
          aggr, parcel = self._construct_aggregate()
          e_set.append(parcel[0])
          E_nuc_set.append(parcel[1])
          DIP_nuc_set.append(parcel[2])

          # Monomer with point charges
          if self._use_external_field_model:
             F_extField, V_extField, wfn_extField, E_nuc_field = self._extField()

          # Set sources of perturbation
          self._determine_perturbing_densities()

          # Overlap integrals
          S_AB = self._dimer_wfn.S().to_array(dense=True)[:self._nbf,self._nbf:]

          # Fock matrix
          F = self._dimer_wfn.Fa().to_array(dense=True)

          # Basis set of the dimer
          #self._bfs_dimer_set.append(self._dimer_wfn.basisset())

          # Core Hamiltonian
          H = self._dimer_wfn.H().to_array(dense=True)
          T = self._dimer_mints.ao_kinetic().to_array(dense=True)
          V = self._dimer_mints.ao_potential().to_array(dense=True)
          U = V.copy()
          mints = psi4.core.MintsHelper(self._bfs_0)
          V_0= mints.ao_potential()
          U[:self._nbf,:self._nbf] -= V_0
          U[self._nbf:,self._nbf:] -= V_0

          # OPDM (dimer)
          dD     = self._dimer_wfn.Da().to_array(dense=True)
          dD[:self._nbf,:self._nbf]-= self._D0
          dD[self._nbf:,self._nbf:]-= self._D0

          dD_AA = dD[:self._nbf,:self._nbf].copy()
          dD_BB = dD[self._nbf:,self._nbf:].copy()
          dD_AB = dD[:self._nbf,self._nbf:].copy()

          self._save_dD(dD, prefix='dd')

          # G tensor
          if   'e' in self._dms_types: # G = 2J - K
                G_ref   = F - H
          elif 'g' in self._dms_types: # G = 2J - K + U
                G_ref = F - H + U
          #      T_mon = T[:self._nbf,:self._nbf].copy()
          #      dG   = self._dimer_wfn.Fa().to_array(dense=True) - T
          #elif 'f' in self._dms_types: # G = T + V + 2J - K = F
          #      dG   = self._dimer_wfn.Fa().to_array(dense=True)
          #elif 'g1' in self._dms_types: # G = 2V + 2J - K #TODO
          #      dG   = self._dimer_wfn.Fa().to_array(dense=True) - T + V
          #      raise NotImplementedError # check it (but maybe not needed is this condition)

          dG = G_ref.copy()
          dG[:self._nbf,:self._nbf]-= self._G0
          dG[self._nbf:,self._nbf:]-= self._G0

          dG_AA = dG[:self._nbf,:self._nbf].copy()
          dG_BB = dG[self._nbf:,self._nbf:].copy()
          dG_AB = dG[:self._nbf,self._nbf:].copy()

          self._save_dD(dG, prefix='gg')

          if self._use_external_field_model:
             dD_extField = wfn_extField.Da().to_array(dense=True).copy()
             dD_extField-= self._D0

             dG_extField = wfn_extField.Fa().to_array(dense=True) - wfn_extField.H().to_array(dense=True) # 2J - K 
             dG_extField-= self._dms_extField_g._M.copy()
                                                                                                                   

          # Electric field
          F_A, F_B, F_A_mat, F_B_mat = self._compute_efield()

          # Dipole integrals
          DIP = [x.to_array(dense=True) for x in self._dimer_mints.ao_dipole()]
          DIP_set.append(DIP)

          # Auxiliary matrices
          W_AB = self._DA @ S_AB
          W_BA = self._DB @ S_AB.T
         #W_AB = S_AB.copy()
         #W_BA = S_AB.copy().T

         #ss = 1000.0
         #W_AB*= ss
         #W_BA*= ss

          A_AB = W_AB.T @ W_AB
          A_BA = W_BA.T @ W_BA

          K    = W_AB.T @ dD_AB
          L    = W_BA.T @ dD_AB.T

          w_AB = self._GA @ S_AB
          w_BA = self._GB @ S_AB.T

          a_AB = w_AB.T @ w_AB
          a_BA = w_BA.T @ w_BA

          k    = w_AB.T @ dG_AB
          l    = w_BA.T @ dG_AB.T

          # Accumulate
          S_AB_set.append(S_AB)
          H_set.append(H)
          T_set.append(T)
          V_set.append(V)

          dD_AA_set_ref.append(dD_AA)
          dD_BB_set_ref.append(dD_BB)
          dD_AB_set_ref.append(dD_AB)

          dG_AA_set_ref.append(dG_AA)
          dG_BB_set_ref.append(dG_BB)
          dG_AB_set_ref.append(dG_AB)

          F_A_set.append(F_A)
          F_B_set.append(F_B)

          if self._use_external_field_model:
             dD_extField_set_ref.append(dD_extField) 
             dG_extField_set_ref.append(dG_extField)
             F_extField_set.append(F_extField)
             V_extField_set.append(V_extField)
             E_nuc_field_set.append(E_nuc_field)

          W_AB_set.append(W_AB)
          W_BA_set.append(W_BA)

          A_AB_set.append(A_AB)
          A_BA_set.append(A_BA)

          K_set.append(K)
          L_set.append(L)

          w_AB_set.append(w_AB)
          w_BA_set.append(w_BA)

          a_AB_set.append(a_AB)
          a_BA_set.append(a_BA)

          k_set.append(k)
          l_set.append(l)
          #

      e_set = numpy.array(e_set)
      E_nuc_set = numpy.array(E_nuc_set)
      DIP_nuc_set = numpy.array(DIP_nuc_set)
      #
      S_AB_set= numpy.array(S_AB_set)
      H_set= numpy.array(H_set)
      T_set= numpy.array(T_set)
      V_set= numpy.array(V_set)
      DIP_set= numpy.array(DIP_set)
      #
      dD_AA_set_ref= numpy.array(dD_AA_set_ref)
      dD_BB_set_ref= numpy.array(dD_BB_set_ref)
      dD_AB_set_ref= numpy.array(dD_AB_set_ref)
      #
      dG_AA_set_ref= numpy.array(dG_AA_set_ref)
      dG_BB_set_ref= numpy.array(dG_BB_set_ref)
      dG_AB_set_ref= numpy.array(dG_AB_set_ref)
      #
      F_A_set = numpy.array(F_A_set)
      F_B_set = numpy.array(F_B_set)
      #
      if self._use_external_field_model:
         dD_extField_set_ref = numpy.array(dD_extField_set_ref) 
         dG_extField_set_ref = numpy.array(dG_extField_set_ref)
         F_extField_set = numpy.array(F_extField_set)
         V_extField_set = numpy.array(V_extField_set)
         E_nuc_field_set= numpy.array(E_nuc_field_set)
      #
      W_AB_set= numpy.array(W_AB_set)
      W_BA_set= numpy.array(W_BA_set)
      #
      A_AB_set= numpy.array(A_AB_set)
      A_BA_set= numpy.array(A_BA_set)
      #
      K_set   = numpy.array(K_set)
      L_set   = numpy.array(L_set)
      #
      w_AB_set= numpy.array(w_AB_set)
      w_BA_set= numpy.array(w_BA_set)
      #
      a_AB_set= numpy.array(a_AB_set)
      a_BA_set= numpy.array(a_BA_set)
      #
      k_set   = numpy.array(k_set)
      l_set   = numpy.array(l_set)

      # Average electric field in the sets
      fx = numpy.abs(F_A_set[:,:,0]).mean()
      fy = numpy.abs(F_A_set[:,:,1]).mean()
      fz = numpy.abs(F_A_set[:,:,2]).mean()
      psi4.core.print_out(" DMSFit: Average E-Field in clusters on A: FA_X FA_Y FA_Z = [%14.6f %14.6f %14.6f] A.U.\n" % (fx,fy,fz))

      fx = numpy.abs(F_B_set[:,:,0]).mean()
      fy = numpy.abs(F_B_set[:,:,1]).mean()
      fz = numpy.abs(F_B_set[:,:,2]).mean()
      psi4.core.print_out(" DMSFit: Average E-Field in clusters on B: FB_X FB_Y FB_Z = [%14.6f %14.6f %14.6f] A.U.\n" % (fx,fy,fz))

      # Save on disk
      e_set.tofile('temp_e_set.dat')
      E_nuc_set.tofile('temp_E_nuc_set.dat')
      DIP_nuc_set.tofile('temp_DIP_nuc_set.dat')
      #
      S_AB_set.tofile('temp_S_AB_set.dat')
      H_set.tofile('temp_H_set.dat')
      T_set.tofile('temp_T_set.dat')
      V_set.tofile('temp_V_set.dat')
      DIP_set.tofile('temp_DIP_set.dat')
      #
      dD_AA_set_ref.tofile('temp_dD_AA_ref_set.dat')
      dD_BB_set_ref.tofile('temp_dD_BB_ref_set.dat')
      dD_AB_set_ref.tofile('temp_dD_AB_ref_set.dat')
      #
      dG_AA_set_ref.tofile('temp_dG_AA_ref_set.dat')
      dG_BB_set_ref.tofile('temp_dG_BB_ref_set.dat')
      dG_AB_set_ref.tofile('temp_dG_AB_ref_set.dat')
      #
      F_A_set .tofile('temp_F_A_set.dat')
      F_B_set .tofile('temp_F_B_set.dat')
      #
      if self._use_external_field_model:
         dD_extField_set_ref.tofile('temp_dD_extField_ref_set.dat')
         dG_extField_set_ref.tofile('temp_dG_extField_ref_set.dat')
         F_extField_set .tofile('temp_F_extField_set.dat')
         V_extField_set .tofile('temp_V_extField_set.dat')
         E_nuc_field_set.tofile('temp_E_nuc_field_set.dat')
      #
      W_AB_set.tofile('temp_W_AB_set.dat')
      W_BA_set.tofile('temp_W_BA_set.dat')
      #
      A_AB_set.tofile('temp_A_AB_set.dat')
      A_BA_set.tofile('temp_A_BA_set.dat')
      #
      K_set   .tofile('temp_K_set.dat')
      L_set   .tofile('temp_L_set.dat')
      #
      w_AB_set.tofile('temp_w_AB_set.dat')
      w_BA_set.tofile('temp_w_BA_set.dat')
      #
      a_AB_set.tofile('temp_a_AB_set.dat')
      a_BA_set.tofile('temp_a_BA_set.dat')
      #
      k_set   .tofile('temp_k_set.dat')
      l_set   .tofile('temp_l_set.dat')
     

  def _compute_group_1(self, dms, dms_field,
                                  F_A_set, F_B_set, dM_AA_ref_set, dM_BB_ref_set, A_AB_set, A_BA_set,
                                  W_AB_set, W_BA_set, adjust=True):#TODO
      "Compute 1st group of parameters"
      psi4.core.print_out(" ---> Computing DMS for Z1 group of type %s <---\n\n" % dms._type_long)

      n = dms.n()
      N = dms.N()

      # ---> Fit 10 and 20 susceptibilities intrinsically <--- #
      if not self._use_external_field_model:

         # Hessian                                                                                                
         dim_1 = N*3
         dim_2 = N*6
         H = numpy.zeros((dim_1 + dim_2, dim_1 + dim_2))
                                                                                               
         H_10_10 = numpy.einsum("niu,njw->iujw", F_B_set, F_B_set) 
         H_10_10+= numpy.einsum("niu,njw->iujw", F_A_set, F_A_set) 
                                                                                               
         H[:dim_1,:dim_1] = H_10_10.copy().reshape(dim_1, dim_1)
         del H_10_10
                                                                                               
         H_20_20 = numpy.einsum("niu,nix,njw,njy->iuxjwy", F_B_set, F_B_set, F_B_set, F_B_set)
         H_20_20+= numpy.einsum("niu,nix,njw,njy->iuxjwy", F_A_set, F_A_set, F_A_set, F_A_set)
                                                                                               
         u = numpy.zeros((N,6,N,6))
         u[:,0,:,0] = H_20_20[:,0,0,:,0,0].copy() * 1.0
         u[:,0,:,1] = H_20_20[:,0,0,:,0,1].copy() * 2.0
         u[:,0,:,2] = H_20_20[:,0,0,:,0,2].copy() * 2.0
         u[:,0,:,3] = H_20_20[:,0,0,:,1,1].copy() * 1.0
         u[:,0,:,4] = H_20_20[:,0,0,:,1,2].copy() * 2.0
         u[:,0,:,5] = H_20_20[:,0,0,:,2,2].copy() * 1.0
                                                                                               
         u[:,1,:,0] = H_20_20[:,0,1,:,0,0].copy() * 2.0
         u[:,1,:,1] = H_20_20[:,0,1,:,0,1].copy() * 4.0
         u[:,1,:,2] = H_20_20[:,0,1,:,0,2].copy() * 4.0
         u[:,1,:,3] = H_20_20[:,0,1,:,1,1].copy() * 2.0
         u[:,1,:,4] = H_20_20[:,0,1,:,1,2].copy() * 4.0
         u[:,1,:,5] = H_20_20[:,0,1,:,2,2].copy() * 2.0
                                                                                               
         u[:,2,:,0] = H_20_20[:,0,2,:,0,0].copy() * 2.0
         u[:,2,:,1] = H_20_20[:,0,2,:,0,1].copy() * 4.0
         u[:,2,:,2] = H_20_20[:,0,2,:,0,2].copy() * 4.0
         u[:,2,:,3] = H_20_20[:,0,2,:,1,1].copy() * 2.0
         u[:,2,:,4] = H_20_20[:,0,2,:,1,2].copy() * 4.0
         u[:,2,:,5] = H_20_20[:,0,2,:,2,2].copy() * 2.0
                                                                                               
         u[:,3,:,0] = H_20_20[:,1,1,:,0,0].copy() * 1.0
         u[:,3,:,1] = H_20_20[:,1,1,:,0,1].copy() * 2.0
         u[:,3,:,2] = H_20_20[:,1,1,:,0,2].copy() * 2.0
         u[:,3,:,3] = H_20_20[:,1,1,:,1,1].copy() * 1.0
         u[:,3,:,4] = H_20_20[:,1,1,:,1,2].copy() * 2.0
         u[:,3,:,5] = H_20_20[:,1,1,:,2,2].copy() * 1.0
                                                                                               
         u[:,4,:,0] = H_20_20[:,1,2,:,0,0].copy() * 2.0
         u[:,4,:,1] = H_20_20[:,1,2,:,0,1].copy() * 4.0
         u[:,4,:,2] = H_20_20[:,1,2,:,0,2].copy() * 4.0
         u[:,4,:,3] = H_20_20[:,1,2,:,1,1].copy() * 2.0
         u[:,4,:,4] = H_20_20[:,1,2,:,1,2].copy() * 4.0
         u[:,4,:,5] = H_20_20[:,1,2,:,2,2].copy() * 2.0
                                                                                               
         u[:,5,:,0] = H_20_20[:,2,2,:,0,0].copy() * 1.0
         u[:,5,:,1] = H_20_20[:,2,2,:,0,1].copy() * 2.0
         u[:,5,:,2] = H_20_20[:,2,2,:,0,2].copy() * 2.0
         u[:,5,:,3] = H_20_20[:,2,2,:,1,1].copy() * 1.0
         u[:,5,:,4] = H_20_20[:,2,2,:,1,2].copy() * 2.0
         u[:,5,:,5] = H_20_20[:,2,2,:,2,2].copy() * 1.0
                                                                                               
         H[dim_1:,dim_1:] = u.copy().reshape(dim_2, dim_2)
         del u, H_20_20
                                                                                               
         H_10_20 = numpy.einsum("niu,njw,njy->iujwy", F_B_set, F_B_set, F_B_set)
         H_10_20+= numpy.einsum("niu,njw,njy->iujwy", F_A_set, F_A_set, F_A_set)
                                                                                               
         u = numpy.zeros((N,3,N,6))
         u[:,:,:,0] = H_10_20[:,:,:,0,0].copy() * 1.0
         u[:,:,:,1] = H_10_20[:,:,:,0,1].copy() * 2.0
         u[:,:,:,2] = H_10_20[:,:,:,0,2].copy() * 2.0
         u[:,:,:,3] = H_10_20[:,:,:,1,1].copy() * 1.0
         u[:,:,:,4] = H_10_20[:,:,:,1,2].copy() * 2.0
         u[:,:,:,5] = H_10_20[:,:,:,2,2].copy() * 1.0
                                                                                               
         H[:dim_1,dim_1:] = u.copy().reshape(dim_1, dim_2)
         H[dim_1:,:dim_1] = H[:dim_1,dim_1:].copy().T
         del u, H_10_20
                                                                                               
         H *= 2.0
                                                                                                                  
                                                                                                                  
         # Only (1,0) and (2,0) susceptibilities: fitting for each target matrix element separately
         if (0,2) not in self._dms_da.available_orders():
                                                                                                                  
             # Fit
             psi4.core.print_out(" * Computing Hessian Inverse...\n")
             Hi = self._invert_hessian(H)
                                                                                                   
             S_1 = numpy.zeros((n,n,N,3))
             S_2 = numpy.zeros((n,n,N,6))
                                                                                                   
             for i in range(n):
                 for j in range(n):
                    #psi4.core.print_out(" * Computing the DMS tensors for (%i %i)...\n" % (i,j))
                                                                                                   
                     # gradient                                                           
                     DIM_1 = N*3
                     DIM_2 = N*6
                                                                                          
                     g_10 = numpy.einsum("n,niu->iu", dM_AA_ref_set[:,i,j], F_B_set)
                     g_10+= numpy.einsum("n,niu->iu", dM_BB_ref_set[:,i,j], F_A_set)
                                                                                          
                     u = numpy.zeros((N,6))
                     g_20 = numpy.einsum("n,niu,niw->iuw", dM_AA_ref_set[:,i,j], F_B_set, F_B_set)
                     g_20+= numpy.einsum("n,niu,niw->iuw", dM_BB_ref_set[:,i,j], F_A_set, F_A_set)
                     u[:,0] = g_20[:,0,0].copy()
                     u[:,1] = g_20[:,0,1].copy() * 2.0
                     u[:,2] = g_20[:,0,2].copy() * 2.0
                     u[:,3] = g_20[:,1,1].copy()
                     u[:,4] = g_20[:,1,2].copy() * 2.0
                     u[:,5] = g_20[:,2,2].copy()
                                                                                          
                     g = numpy.zeros(dim_1 + dim_2)
                                                                                          
                     g[:dim_1] = g_10.ravel().copy()
                     g[dim_1:] = u.ravel()
                     del u, g_10, g_20
                     g *= -2.0
                                                                                          
                     s = - g @ Hi
                                                                                                   
                     S_1[i,j] = s[:dim_1].copy().reshape(N,3)
                     S_2[i,j] = s[dim_1:].copy().reshape(N,6)
                                                                                                                  
             S_1 = composite.form_futf(S_1)
             S_2 = composite.form_futf(S_2)
                                                                                                   
             s = numpy.hstack([S_1.ravel(), S_2.ravel()])
                                                                                                                  
         # (0,2) and higher order susceptibilities present
         else:
            h_10_10 = H[:dim_1,:dim_1].copy() / 2.0
            h_10_20 = H[:dim_1,dim_1:].copy() / 2.0
            h_20_20 = H[dim_1:,dim_1:].copy() / 2.0
                                                                                                                  
            FF_A_set = numpy.einsum("niu,niw->uwni", F_A_set, F_A_set)
            FF_B_set = numpy.einsum("niu,niw->uwni", F_B_set, F_B_set)
            FF_A_set = numpy.einsum("uw,uwni->uwni", composite.symmetry_matrix(3), FF_A_set)
            FF_B_set = numpy.einsum("uw,uwni->uwni", composite.symmetry_matrix(3), FF_B_set)
            FF_A_set = composite.form_futf(FF_A_set).transpose(1,2,0)
            FF_B_set = composite.form_futf(FF_B_set).transpose(1,2,0)
                                                                                                                  
            I = numpy.identity(n)
            r = composite.symmetry_matrix(n)
            t = numpy.triu(numpy.ones(n))
                                                                                                                  
            del H # --> move it above!
                                                                                                                  
            n_ = composite.number_of_elements(n)
            DIM_1 = n_ * dim_1
            DIM_2 = n_ * dim_2
            DIM_3 = n*n
            DIM = DIM_1 + DIM_2 + DIM_3
            g = numpy.zeros(DIM)
            H = numpy.zeros((DIM, DIM))
            O1 = DIM_1 + DIM_2
            O2 = O1 + DIM_3
                                                                                                                  
            # Construct gradient
            # (10)
            g_10 = numpy.zeros((n,n,N,3))
            for s in range(self._nsamples):
                dM_AA_ref = dM_AA_ref_set[s]
                dM_BB_ref = dM_BB_ref_set[s]
                F_A = F_A_set[s]
                F_B = F_B_set[s]
                X_AA = numpy.triu(dM_AA_ref) #composite.partial_contraction_with_trace(I, dM_AA_ref)
                X_BB = numpy.triu(dM_BB_ref) #composite.partial_contraction_with_trace(I, dM_BB_ref)
                g_10+= numpy.einsum("ab,iu->abiu", X_AA, F_B)
                g_10+= numpy.einsum("ab,iu->abiu", X_BB, F_A)
            g[:DIM_1] = composite.form_futf(g_10).reshape(DIM_1).copy()
           #g_10 = numpy.einsum("nab,niu->abiu", dM_AA_ref_set, F_B_set)
           #g_10+= numpy.einsum("nab,niu->abiu", dM_BB_ref_set, F_A_set)
            del g_10
                                                                                                                  
            # (20)
            g_20 = numpy.zeros((n,n,N,6))
            for s in range(self._nsamples):
                dM_AA_ref = dM_AA_ref_set[s]
                dM_BB_ref = dM_BB_ref_set[s]
                FF_A = FF_A_set[s]
                FF_B = FF_B_set[s]
                X_AA = numpy.triu(dM_AA_ref) #composite.partial_contraction_with_trace(I, dM_AA_ref)
                X_BB = numpy.triu(dM_BB_ref) #composite.partial_contraction_with_trace(I, dM_BB_ref)
                g_20+= numpy.einsum("ab,iu->abiu", X_AA, FF_B)
                g_20+= numpy.einsum("ab,iu->abiu", X_BB, FF_A)
            g[DIM_1:O1] = composite.form_futf(g_20).reshape(DIM_2).copy()
            del g_20
                                                                                                                  
            # (02)
            S_AB_set = numpy.fromfile("temp_S_AB_set.dat").reshape(self._nsamples,n,n) 
            g_02 = numpy.zeros((n,n))
            tau = numpy.ones((n,n)) - I
            t   = numpy.triu(numpy.ones(n))
            for s in range(self._nsamples):
                W_AB = W_AB_set[s]
                W_BA = W_BA_set[s]
                S_AB = S_AB_set[s]
                S_BA = S_AB.T.copy()
                W_ABA= W_AB @ S_BA
                W_BAB= W_BA @ S_AB
                dM_AA_ref = dM_AA_ref_set[s]
                dM_BB_ref = dM_BB_ref_set[s]
                X_ABA= composite.partial_contraction_with_trace(W_ABA, dM_AA_ref)
                X_BAB= composite.partial_contraction_with_trace(W_BAB, dM_BB_ref)
                Y_ABA= composite.partial_contraction_with_trace(I, dM_AA_ref)
                Y_BAB= composite.partial_contraction_with_trace(I, dM_BB_ref)
                g_02 += W_ABA.T @ Y_ABA + X_ABA.T
                g_02 += W_BAB.T @ Y_BAB + X_BAB.T

            g[O1:O2] = g_02.reshape(DIM_3).copy()
            del g_02
                                                                                                                  
            g *= -2.0
                                                                                                                  
            # Construct Hessian
            H_4 = numpy.einsum("ac,abd->abcd",I,composite.partial_contraction(I,I))
                                                                                                                  
            # (10,10)
            H_10_10 = numpy.einsum("abcd,iujw->abiucdjw", H_4, h_10_10.reshape(N,3,N,3))
           #H_10_10 = numpy.einsum("ac,bd,iujw->abiucdjw",I,I,h_10_10.reshape(N,3,N,3))
            H[     :DIM_1,     :DIM_1] = composite.form_superfutf(H_10_10, m1=4).reshape(DIM_1, DIM_1).copy()
            del H_10_10
                                                                                                                  
            # (20,20)
            H_20_20 = numpy.einsum("abcd,iujw->abiucdjw", H_4, h_20_20.reshape(N,6,N,6))
           #H_20_20 = numpy.einsum("ac,bd,iujw->abiucdjw",I,I,h_20_20.reshape(N,6,N,6))
            H[DIM_1:O1   ,DIM_1:O1   ] = composite.form_superfutf(H_20_20, m1=4).reshape(DIM_2, DIM_2).copy()
            del H_20_20
                                                                                                                  
            # (10,20)
            H_10_20 = numpy.einsum("abcd,iujw->abiucdjw", H_4, h_10_20.reshape(N,3,N,6))
           #H_10_20 = numpy.einsum("ac,bd,iujw->abiucdjw",I,I,h_10_20.reshape(N,3,N,6))
            H[:DIM_1, DIM_1:O1] = composite.form_superfutf(H_10_20, m1=4).reshape(DIM_1, DIM_2).copy()
            del H_10_20, H_4
            H[DIM_1:O1, :DIM_1] = H[:DIM_1, DIM_1:O1].T.copy()
                                                                                                                  
            # (02,02)
            H_02_02 = numpy.zeros((n,n,n,n))
            dd = composite.partial_contraction(I,I)
            for s in range(self._nsamples):
                W_AB = W_AB_set[s]
                W_BA = W_BA_set[s]
                S_AB = S_AB_set[s]
                S_BA = S_AB.T.copy()
                W_ABA= W_AB @ S_BA
                W_BAB= W_BA @ S_AB

                X_ABA = composite.partial_contraction(W_ABA, W_ABA)
                X_BAB = composite.partial_contraction(W_BAB, W_BAB)
                Z_ABA = composite.partial_contraction(W_ABA, I    )
                Z_BAB = composite.partial_contraction(W_BAB, I    )

                H_02_02 += numpy.einsum("pa,pc,pbd->abcd", W_ABA, W_ABA, dd)
                H_02_02 += numpy.einsum("dac,bd->abcd", X_ABA, I)
                H_02_02 += numpy.einsum("bc,bad->abcd", W_ABA, Z_ABA)
                H_02_02 += numpy.einsum("da,dcb->abcd", W_ABA, Z_ABA)
                #
                H_02_02 += numpy.einsum("pa,pc,pbd->abcd", W_BAB, W_BAB, dd)
                H_02_02 += numpy.einsum("dac,bd->abcd", X_BAB, I)
                H_02_02 += numpy.einsum("bc,bad->abcd", W_BAB, Z_BAB)
                H_02_02 += numpy.einsum("da,dcb->abcd", W_BAB, Z_BAB)

            H[O1:O2   ,O1:O2   ] = H_02_02.reshape(DIM_3, DIM_3).copy()
            del H_02_02
                                                                                                                  
            # (10,02)
            H_10_02 = numpy.zeros((n,n,N*3,n,n))
            for s in range(self._nsamples):
                W_AB = W_AB_set[s] ; F_A = F_A_set[s].reshape(3*N)
                W_BA = W_BA_set[s] ; F_B = F_B_set[s].reshape(3*N)
                S_AB = S_AB_set[s]
                S_BA = S_AB.T.copy()

                W_ABA= W_AB @ S_BA
                W_BAB= W_BA @ S_AB

                Z_ABA = composite.partial_contraction(W_ABA, I    )
                Z_BAB = composite.partial_contraction(W_BAB, I    )

                H_10_02 += numpy.einsum("ac,abd,q->abqcd", W_ABA, dd, F_B)
                H_10_02 += numpy.einsum("ad,acb,q->abqcd", I, Z_ABA, F_B)
                #
                H_10_02 += numpy.einsum("ac,abd,q->abqcd", W_BAB, dd, F_A)
                H_10_02 += numpy.einsum("ad,acb,q->abqcd", I, Z_BAB, F_A)

            H[:DIM_1,O1:O2] = composite.form_futf(H_10_02).reshape(DIM_1, DIM_3).copy()
            del H_10_02
            H[O1:O2,:DIM_1] = H[:DIM_1,O1:O2].T.copy()
                                                                                                                  
            # (20,02)
            H_20_02 = numpy.zeros((n,n,N*6,n,n))
            for s in range(self._nsamples):
                W_AB = W_AB_set[s] ; FF_A = FF_A_set[s].reshape(6*N)
                W_BA = W_BA_set[s] ; FF_B = FF_B_set[s].reshape(6*N)
                S_AB = S_AB_set[s]
                S_BA = S_AB.T.copy()

                W_ABA= W_AB @ S_BA
                W_BAB= W_BA @ S_AB

                Z_ABA = composite.partial_contraction(W_ABA, I    )
                Z_BAB = composite.partial_contraction(W_BAB, I    )

                H_20_02 += numpy.einsum("ac,abd,q->abqcd", W_ABA, dd, FF_B)
                H_20_02 += numpy.einsum("ad,acb,q->abqcd", I, Z_ABA, FF_B)
                #
                H_20_02 += numpy.einsum("ac,abd,q->abqcd", W_BAB, dd, FF_A)
                H_20_02 += numpy.einsum("ad,acb,q->abqcd", I, Z_BAB, FF_A)

            H[DIM_1:O1,O1:O2] = composite.form_futf(H_20_02).reshape(DIM_2, DIM_3).copy()
            del H_20_02
            H[O1:O2,DIM_1:O1] = H[DIM_1:O1,O1:O2].T.copy()
                                                                                                                  
            H *= 2.0
            Hi = self._invert_hessian(H)
                                                                                                                  
            # Fit
            s = - g @ Hi
                                                                                                                  
                                                                                                                  
         # Extract susceptibilities
         psi4.core.print_out(" * Setting the DMS tensors...\n")
         dms.set_s1(s)
         for order in [(1,0),(2,0),(0,2)]:
             if order in dms.available_orders(): dms._generate_B_from_s(order)

      # ---> External Field Model <--- #
      else:
         n = dms.n()
         N = dms.N()
         s = self._nsamples
         #
         OFFs= [0,]
         DIM = 0
         FIELD_ONLY = True
         # (0,2)
         if (0,2) in dms.available_orders():
             DIM_1 = n*n
             DIM+= DIM_1
             OFFs.append(DIM)
             FIELD_ONLY = False
         # (0,4)
         if (0,4) in dms.available_orders():
             DIM_2 = n*n
             DIM+= DIM_2
             OFFs.append(DIM)
             FIELD_ONLY = False
         # (0,6)
         if (0,6) in dms.available_orders():
             DIM_3 = n*n
             DIM+= DIM_3
             OFFs.append(DIM)
             FIELD_ONLY = False

         # allocate 
         g = numpy.zeros(DIM)
         H = numpy.zeros((DIM, DIM))
                                                                                  
         # gradient
         if not FIELD_ONLY:
            S_AB_set = numpy.fromfile("temp_S_AB_set.dat").reshape(s,n,n)
            S_BA_set = S_AB_set.transpose(0,2,1)
            W_ABA_set = numpy.einsum("nab,nbc->nac", W_AB_set, S_BA_set)
            W_BAB_set = numpy.einsum("nab,nbc->nac", W_BA_set, S_AB_set)
            W_ABABA_set = numpy.einsum("nab,nbd,ndc->nac", W_ABA_set, S_AB_set, S_BA_set)
            W_BABAB_set = numpy.einsum("nab,nbd,ndc->nac", W_BAB_set, S_BA_set, S_AB_set)
            W_ABABABA_set = numpy.einsum("nab,nbd,ndc->nac", W_ABABA_set, S_AB_set, S_BA_set)
            W_BABABAB_set = numpy.einsum("nab,nbd,ndc->nac", W_BABAB_set, S_BA_set, S_AB_set)

            B_10_field= dms_field.B(1,0)
            B_20_field= dms_field.B(2,0)
            F_A_set_2 = numpy.einsum("niu,niw->niuw", F_A_set, F_A_set)
            F_B_set_2 = numpy.einsum("niu,niw->niuw", F_B_set, F_B_set)
            FF_A_set = composite.form_futf(numpy.einsum("uw,niuw->uwni", composite.symmetry_matrix(3), F_A_set_2)).transpose(1,2,0)
            FF_B_set = composite.form_futf(numpy.einsum("uw,niuw->uwni", composite.symmetry_matrix(3), F_B_set_2)).transpose(1,2,0)

            dM_AA_start_set = numpy.einsum("abiu,niu->nab", B_10_field, F_B_set)
            dM_AA_start_set+= numpy.einsum("abiuw,niuw->nab", B_20_field, F_B_set_2)
            dM_BB_start_set = numpy.einsum("abiu,niu->nab", B_10_field, F_A_set)
            dM_BB_start_set+= numpy.einsum("abiuw,niuw->nab", B_20_field, F_A_set_2)

            ddM_AA_set = dM_AA_ref_set.copy()
            ddM_BB_set = dM_BB_ref_set.copy()
            if adjust:
               ddM_AA_set -= dM_AA_start_set
               ddM_BB_set -= dM_BB_start_set

            # (02) block
            psi4.core.print_out(" * Computing Gradient (0,2)...\n")
            START = OFFs[0]
            END   = OFFs[1]

            g_02 = 2.0 * numpy.einsum("nad,nac->cd", ddM_AA_set, W_ABA_set)
            g_02+= 2.0 * numpy.einsum("nad,nac->cd", ddM_BB_set, W_BAB_set)

            g[START:END] = g_02.ravel().copy(); del g_02
 
            # (04) block
            if (0,4) in dms.available_orders():
                psi4.core.print_out(" * Computing Gradient (0,4)...\n")
                START = OFFs[1]
                END   = OFFs[2]
                #
                g_04 = 2.0 * numpy.einsum("nad,nac->cd", ddM_AA_set, W_ABABA_set)
                g_04+= 2.0 * numpy.einsum("nad,nac->cd", ddM_BB_set, W_BABAB_set)

                g[START:END] = g_04.reshape(DIM_2).copy(); del g_04

            # (06) block
            if (0,6) in dms.available_orders():
                psi4.core.print_out(" * Computing Gradient (0,6)...\n")
                START = OFFs[2]
                END   = OFFs[3]
                #
                g_06 = 2.0 * numpy.einsum("nad,nac->cd", ddM_AA_set, W_ABABABA_set)
                g_06+= 2.0 * numpy.einsum("nad,nac->cd", ddM_BB_set, W_BABABAB_set)

                g[START:END] = g_06.reshape(DIM_3).copy(); del g_06

            g *= -2.0

            # Hessian                                                         
            I = numpy.identity(n)

            # (02) susceptibility
            START = OFFs[0]
            END   = OFFs[1]
            #
            L = END - START
            #
            psi4.core.print_out(" * Computing Hessian (0,2) (0,2)...\n")
            W_ABA_ABA = numpy.einsum("nab,nac->nbc", W_ABA_set, W_ABA_set)
            W_BAB_BAB = numpy.einsum("nab,nac->nbc", W_BAB_set, W_BAB_set)
            A_ABA= W_ABA_ABA.sum(axis=0)
            A_BAB= W_BAB_BAB.sum(axis=0)

            H_02_02 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABA    , I        )
            H_02_02+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABA_set, W_ABA_set)
            H_02_02+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BAB    , I        )
            H_02_02+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BAB_set, W_BAB_set)
            #
            H[START:END,START:END] = H_02_02.reshape(L,L).copy()
            del H_02_02

            # (04) susceptibility
            if (0,4) in dms.available_orders():
                PREV  = OFFs[0]
                START = OFFs[1]
                END   = OFFs[2]
                #
                # (0,4) (0,4) block
                psi4.core.print_out(" * Computing Hessian (0,4) (0,4)...\n")
                W_ABABA_ABABA = numpy.einsum("nab,nac->nbc", W_ABABA_set, W_ABABA_set)
                W_BABAB_BABAB = numpy.einsum("nab,nac->nbc", W_BABAB_set, W_BABAB_set)
                A_ABABA= W_ABABA_ABABA.sum(axis=0)
                A_BABAB= W_BABAB_BABAB.sum(axis=0)

                H_04_04 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABABA    , I        )
                H_04_04+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABABA_set, W_ABABA_set)
                H_04_04+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BABAB    , I        )
                H_04_04+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BABAB_set, W_BABAB_set)
                #
                H[START:END,START:END] = H_04_04.reshape(DIM_2,DIM_2).copy(); del H_04_04

                # (0,2) (0,4)
                psi4.core.print_out(" * Computing Hessian (0,4) (0,4)...\n")
                W_ABA_ABABA = numpy.einsum("nab,nac->nbc", W_ABA_set, W_ABABA_set)
                W_BAB_BABAB = numpy.einsum("nab,nac->nbc", W_BAB_set, W_BABAB_set)
                A_ABA= W_ABA_ABABA.sum(axis=0)
                A_BAB= W_BAB_BABAB.sum(axis=0)

                H_02_04 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABA    , I        )
                H_02_04+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABA_set, W_ABABA_set)
                H_02_04+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BAB    , I        )
                H_02_04+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BAB_set, W_BABAB_set)
                #
                H[PREV:START,START:END] = H_02_04.reshape(DIM_1,DIM_2).copy(); del H_02_04
                H[START:END,PREV:START] = H[PREV:START,START:END].T.copy()

            # (06) susceptibility
            if (0,6) in dms.available_orders():
                DPREV = OFFs[0]
                PREV  = OFFs[1]
                START = OFFs[2]
                END   = OFFs[3]
                #
                # (0,6) (0,6) block
                psi4.core.print_out(" * Computing Hessian (0,6) (0,6)...\n")
                W_ABABABA_ABABABA = numpy.einsum("nab,nac->nbc", W_ABABABA_set, W_ABABABA_set)
                W_BABABAB_BABABAB = numpy.einsum("nab,nac->nbc", W_BABABAB_set, W_BABABAB_set)
                A_ABABABA= W_ABABABA_ABABABA.sum(axis=0)
                A_BABABAB= W_BABABAB_BABABAB.sum(axis=0)

                H_06_06 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABABABA    , I            )
                H_06_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABABABA_set, W_ABABABA_set)
                H_06_06+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BABABAB    , I            )
                H_06_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BABABAB_set, W_BABABAB_set)
                #
                H[START:END,START:END] = H_06_06.reshape(DIM_3,DIM_3).copy(); del H_06_06

                # (0,4) (0,6)
                psi4.core.print_out(" * Computing Hessian (0,4) (0,6)...\n")
                W_ABABA_ABABABA = numpy.einsum("nab,nac->nbc", W_ABABA_set, W_ABABABA_set)
                W_BABAB_BABABAB = numpy.einsum("nab,nac->nbc", W_BABAB_set, W_BABABAB_set)
                A_ABABA= W_ABABA_ABABABA.sum(axis=0)
                A_BABAB= W_BABAB_BABABAB.sum(axis=0)

                H_04_06 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABABA    , I            )
                H_04_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABABA_set, W_ABABABA_set)
                H_04_06+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BABAB    , I            )
                H_04_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BABAB_set, W_BABABAB_set)
                #
                H[PREV:START,START:END] = H_04_06.reshape(DIM_2,DIM_3).copy(); del H_04_06 
                H[START:END,PREV:START] = H[PREV:START,START:END].T.copy()

                # (0,2) (0,6)
                psi4.core.print_out(" * Computing Hessian (0,2) (0,6)...\n")
                W_ABA_ABABABA = numpy.einsum("nab,nac->nbc", W_ABA_set, W_ABABABA_set)
                W_BAB_BABABAB = numpy.einsum("nab,nac->nbc", W_BAB_set, W_BABABAB_set)
                A_ABA= W_ABA_ABABABA.sum(axis=0)
                A_BAB= W_BAB_BABABAB.sum(axis=0)

                H_02_06 = 2.0 * numpy.einsum("ac,bd->abcd",   A_ABA    , I            )
                H_02_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_ABA_set, W_ABABABA_set)
                H_02_06+= 2.0 * numpy.einsum("ac,bd->abcd",   A_BAB    , I            )
                H_02_06+= 2.0 * numpy.einsum("nda,nbc->abcd", W_BAB_set, W_BABABAB_set)
                #
                H[DPREV:PREV,START:END] = H_02_06.reshape(DIM_1,DIM_3).copy(); del H_02_06 
                H[START:END,DPREV:PREV] = H[DPREV:PREV,START:END].T.copy()

            H *= 2.0
            Hi = self._invert_hessian(H)
                                                                                                                  
            # Fit
            s = - g @ Hi
                                                                                                                                    
            # combine with external field DMS's
            s_field = dms_field._s1.copy()
            s = numpy.hstack([s_field, s])

         else: 
            s = dms_field._s1.copy()
                                                                                                                  
         # Extract susceptibilities
         psi4.core.print_out(" * Setting the DMS tensors...\n")
         dms.set_s1(s)
         for order in [(1,0),(2,0),(0,2),(0,4),(0,6)]:
             if order in dms.available_orders(): dms._generate_B_from_s(order)



  def _compute_group_2(self, dms, K_set, L_set, dM_AB_set, W_AB_set, W_BA_set, A_AB_set, A_BA_set):
      "Compute 2nd group of parameters"
      psi4.core.print_out(" ---> Computing DMS for Z2 group of type %s <---\n\n" % dms._type_long)
      n = dms.n()
      N = dms.N()
      s = self._nsamples
      #
      OFFs= [0,]
      # (0,1)
      DIM = n*n
      OFFs.append(DIM)
      # (0,3)
      if (0,3) in dms.available_orders():
          DIM += n*n
          OFFs.append(DIM)
      # (0,5)
      if (0,5) in dms.available_orders():
          DIM += n*n
          OFFs.append(DIM)
      DIM_1 = n*n
      DIM_2 = DIM_1
      DIM_3 = DIM_2
   
      # allocate 
      g = numpy.zeros(DIM)
      H = numpy.zeros((DIM, DIM))

      # gradient
      # (01) block
      psi4.core.print_out(" * Computing Gradient (0,1)...\n")
      START = OFFs[0]
      END   = OFFs[1]
      KL = K_set + L_set
      g[START:END] = KL.sum(axis=0).ravel()

      # (03) block
      if (0,3) in dms.available_orders():
          psi4.core.print_out(" * Computing Gradient (0,3)...\n")
          START = OFFs[1]
          END   = OFFs[2]
          #
          S_AB_set = numpy.fromfile("temp_S_AB_set.dat").reshape(s,n,n)
          S_BA_set = S_AB_set.transpose(0,2,1)
          W_ABAB_set = numpy.einsum("nab,nbc,ncd->nad", W_AB_set, S_BA_set, S_AB_set)
          W_BABA_set = numpy.einsum("nab,nbc,ncd->nad", W_BA_set, S_AB_set, S_BA_set)
          dM_BA_set = dM_AB_set.transpose(0,2,1)
          #
          K_set = numpy.einsum("nad,nab->ndb", W_ABAB_set, dM_AB_set)
          L_set = numpy.einsum("nad,nab->ndb", W_BABA_set, dM_BA_set)
          KL = K_set + L_set
          #
          g[START:END] = KL.sum(axis=0).ravel()

      # (05) block
      if (0,5) in dms.available_orders():
          psi4.core.print_out(" * Computing Gradient (0,5)...\n")
          START = OFFs[2]
          END   = OFFs[3]
          #
          S_AB_set = numpy.fromfile("temp_S_AB_set.dat").reshape(s,n,n)
          S_BA_set = S_AB_set.transpose(0,2,1)
          W_ABABAB_set = numpy.einsum("nab,nbc,ncd->nad", W_ABAB_set, S_BA_set, S_AB_set)
          W_BABABA_set = numpy.einsum("nab,nbc,ncd->nad", W_BABA_set, S_AB_set, S_BA_set)
          dM_BA_set = dM_AB_set.transpose(0,2,1)
          #
          K_set = numpy.einsum("nad,nab->ndb", W_ABABAB_set, dM_AB_set)
          L_set = numpy.einsum("nad,nab->ndb", W_BABABA_set, dM_BA_set)
          KL = K_set + L_set
          #
          g[START:END] = KL.sum(axis=0).ravel()


      g *= -2.0

      # Hessian
      I = numpy.identity(n)
      A = A_AB_set + A_BA_set

      # (01) susceptibility
      START = OFFs[0]
      END   = OFFs[1]
      #
      L = END - START
      #
      psi4.core.print_out(" * Computing Hessian (0,1) (0,1)...\n")
      H_01_01 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
      H_01_01+= numpy.einsum("nda,nbc->abcd", W_AB_set, W_BA_set)
      H_01_01+= numpy.einsum("nda,nbc->abcd", W_BA_set, W_AB_set)
      #
      H[START:END,START:END] = H_01_01.reshape(L,L).copy()
      del H_01_01

      # (03) susceptibility
      if (0,3) in dms.available_orders():
          PREV  = OFFs[0]
          START = OFFs[1]
          END   = OFFs[2]
          #
          psi4.core.print_out(" * Computing Hessian (0,3) (0,3)...\n")
          #
          A_AB_set = numpy.einsum("nab,nac->nbc", W_ABAB_set, W_ABAB_set)
          A_BA_set = numpy.einsum("nab,nac->nbc", W_BABA_set, W_BABA_set)
          A = A_AB_set + A_BA_set
          #
          H_03_03 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
          H_03_03+= numpy.einsum("nda,nbc->abcd", W_ABAB_set, W_BABA_set)
          H_03_03+= numpy.einsum("nda,nbc->abcd", W_BABA_set, W_ABAB_set)
          #
          H[START:END,START:END] = H_03_03.reshape(DIM_2,DIM_2).copy(); del H_03_03

          psi4.core.print_out(" * Computing Hessian (0,1) (0,3)...\n")
          #
          A_AB_set = numpy.einsum("nab,nac->nbc", W_AB_set, W_ABAB_set)
          A_BA_set = numpy.einsum("nab,nac->nbc", W_BA_set, W_BABA_set)
          A = A_AB_set + A_BA_set
          #
          H_01_03 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
          H_01_03+= numpy.einsum("nda,nbc->abcd", W_AB_set, W_BABA_set)
          H_01_03+= numpy.einsum("nda,nbc->abcd", W_BA_set, W_ABAB_set)
          #
          H[PREV:START,START:END] = H_01_03.reshape(DIM_1,DIM_2).copy(); del H_01_03
          H[START:END,PREV:START] = H[PREV:START,START:END].T.copy()

      # (05) susceptibility
      if (0,5) in dms.available_orders():
          DPREV = OFFs[0]
          PREV  = OFFs[1]
          START = OFFs[2]
          END   = OFFs[3]
          #
          psi4.core.print_out(" * Computing Hessian (0,5) (0,5)...\n")
          #
          A_AB_set = numpy.einsum("nab,nac->nbc", W_ABABAB_set, W_ABABAB_set)
          A_BA_set = numpy.einsum("nab,nac->nbc", W_BABABA_set, W_BABABA_set)
          A = A_AB_set + A_BA_set
          #
          H_05_05 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
          H_05_05+= numpy.einsum("nda,nbc->abcd", W_ABABAB_set, W_BABABA_set)
          H_05_05+= numpy.einsum("nda,nbc->abcd", W_BABABA_set, W_ABABAB_set)
          #
          H[START:END,START:END] = H_05_05.reshape(DIM_3,DIM_3).copy(); del H_05_05

          psi4.core.print_out(" * Computing Hessian (0,3) (0,5)...\n")
          #
          A_AB_set = numpy.einsum("nab,nac->nbc", W_ABAB_set, W_ABABAB_set)
          A_BA_set = numpy.einsum("nab,nac->nbc", W_BABA_set, W_BABABA_set)
          A = A_AB_set + A_BA_set
          #
          H_03_05 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
          H_03_05+= numpy.einsum("nda,nbc->abcd", W_ABAB_set, W_BABABA_set)
          H_03_05+= numpy.einsum("nda,nbc->abcd", W_BABA_set, W_ABABAB_set)
          #
          H[PREV:START,START:END] = H_03_05.reshape(DIM_2,DIM_3).copy(); del H_03_05
          H[START:END,PREV:START] = H[PREV:START,START:END].T.copy()

          psi4.core.print_out(" * Computing Hessian (0,1) (0,5)...\n")
          #
          A_AB_set = numpy.einsum("nab,nac->nbc", W_AB_set, W_ABABAB_set)
          A_BA_set = numpy.einsum("nab,nac->nbc", W_BA_set, W_BABABA_set)
          A = A_AB_set + A_BA_set
          #
          H_01_05 = numpy.einsum("bd,ac->abcd", I, A.sum(axis=0)); del A
          H_01_05+= numpy.einsum("nda,nbc->abcd", W_AB_set, W_BABABA_set)
          H_01_05+= numpy.einsum("nda,nbc->abcd", W_BA_set, W_ABABAB_set)
          #
          H[DPREV:PREV,START:END] = H_01_05.reshape(DIM_1,DIM_3).copy(); del H_01_05
          H[START:END,DPREV:PREV] = H[DPREV:PREV,START:END].T.copy()

      H *= 2.0

      # Fit
      psi4.core.print_out(" * Computing Hessian Inverse...\n")
      Hi = self._invert_hessian(H)

      psi4.core.print_out(" * Computing the DMS tensors...\n")
      s = - g @ Hi

      # Extract susceptibilities
      psi4.core.print_out(" * Setting the DMS tensors...\n")
      dms.set_s2(s)
      for order in [(0,1),(0,3),(0,5)]:
          if order in dms.available_orders(): dms._generate_B_from_s(order)


  def _check(self):
      "Check the quality of fitting on training set"
      s = self._nsamples
      n = self._dms_da.n()
      N = self._dms_da.N()

      t = numpy.triu(numpy.ones(n))
      r = composite.symmetry_matrix(n)

      # compulsory DMS tensors
      B_10 = self.B(1,0,'da')
      B_20 = self.B(2,0,'da')
      b_10 = self.B(1,0,'g')
      b_20 = self.B(2,0,'g')

      if DMSFit.compute_dmatpol_susceptibilities:
         B_10_extField_dmatpol = self.B(1,0,'dmatpol')
         B_20_extField_dmatpol = self.B(2,0,'dmatpol')


      # additional DMS tensors
      if (0,1) in self._dms_da.available_orders():
          B_01 = self.B(0,1,'da')
          b_01 = self.B(0,1,'g')
      if (0,2) in self._dms_da.available_orders():
          B_02 = self.B(0,2,'da')
          b_02 = self.B(0,2,'g')
      if (0,3) in self._dms_da.available_orders():
          B_03 = self.B(0,3,'da')
          b_03 = self.B(0,3,'g')
      if (0,4) in self._dms_da.available_orders():
          B_04 = self.B(0,4,'da')
          b_04 = self.B(0,4,'g')
      if (0,5) in self._dms_da.available_orders():
          B_05 = self.B(0,5,'da')
          b_05 = self.B(0,5,'g')
      if (0,6) in self._dms_da.available_orders():
          B_06 = self.B(0,6,'da')
          b_06 = self.B(0,6,'g')


      # read perturbations
      W_AB_set = numpy.fromfile('temp_W_AB_set.dat').reshape(s,n,n)
      W_BA_set = numpy.fromfile('temp_W_BA_set.dat').reshape(s,n,n)
      w_AB_set = numpy.fromfile('temp_w_AB_set.dat').reshape(s,n,n)
      w_BA_set = numpy.fromfile('temp_w_BA_set.dat').reshape(s,n,n)
      F_A_set = numpy.fromfile('temp_F_A_set.dat').reshape(s,N,3)
      F_B_set = numpy.fromfile('temp_F_B_set.dat').reshape(s,N,3)

      if self._use_external_field_model:
         F_extField_set = numpy.fromfile('temp_F_extField_set.dat').reshape(s,N,3)           
         V_extField_set = numpy.fromfile('temp_V_extField_set.dat').reshape(s,n,n)
         E_nuc_field_set= numpy.fromfile('temp_E_nuc_field_set.dat')
         dD_extField_set_ref = numpy.fromfile('temp_dD_extField_ref_set.dat').reshape(s,n,n)
         dG_extField_set_ref = numpy.fromfile('temp_dG_extField_ref_set.dat').reshape(s,n,n)

      S_AB_set = numpy.fromfile('temp_S_AB_set.dat').reshape(s,n,n)

      H_set = numpy.fromfile('temp_H_set.dat').reshape(s,n*2,n*2)
      T_set = numpy.fromfile('temp_T_set.dat').reshape(s,n*2,n*2)
      V_set = numpy.fromfile('temp_V_set.dat').reshape(s,n*2,n*2)
      DIP_set = numpy.fromfile('temp_DIP_set.dat').reshape(s,3,n*2,n*2)
      e_set = numpy.fromfile('temp_e_set.dat')
      E_nuc_set = numpy.fromfile('temp_E_nuc_set.dat')
      DIP_nuc_set = numpy.fromfile('temp_DIP_nuc_set.dat').reshape(s,3)

      dD_AA_set_ref = numpy.fromfile('temp_dD_AA_ref_set.dat').reshape(s,n,n)
      dD_BB_set_ref = numpy.fromfile('temp_dD_BB_ref_set.dat').reshape(s,n,n)
      dD_AB_set_ref = numpy.fromfile('temp_dD_AB_ref_set.dat').reshape(s,n,n)

      dG_AA_set_ref = numpy.fromfile('temp_dG_AA_ref_set.dat').reshape(s,n,n)
      dG_BB_set_ref = numpy.fromfile('temp_dG_BB_ref_set.dat').reshape(s,n,n)
      dG_AB_set_ref = numpy.fromfile('temp_dG_AB_ref_set.dat').reshape(s,n,n)

      psi4.core.print_out(" ---> DMSFit: Check - training set <---\n\n")
      for s in range(self._nsamples):

          self._i = s+1

          S_AB = S_AB_set[s]; S_BA = S_AB.T.copy()
          W_AB = W_AB_set[s]
          W_BA = W_BA_set[s]
          w_AB = w_AB_set[s]
          w_BA = w_BA_set[s]
          #
          W_ABA= W_AB @ S_BA 
          W_BAB= W_BA @ S_AB
          w_ABA= w_AB @ S_BA
          w_BAB= w_BA @ S_AB
          #
          W_ABAB= W_ABA @ S_AB
          W_BABA= W_BAB @ S_BA
          w_ABAB= w_ABA @ S_AB
          w_BABA= w_BAB @ S_BA
          #
          W_ABABA= W_ABAB @ S_BA
          W_BABAB= W_BABA @ S_AB
          w_ABABA= w_ABAB @ S_BA
          w_BABAB= w_BABA @ S_AB
          #
          W_ABABAB= W_ABABA @ S_AB
          W_BABABA= W_BABAB @ S_BA
          w_ABABAB= w_ABABA @ S_AB
          w_BABABA= w_BABAB @ S_BA
          #
          W_ABABABA= W_ABABAB @ S_BA
          W_BABABAB= W_BABABA @ S_AB
          w_ABABABA= w_ABABAB @ S_BA
          w_BABABAB= w_BABABA @ S_AB

          F_A  = F_A_set[s]
          F_B  = F_B_set[s]


          # One-electron Fock matrix components
          H = H_set[s].copy()
          T = T_set[s].copy()
          V = V_set[s].copy()
          #
          H_AA = H_set[s][:n,:n]
          H_BB = H_set[s][n:,n:]
          H_AB = H_set[s][:n,n:]

          # Monomer Fock matrix and its components
          F_0   = self._wfn_0.Fa().to_array(dense=True)
          H_0   = self._wfn_0.H ().to_array(dense=True)
          T_0   = T[:n,:n]
          V_0   = H_0 - T_0
          G_0   = F_0 - H_0

          # Reference difference matrices
          dD_AA_ref = dD_AA_set_ref[s]
          dD_BB_ref = dD_BB_set_ref[s]
          dD_AB_ref = dD_AB_set_ref[s]

          dG_AA_ref = dG_AA_set_ref[s]
          dG_BB_ref = dG_BB_set_ref[s]
          dG_AB_ref = dG_AB_set_ref[s]
 
          #
          if self._use_external_field_model:
             F_extField  = F_extField_set[s]                    
             V_extField  = V_extField_set[s]
             E_nuc_field_mon = E_nuc_field_set[s]
             E_nuc_mon_0 = self._mol.nuclear_repulsion_energy()
             dD_extField_ref = dD_extField_set_ref[s]
             dG_extField_ref = dG_extField_set_ref[s]

             dD_extField_com = numpy.einsum("acix,ix->ac", B_10, F_extField)
             dD_extField_com+= numpy.einsum("acixy,ix,iy->ac", B_20, F_extField, F_extField)

             dG_extField_com = numpy.einsum("acix,ix->ac", b_10, F_extField)
             dG_extField_com+= numpy.einsum("acixy,ix,iy->ac", b_20, F_extField, F_extField)

             self._save_dD(dD_extField_com, prefix='d_extfield_com') 
             self._save_dD(dD_extField_ref, prefix='d_extfield_ref')
                                                                     
             self._save_dD(dG_extField_com, prefix='g_extfield_com')
             self._save_dD(dG_extField_ref, prefix='g_extfield_ref')

             D_extField_com = dD_extField_com + self._D0                       
             G_extField_com = dG_extField_com + self._dms_extField_g._M.copy()
             D_extField_ref = dD_extField_ref + self._D0
             G_extField_ref = dG_extField_ref + self._dms_extField_g._M.copy()
                                                                               
             mints = psi4.core.MintsHelper(self._bfs_0)
             V_mon = mints.ao_potential().to_array(dense=True)
             T_mon = mints.ao_kinetic().to_array(dense=True)
             U_mon = V_extField_set[s]
             H_mon = T_mon + V_mon + U_mon
                                                                               
             F_extField_com = G_extField_com + H_mon 
             F_extField_ref = G_extField_ref + H_mon 

             E_extField_com = (D_extField_com @ (H_mon + F_extField_com)).trace() + E_nuc_mon_0 + E_nuc_field_mon      
             E_extField_ref = (D_extField_ref @ (H_mon + F_extField_ref)).trace() + E_nuc_mon_0 + E_nuc_field_mon
                                                                                                                       
             dg_extField_com = F_extField_com - U_mon - F_0 # dg = F' - U - F_0
             dg_extField_ref = F_extField_ref - U_mon - F_0 # dg = F' - U - F_0
                                                                                                                       
             dE_extField_pol_com = (self._D0 @ dg_extField_com).trace() \
                                  +((D_extField_com - self._D0) @ (H_mon + F_extField_com)).trace()
             dE_extField_pol_ref = (self._D0 @ dg_extField_ref).trace() \
                                  +((D_extField_ref - self._D0) @ (H_mon + F_extField_ref)).trace()
                                                                                                                       
             # compute Fock matrices only from OPDM and ERIs (for test)
             I = numpy.identity(n)
             jk = psi4.core.JK.build(self._bfs_0, jk_type="direct")
             jk.set_memory(int(5e8))
             jk.initialize()
             jk.C_clear()                                           
             jk.C_left_add(psi4.core.Matrix.from_array(D_extField_com, ""))
             jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
             jk.compute()
             J = jk.J()[0].to_array(dense=True)
             K = jk.K()[0].to_array(dense=True)
             jk.finalize()
             psi4.core.clean()
             F_extField_test = H_mon + 2.0 * J - K
                                                                                                                       
             dg_extField_test = F_extField_test - U_mon - F_0 # dg = F' - U - F_0
             dE_extField_pol_test= (self._D0 @ dg_extField_test).trace() \
                                  +((D_extField_com - self._D0) @ (H_mon + F_extField_test)).trace()
                                                                                                                       
             psi4.core.print_out("ExtField: E_com    = %14.6f E_ref    = %14.6f\n" % (E_extField_com, E_extField_ref))
             psi4.core.print_out("ExtField: E_pol_com= %14.6E E_pol_tes= %14.6E E_pol_ref= %14.6E\n" % \
                       (dE_extField_pol_com, dE_extField_pol_test, dE_extField_pol_ref))
             #


          #
          dD_AA_com = numpy.einsum("acix,ix->ac", B_10, F_B)
          dD_AA_com+= numpy.einsum("acixy,ix,iy->ac", B_20, F_B, F_B)
          if (0,2) in self._dms_da.available_orders():
              dD_AA_com+= W_ABA @ B_02 + B_02.T @ W_ABA.T 
          if (0,4) in self._dms_da.available_orders():
              dD_AA_com+= W_ABABA @ B_04 + B_04.T @ W_ABABA.T 
          if (0,6) in self._dms_da.available_orders():
              dD_AA_com+= W_ABABABA @ B_06 + B_06.T @ W_ABABABA.T 

          dD_BB_com = numpy.einsum("acix,ix->ac", B_10, F_A)
          dD_BB_com+= numpy.einsum("acixy,ix,iy->ac", B_20, F_A, F_A)
          if (0,2) in self._dms_da.available_orders():
              dD_BB_com+= W_BAB @ B_02 + B_02.T @ W_BAB.T 
          if (0,4) in self._dms_da.available_orders():
              dD_BB_com+= W_BABAB @ B_04 + B_04.T @ W_BABAB.T 
          if (0,6) in self._dms_da.available_orders():
              dD_BB_com+= W_BABABAB @ B_06 + B_06.T @ W_BABABAB.T 

          #
          dG_AA_com = numpy.einsum("acix,ix->ac", b_10, F_B)
          dG_AA_com+= numpy.einsum("acixy,ix,iy->ac", b_20, F_B, F_B)
          if (0,2) in self._dms_da.available_orders():
              dG_AA_com+= w_ABA @ b_02 + b_02.T @ w_ABA.T 
          if (0,4) in self._dms_da.available_orders():
              dG_AA_com+= w_ABABA @ b_04 + b_04.T @ w_ABABA.T 
          if (0,6) in self._dms_da.available_orders():
              dG_AA_com+= w_ABABABA @ b_06 + b_06.T @ w_ABABABA.T 


          dG_BB_com = numpy.einsum("acix,ix->ac", b_10, F_A)
          dG_BB_com+= numpy.einsum("acixy,ix,iy->ac", b_20, F_A, F_A)
          if (0,2) in self._dms_da.available_orders():
              dG_BB_com+= w_BAB @ b_02 + b_02.T @ w_BAB.T 
          if (0,4) in self._dms_da.available_orders():
              dG_BB_com+= w_BABAB @ b_04 + b_04.T @ w_BABAB.T 
          if (0,6) in self._dms_da.available_orders():
              dG_BB_com+= w_BABABAB @ b_06 + b_06.T @ w_BABABAB.T 

          #
          #
          dD_AB_com = W_AB @ B_01 + B_01.T @ W_BA.T
          dG_AB_com = w_AB @ b_01 + b_01.T @ w_BA.T
          #
          #
          if (0,3) in self._dms_da.available_orders():
              dD_AB_com += W_ABAB @ B_03 + B_03.T @ W_BABA.T
              dG_AB_com += w_ABAB @ b_03 + b_03.T @ w_BABA.T
          if (0,5) in self._dms_da.available_orders():
              dD_AB_com += W_ABABAB @ B_05 + B_05.T @ W_BABABA.T
              dG_AB_com += w_ABABAB @ b_05 + b_05.T @ w_BABABA.T

          if 0:
          #dD_AA_com = dD_AA_ref.copy()
          #dD_BB_com = dD_BB_ref.copy()
           dG_AA_com = dG_AA_ref.copy()
           dG_BB_com = dG_BB_ref.copy()
          #dD_AB_com.fill(0.0)
          #dG_AB_com.fill(0.0)


          # block AA
          rms_da = self._rms(dD_AA_ref, dD_AA_com)
          rms_g  = self._rms(dG_AA_ref, dG_AA_com)

          t1_da = (dD_AA_ref @ H_AA).trace()
          t2_da = (dD_AA_com @ H_AA).trace()
          t1_g  = (dG_AA_ref @ H_AA).trace()
          t2_g  = (dG_AA_com @ H_AA).trace()

          psi4.core.print_out(" Sample=%03d AA Da RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_da, dD_AA_ref[0,0], dD_AA_com[0,0], t1_da, t2_da))
          psi4.core.print_out(" Sample=%03d AA G  RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_g , dG_AA_ref[0,0], dG_AA_com[0,0], t1_g , t2_g ))


          # block BB
          rms_da = self._rms(dD_BB_ref, dD_BB_com)
          rms_g  = self._rms(dG_BB_ref, dG_BB_com)

          t1_da = (dD_BB_ref @ H_BB).trace()
          t2_da = (dD_BB_com @ H_BB).trace()
          t1_g  = (dG_BB_ref @ H_BB).trace()
          t2_g  = (dG_BB_com @ H_BB).trace()

          psi4.core.print_out(" Sample=%03d BB Da RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_da, dD_BB_ref[0,0], dD_BB_com[0,0], t1_da, t2_da))
          psi4.core.print_out(" Sample=%03d BB G  RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_g , dG_BB_ref[0,0], dG_BB_com[0,0], t1_g , t2_g ))


          # block AB
          rms_da = self._rms(dD_AB_ref, dD_AB_com)
          rms_g  = self._rms(dG_AB_ref, dG_AB_com)

          t1_da = (dD_AB_ref @ H_AB).trace()
          t2_da = (dD_AB_com @ H_AB).trace()
          t1_g  = (dG_AB_ref @ H_AB).trace()
          t2_g  = (dG_AB_com @ H_AB).trace()

          psi4.core.print_out(" Sample=%03d AB Da RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_da, dD_AB_ref[0,0], dD_AB_com[0,0], t1_da, t2_da))
          psi4.core.print_out(" Sample=%03d AB G  RMS=%13.5E  %14.5f %14.5f  Tr[dD,H]=%14.5f  %14.5f\n" \
                   % (s+1, rms_g , dG_AB_ref[0,0], dG_AB_com[0,0], t1_g , t2_g ))


          self._save_dD(dD_AA_ref, prefix='d_aa_ref')
          self._save_dD(dD_AA_com, prefix='d_aa_com')
          self._save_dD(dG_AA_ref, prefix='g_aa_ref')
          self._save_dD(dG_AA_com, prefix='g_aa_com')

          self._save_dD(dD_BB_ref, prefix='d_bb_ref')
          self._save_dD(dD_BB_com, prefix='d_bb_com')
          self._save_dD(dG_BB_ref, prefix='g_bb_ref')
          self._save_dD(dG_BB_com, prefix='g_bb_com')

          self._save_dD(dD_AB_ref, prefix='d_ab_ref')
          self._save_dD(dD_AB_com, prefix='d_ab_com')
          self._save_dD(dG_AB_ref, prefix='g_ab_ref')
          self._save_dD(dG_AB_com, prefix='g_ab_com')

          # reconstruct whole matrices
          dD_com = numpy.zeros((n+n,n+n))
          dG_com = numpy.zeros((n+n,n+n))
          dD_ref = numpy.zeros((n+n,n+n))
          dG_ref = numpy.zeros((n+n,n+n))

          dD_com[:n,:n] = dD_AA_com.copy() 
          dD_com[n:,n:] = dD_BB_com.copy()
          dD_com[:n,n:] = dD_AB_com.copy()
          dD_com[n:,:n] = dD_AB_com.copy().T

          dD_ref[:n,:n] = dD_AA_ref.copy()
          dD_ref[n:,n:] = dD_BB_ref.copy()
          dD_ref[:n,n:] = dD_AB_ref.copy()
          dD_ref[n:,:n] = dD_AB_ref.copy().T

          dG_com[:n,:n] = dG_AA_com.copy()
          dG_com[n:,n:] = dG_BB_com.copy()
          dG_com[:n,n:] = dG_AB_com.copy()
          dG_com[n:,:n] = dG_AB_com.copy().T

          dG_ref[:n,:n] = dG_AA_ref.copy()
          dG_ref[n:,n:] = dG_BB_ref.copy()
          dG_ref[:n,n:] = dG_AB_ref.copy()
          dG_ref[n:,:n] = dG_AB_ref.copy().T

          #
          D_com = dD_com.copy()
          D_com[:n,:n]+= self._D0
          D_com[n:,n:]+= self._D0


          G_com = dG_com.copy()
          G_com[:n,:n]+= self._G0
          G_com[n:,n:]+= self._G0

          D_ref = dD_ref.copy()
          D_ref[:n,:n]+= self._D0
          D_ref[n:,n:]+= self._D0

          G_ref = dG_ref.copy()
          G_ref[:n,:n]+= self._G0
          G_ref[n:,n:]+= self._G0

          # test only off-diagonals
          if 0:
             D_com[:n,:n].fill(0.0) 
             D_com[n:,n:].fill(0.0)
             D_ref[:n,:n].fill(0.0)
             D_ref[n:,n:].fill(0.0)
                                    
             G_com[:n,:n].fill(0.0)
             G_com[n:,n:].fill(0.0)
             G_ref[:n,:n].fill(0.0)
             G_ref[n:,n:].fill(0.0)

             if 0:
                D_com[:n,:n] = self._D0 
                D_com[n:,n:] = self._D0
                D_ref[:n,:n] = self._D0
                D_ref[n:,n:] = self._D0
                                       
                G_com[:n,:n] = self._G0
                G_com[n:,n:] = self._G0
                G_ref[:n,:n] = self._G0
                G_ref[n:,n:] = self._G0


          # compute U_diff
          U_AA = V[:n,:n] - V_0
          U_BB = V[n:,n:] - V_0
          U_diff = numpy.zeros((2*n,2*n))
          U_diff[:n,:n] = U_AA
          U_diff[n:,n:] = U_BB
         #U_diff[:n,n:] = V[:n,n:].copy()
         #U_diff[n:,:n] = V[n:,:n].copy()
          H_0_diag = numpy.zeros((2*n,2*n))
          H_0_diag[:n,:n] = H_0.copy()
          H_0_diag[n:,n:] = H_0.copy()

          if 'e' in self._dms_types: # G = 2J - K
             #F_com = G_com + H
             #F_ref = G_ref + H
              F_com = G_com + H #U_diff + H_0_diag
              F_ref = G_ref + H #U_diff + H_0_diag
          elif 'g' in self._dms_types: # G = 2J - K + U
              F_ref = G_ref + H - U_diff
              F_com = G_com + H - U_diff
          #    F_com = G_com + T
          #    F_ref = G_ref + T 
          #    # G = [V + U] + 2J - K
          #   #F_extField_com = G_extField_com + T_mon
          #   #F_extField_ref = G_extField_ref + T_mon
          #   #F_extField_com = F_extField_com - U_mon
          #   #F_extField_ref = F_extField_ref - U_mon
          #elif 'f' in self._dms_types: # G = T + V + 2J - K = F
          #    F_com = G_com  
          #    F_ref = G_ref  
          #    # G = T + [V + U] + 2J - K = F'
          #   #F_extField_com = G_extField_com 
          #   #F_extField_ref = G_extField_ref
          #   #F_extField_com = F_extField_com - U_mon
          #   #F_extField_ref = F_extField_ref - U_mon
          #elif 'g1' in self._dms_types: # G = 2V + 2J - K
          #    F_com = G_com - V + T
          #    F_ref = G_ref - V + T
          #    # G = #TODO
          #   #F_extField_com = G_extField_com - V_mon + T_mon
          #   #F_extField_ref = G_extField_ref - V_mon + T_mon
          #   #F_extField_com = F_extField_com - U_mon
          #   #F_extField_ref = F_extField_ref - U_mon

          if 1:
             # test total Fock matrix from dD_com and ERI's
             I = numpy.identity(2*n)
             #bfs_dimer = self._bfs_dimer_set[s]
             dimer = psi_molecule_from_file("geom_%03d.xyz" % self._i)
             bfs_dimer = psi4.core.BasisSet.build(dimer, "BASIS",psi4.core.get_global_option("BASIS"), puream=False)
             jk = psi4.core.JK.build(bfs_dimer, jk_type="direct")
             jk.set_memory(int(5e8))
             jk.initialize()
             jk.C_clear()                                           
             jk.C_left_add(psi4.core.Matrix.from_array(D_com, ""))
             jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
             jk.compute()
             J = jk.J()[0].to_array(dense=True)
             K = jk.K()[0].to_array(dense=True)
             jk.finalize()
             psi4.core.clean()
             F_com_test = H + 2.0 * J - K


          # compute total energy (MCBS)
          dE_com = (D_com @ (H + F_com)).trace() + E_nuc_set[s] - 2.0 * self._e0
          dE_ref = (D_ref @ (H + F_ref)).trace() + E_nuc_set[s] - 2.0 * self._e0
          dE_ref_test = e_set[s] - 2.0 * self._e0
          if 1: 
             dE_com_test = (D_com @ (H + F_com_test)).trace() + E_nuc_set[s] - 2.0 * self._e0


          psi4.core.print_out(" Sample=%03d Dimer Energy %14.8f  %14.8f  %14.8f  %14.8f\n" \
                   % (s+1, dE_com, dE_com_test, dE_ref, dE_ref_test))

          self._save_dD(dD_com, prefix='dd_com')
          self._save_dD(dD_ref, prefix='dd_ref')
          self._save_dD(dG_com, prefix='dg_com')
          self._save_dD(dG_ref, prefix='dg_ref')

          # compute dipole moment
          DIP = DIP_set[s]
          DIP_nuc = DIP_nuc_set[s]
          mu_x_ref = 2.0 * (DIP[0] @ D_ref).trace() + DIP_nuc[0]
          mu_y_ref = 2.0 * (DIP[1] @ D_ref).trace() + DIP_nuc[1]
          mu_z_ref = 2.0 * (DIP[2] @ D_ref).trace() + DIP_nuc[2]
          mu_x_com = 2.0 * (DIP[0] @ D_com).trace() + DIP_nuc[0]
          mu_y_com = 2.0 * (DIP[1] @ D_com).trace() + DIP_nuc[1]
          mu_z_com = 2.0 * (DIP[2] @ D_com).trace() + DIP_nuc[2]

          psi4.core.print_out(" Sample=%03d Mu X %14.5f  %14.5f\n" \
                   % (s+1, mu_x_com, mu_x_ref))
          psi4.core.print_out(" Sample=%03d Mu Y %14.5f  %14.5f\n" \
                   % (s+1, mu_y_com, mu_y_ref))
          psi4.core.print_out(" Sample=%03d Mu Z %14.5f  %14.5f\n" \
                   % (s+1, mu_z_com, mu_z_ref))




  # -----> Private Utility methods <----- #

  def __draw_translation(self):
      theta = numpy.arccos(self._random_double())
      phi   = 2.0 * numpy.pi * self._random_double()
      r     = self._start + self._range * numpy.cbrt(self._random_double())
      x     = r * numpy.sin(theta) * numpy.cos(phi)
      y     = r * numpy.sin(theta) * numpy.sin(phi)
      z     = r * numpy.cos(theta)                 

      t     = numpy.array([x,y,z])
      return t

  def __new_translation(self):
      "Select new random translation"
      done = False
      geom = self._mol.geometry().to_array(dense=True) 
      if DMSFit.generate_random_samples:
         while not done:                      
            t = self.__draw_translation()
            if not self._clash(geom, geom+t):
               done = True 
      else:
          r = DMSFit.minimum_atom_atom_distance + self._i * self._range/self._nsamples
          t = r * numpy.array([1.0,0.0,0.0])
      return t



class Rotation_DMSFit(ExternalField_EFP_DMSFit):
  """
 Rotation method to fit DMS tensors.
"""
  def __init__(self, mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model,
                     start=1.0, srange=2.0):
      super().__init__(mol, method, nsamples, dms_types, order_type, use_iterative_model, use_external_field_model)

      # translation parameters
      self._start = start
      self._range = srange

      raise NotImplementedError("Rotation DMSFit method is not developed yet")

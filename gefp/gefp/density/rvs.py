#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Reduced Variational Scheme of Stevens and Fink module.
 Bartosz Błasiak, Gundelfingen, 14 Feb 2020
"""
import psi4
import numpy
import gefp
from ..math.orthonorm import GrammSchmidt

__all__ = ["RVS"]

class RVS:
  """
 -----------------------------------------------------------------------------
 Reduced Variational Scheme (RVS)
 Ref.: Stevens & Fink, Chem. Phys. Lett., Vol. 139, pp. 15-22 (1987)

 Implementation for closed-shell systems.
 -----------------------------------------------------------------------------

 Notes:

 o This code contains general implementation of the RVS-SCF method.
   Currently, only dimers are automatically analyzed. To use for
   multimers, define a subclass and redefine the `run` method
   to include n-mers for n>2.

 o Projection of frozen orbitals is achieved by transforming
   Fock matrix to orthogonal MO basis of entire n-fragment 
   aggregate, and setting all the off-diagonal matrix elements
   that are associated to the frozen orbitals to zero 
   [Kairys & Jensen, J. Phys. Chem. A, Vol. 104, No. 28, 2000].

 o Orthogonal MO basis is constructed from the original
   mutually non-orthogonal MOs of isolated fragments by 
   GrammSchmidt orthogonalization with respect to the frozen
   MOs. If no frozen MOs are requested, symmetric Lowdin
   orthogonalization is performed instead.

 -----------------------------------------------------------------------------
 Usage for dimers:

   e, wfn = psi4.energy('scf', return_wfn=True)
   rvs = RVS(wfn)
   rvs.run_dimer(conver=1.0e-7, maxiter=100, ndamp=10)
   print(rvs)
   # access variables:
   print(rvs.vars.keys())

 -----------------------------------------------------------------------------
 Example for redefining for multimers:

 class Trimer_RVS(RVS):
    def __init__(self, wfn):
        super(RVS, self).__init__(wfn)

    def run_trimer(self, conver, maxiter, ndamp):
        "Here there is your implementation for trimer" 

        # example for A-B-C (0-1-2) trimer:

        # energy of Hartree product
        #                        frozen_occ active_vir exclude_occ
        E_Ao_Bo_Co        = self._scf([ ], [   ], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label=' Ao  Bo Co')

        # energy of Bocc frozen, virtual space composed of Cvir and Bvir, and A molecule excluded
        E_fBo_Bv_Cv       = self._scf([1], [1,2], [0], conver, maxiter, damp=0.0, ndamp=ndamp, label='[Bo] Bv Cv')

        # energy of Bocc and Cocc frozen and all virtual space included, Aocc active
        E_fBo_fCo_Av_Bv_Cv= self._scf([1,2], [0,1,2], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label='Ao [Bo] [Co] Av Bv Cv')

    def __repr__(self):
        "Print the contents"
        log = ''
        # ...
        return str(log)

 -----------------------------------------------------------------------------
 B. Błasiak                                       Gundelfingen, 14.02.2020
"""
  def __init__(self, wfn):
      "Initialize RVS solver"

      # sanity check
      if wfn.molecule().nfragments() == 1:
         raise ValueError(" RVS Error: There must be at least 2 fragments in the system.")

      # basic data
      self._wfn = wfn
      self._bfs = wfn.basisset()
      self._mol = wfn.molecule()
      self._nfrag= self._mol.nfragments()

      # RVS components are stored here
      self.vars = {}

      # Orbitals in MCBS per each fragment (constant)
      self._c_occ_0 = []
      self._c_vir_0 = []

      # Isolated fragments (constant)
      self._e_0 = []
      self._frags = []

      # initialize JK object
      self._global_jk = psi4.core.JK.build(self._bfs, jk_type='direct')
      self._global_jk.set_memory(int(5e8))
      self._global_jk.initialize()

      # constant matrices in DCBS
      self._H_core_global = wfn.H().to_array(dense=True)
      self._S_global = self._wfn.S().to_array(dense=True)
      self._X = gefp.math.matrix.matrix_power(self._S_global,-0.5)
      self._Y = gefp.math.matrix.matrix_power(self._S_global, 0.5)

      mints = psi4.core.MintsHelper(self._bfs)
      self._T_kinetic_global = mints.ao_kinetic().to_array(dense=True)


  # ----> Public Interface <----

  def run_dimer(self, conver=1.0e-7, maxiter=100, ndamp=10): 
      "Run the RVS analysis for dimer"
      conv = psi4.constants.hartree2kcalmol

      # compute monomers
      print(" @RVS-SCF. Start.")
      self._isolated_monomers()

      e = sum(self._e_0)
      e1, e2 = self._e_0
      print(" E(1) = %16.6f [AU]" % (e1))
      print(" E(2) = %16.6f [AU]" % (e2))

      # RVS energy components    -Focc-  -Avir-  -Xocc-

     #E_A               = self._scf([ ], [0  ], [1], conver, maxiter, damp=0.0, ndamp=ndamp, label=' A')   
     #E_B               = self._scf([ ], [1  ], [0], conver, maxiter, damp=0.0, ndamp=ndamp, label=' B')   

      E_Ao_Bo           = self._scf([ ], [   ], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label=' Ao  Bo') # HP energy
      # 
      E_fAo_Bo_Bv       = self._scf([0], [1  ], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label='[Ao] Bo Bv')
      E_fAo_Bo_Bv_Av    = self._scf([0], [0,1], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label='[Ao] Bo Bv Av')
      #
      E_fBo_Ao_Av       = self._scf([1], [0  ], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label='[Bo] Ao Av')
      E_fBo_Ao_Av_Bv    = self._scf([1], [0,1], [ ], conver, maxiter, damp=0.0, ndamp=ndamp, label='[Bo] Ao Av Bv')
      #
      E_Ao_Av_Bv        = self._scf([ ], [0,1], [1], conver, maxiter, damp=0.0, ndamp=ndamp, label=' Ao  Av Bv')
      E_Bo_Av_Bv        = self._scf([ ], [0,1], [0], conver, maxiter, damp=0.0, ndamp=ndamp, label=' Bo  Av Bv')

      print(" DE( Ao  Bo         ) = %16.6f [kcal/mol]" % ((E_Ao_Bo       -e) * conv))
      print(" DE([Ao] Bo Bv      ) = %16.6f [kcal/mol]" % ((E_fAo_Bo_Bv   -e) * conv))
      print(" DE([Ao] Bo Bv Av   ) = %16.6f [kcal/mol]" % ((E_fAo_Bo_Bv_Av-e) * conv))
      print(" DE([Bo] Ao Av      ) = %16.6f [kcal/mol]" % ((E_fBo_Ao_Av   -e) * conv))
      print(" DE([Bo] Ao Av Bv   ) = %16.6f [kcal/mol]" % ((E_fBo_Ao_Av_Bv-e) * conv))
      print(" DE( Ao  Av Bv      ) = %16.6f [kcal/mol]" % ((E_Ao_Av_Bv    -e1) * conv))
      print(" DE( Bo  Av Bv      ) = %16.6f [kcal/mol]" % ((E_Bo_Av_Bv    -e2) * conv))

      # BSSE-uncorrected interaction energy
      INT = self._wfn.energy() - e

      # RVS components
      CEX = E_Ao_Bo - e
      POL_B = E_fAo_Bo_Bv         - e - CEX
      SUP_B = E_Bo_Av_Bv          - e2
      CT_BA = E_fAo_Bo_Bv_Av      - e - CEX - POL_B - SUP_B
      POL_A = E_fBo_Ao_Av         - e - CEX
      SUP_A = E_Ao_Av_Bv          - e1
      CT_AB = E_fBo_Ao_Av_Bv      - e - CEX - POL_A - SUP_A

      POL= POL_A + POL_B
      SUP= SUP_A + SUP_B
      CT = CT_AB + CT_BA

      TOT = CEX + POL + SUP + CT

      # save
      self.vars["e1"] = e1
      self.vars["e2"] = e2
      self.vars["etot"] = e
      self.vars["cex"] = CEX
      self.vars["pol_a"] = POL_A
      self.vars["pol_b"] = POL_B
      self.vars["ct_ab"] = CT_AB
      self.vars["ct_ba"] = CT_BA
      self.vars["sup_a"] = SUP_A
      self.vars["sup_b"] = SUP_B
      self.vars["ct"   ] = CT
      self.vars["pol"  ] = POL
      self.vars["sup"  ] = SUP
      self.vars["tot"  ] = TOT
      self.vars["int"  ] = INT

      print(self)

  def __repr__(self):
      "Print the contents. Now programmed only for dimer"
      conv = psi4.constants.hartree2kcalmol

      if self._nfrag == 2:

         CEX   = self.vars["cex"]
         POL_A = self.vars["pol_a"]
         POL_B = self.vars["pol_b"]
         CT_AB = self.vars["ct_ab"]
         CT_BA = self.vars["ct_ba"]
         SUP_A = self.vars["sup_a"]
         SUP_B = self.vars["sup_b"]
         POL   = self.vars["pol"  ]
         CT    = self.vars["ct"   ]
         SUP   = self.vars["sup"  ]
         TOT   = self.vars["tot"  ]
         INT   = self.vars["int"  ]
                                                                        
         log =  "\n"
         log+=  " =========================================\n"
         log+=  "               RVS ANALYSIS               \n"
         log+=  " -----------------------------------------\n"
         log+=  "          [a.u.]         [kcal/mol]       \n"
         log+=  " -----------------------------------------\n"
         log+=  " CEX   = %16.6f %16.6f\n" % (CEX  , CEX * conv)
         log+=  " -----------------------------------------\n"
         log+=  " POL_A = %16.6f %16.6f\n" % (POL_A, POL_A * conv)
         log+=  " POL_B = %16.6f %16.6f\n" % (POL_B, POL_B * conv)
         log+=  " CT_AB = %16.6f %16.6f\n" % (CT_AB, CT_AB * conv)
         log+=  " CT_BA = %16.6f %16.6f\n" % (CT_BA, CT_BA * conv)
         log+=  " SUP_A = %16.6f %16.6f\n" % (SUP_A, SUP_A * conv)
         log+=  " SUP_B = %16.6f %16.6f\n" % (SUP_B, SUP_B * conv)
         log+=  " - - - - - - - - - - - - - - - - - - - - -\n"
         log+=  " POL   = %16.6f %16.6f\n" % (POL  , POL   * conv)
         log+=  " CT    = %16.6f %16.6f\n" % (CT   , CT    * conv)
         log+=  " SUP   = %16.6f %16.6f\n" % (SUP  , SUP   * conv)
         log+=  " -----------------------------------------\n"         
         log+=  " TOT   = %16.6f %16.6f\n" % (TOT  , TOT   * conv)
         log+=  " INT   = %16.6f %16.6f\n" % (INT  , INT   * conv)
         log+=  " =========================================\n\n"

      else: log = " RVS: Display not implemented for %d-mer\n" % self._nfrag

      return str(log)


  # ----> protected interface <----

  def _scf(self, frozen_occ, active_vir, exclude_occ, 
                 conver, maxiter, damp=0.0, ndamp=10, label=None): 
      "Solve RVS SCF problem for a given set of frozen occupied orbitals and active virtual orbitals"

      active_occ = [x for x in range(self._nfrag) if x not in frozen_occ+exclude_occ]
      all_occ    = [x for x in range(self._nfrag) if x not in            exclude_occ]

      # [1] Create AO-MO transformation matrices
      C_occ_active = numpy.hstack([self._c_occ_0[i].copy() for i in active_occ ])
      C_occ_all    = numpy.hstack([self._c_occ_0[i].copy() for i in     all_occ])

      C_act        = C_occ_active.copy() # nMO basis (AO x nMO)
      if active_vir:
         C_vir     = numpy.hstack([self._c_vir_0[i].copy() for i in active_vir ])
         C_act     = numpy.hstack([C_occ_active, C_vir]).copy()
         del C_vir

      naocc_active = C_occ_active.shape[1]
      naocc_frozen = 0

      C_occ_frozen = None
      if frozen_occ:
         C_occ_frozen = numpy.hstack([self._c_occ_0[i].copy() for i in frozen_occ ])
         naocc_frozen = C_occ_frozen.shape[1]
         C_ = self._orthogonalize(C_act, C_occ_frozen, method='gs')
      else:
         C_ = self._orthogonalize(C_act, C_occ_frozen, method='lowdin')

      naocc = naocc_active+naocc_frozen
      del C_act, C_occ_active, C_occ_frozen
        
      # [2] Iteration cycles
      if label is not None:
         print(" @RVS-SCF: Starting iterations of %s" % label)
      E_old = 1.0e+10
      H = self._H_core_ao(exclude_occ)
      F = self._Fock_ao(H, C_occ_all)                   # Fock matrix (AO x AO)         
      D = C_occ_all @ C_occ_all.T                       # OPDM (AO x AO)
      E = self._energy(H, F, D, exclude_occ)
      E_new = E
      F_new = F
      print(" @RVS-SCF Iter %3i.   E= %16.8E" % (0, E_new))

      Iter = 1
      while (abs(E_old-E_new) > conver):

          E_old = E_new
          F_old = F_new
                                                                                            
          f_ = C_.T @ F_old @ C_                            # Fock in oMO basis (oMO x oMO)
          self._project_out(f_, naocc_frozen)
          k_, a_ = numpy.linalg.eigh(f_)                    # Orbitals (oMO x O-MO) 
     #    a = x @ a_                                        # Orbitals (nMO x O-MO)
          A = C_@ a_                                        # Orbitals (AO x O-MO)
          C_occ_all = A[:,:naocc]                           # Orbitals (AO x O-MO_occ)
          D_new = C_occ_all @ C_occ_all.T                   # OPDM (AO x AO)
          F_new = self._Fock_ao(H, C_occ_all)
          if Iter <= ndamp:
             F_new = damp * F_old + (1.0 - damp) * F_new
          E_new = self._energy(H, F_new, D_new, exclude_occ)

          print(" @RVS-SCF Iter %3i.   E= %16.8E  Delta= %16.8E" % (Iter, E_new, E_new - E_old))
          Iter += 1

          if (Iter >= maxiter):
              raise ValueError(" @RVS-SCF Error: Maximum iterations reached!\n")
              break

      if label is not None:
         print(" @RVS-SCF: Done with iterations of %s" % label)

      return E_new


  def _isolated_monomers(self):
      "Solve SCF for all isolated monomers in MCBS"
      nbf_all = self._bfs.nbf()
      nocc_all= 0
      nvir_all= 0

      off = 0
      for i in range(self._nfrag):
          print(" Calculation of isolated monomer %i" % (i+1))

          mol = self._mol.extract_subsets(i+1)
          e, w = psi4.energy('scf', molecule=mol, return_wfn=True)

          c_occ_i= w.Ca_subset("AO","OCC").to_array(dense=True)
          c_vir_i= w.Ca_subset("AO","VIR").to_array(dense=True)

          nbf_i= c_occ_i.shape[0]
          nocc = c_occ_i.shape[1]
          nvir = c_vir_i.shape[1]
          c_occ = numpy.zeros((nbf_all, nocc))
          c_vir = numpy.zeros((nbf_all, nvir))

          c_occ[off:off+nbf_i,:] = c_occ_i.copy()
          c_vir[off:off+nbf_i,:] = c_vir_i.copy()

          self._c_occ_0.append(c_occ)
          self._c_vir_0.append(c_vir)
          self._e_0.append(e)
          self._frags.append(mol)
 
          off += nbf_i
          nocc_all += nocc
          nvir_all += nvir

      self._nbf_all = nbf_all
      self._nocc_all = nocc_all

  def _project_out(self, f, n):
      "Project out first n orbitals from Fock matrix in oMO basis"
      fr= f[:n,:n]
      fr= numpy.diag(numpy.diag(fr))
      f[:n,:n] = fr.copy()
      f[:n,n:].fill(0.0)
      f[n:,:n].fill(0.0)

  def _Fock_ao(self, H_core, C_occ):
      "Compute Fock matrix in AO representation from AO-MO occupied orbital matrix"
      self._global_jk.C_clear()                                           
      self._global_jk.C_left_add(psi4.core.Matrix.from_array(C_occ, ""))
      self._global_jk.compute()
      J = self._global_jk.J()[0].to_array(dense=True)
      K = self._global_jk.K()[0].to_array(dense=True)

      F = H_core + 2.0 * J - K # Fock matrix (AO x AO)
      return F

  def _H_core_ao(self, exclude_occ):
      "Compute H-core Hamiltonian - implemented for multifragment case"
      include_mol = [x for x in range(self._nfrag) if x not in exclude_occ]
      if len(include_mol) == self._nfrag:
         H = self._H_core_global.copy()
      else:
         H = self._T_kinetic_global.copy()

         ep = psi4.core.ExternalPotential()
         BohrToAngstrom = 0.5291772086
         for i in include_mol:
             for a in range(self._frags[i].natom()):
                 q = self._frags[i].Z(a)
                 x = self._frags[i].x(a) * BohrToAngstrom
                 y = self._frags[i].y(a) * BohrToAngstrom
                 z = self._frags[i].z(a) * BohrToAngstrom
                 ep.addCharge(q, x, y, z)
         H += ep.computePotentialMatrix(self._bfs).to_array(dense=True)
      return H

  def _orthogonalize(self, C, C_frozen, method):
      "Orthogonalize sets of orbitals"
      if method == 'lowdin':
         if C_frozen is not None: 
            C = numpy.hstack([C_frozen, C])
         s = C.T @ self._S_global @ C                      # Overlap matrix (nMO x nMO)
         x = gefp.math.matrix.matrix_power(s, -0.5)        # Orthogonalizer (nMO x oMO) 
         C_ = C @ x                                        # LCAO-oMO (AO x oMO)
      elif method == 'gs':
         c_frozen = self._Y @ C_frozen
         c        = self._Y @ C
         schm = GrammSchmidt(c_frozen)
         for i in range(C.shape[1]):
             ci = c[:,i]
             vi = schm.orthogonalize_vector(ci, normalize=True)
             schm.append(vi)
         c_ = schm.V
         C_ = self._X @ c_
      else: raise NotImplementedError("Method %s is not available. Use only gs or lowdin." % method)
      return C_

  def _energy(self, H_core, F, D, exclude_occ, G_core=None):
      "Compute total energy from Fock and ODPM matrices in AO basis"
      if G_core is not None: H_core += G_core
      E = ((H_core + F) @ D).trace()  
      E+= self._energy_nuclear(exclude_occ)
      return E

  def _energy_nuclear(self, exclude_occ):
      "Nuclear repulsion energy"
      include_mol = [x+1 for x in range(self._nfrag) if x not in exclude_occ]
      if len(include_mol) == self._nfrag:
         e = self._mol.nuclear_repulsion_energy()
      else:
         mol = self._mol.extract_subsets(include_mol)
         e = mol.nuclear_repulsion_energy()
      return e

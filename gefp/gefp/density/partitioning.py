#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Partitioning module.
 Bartosz Błasiak, Gundelfingen, Dec 2018
"""

import sys
import math
import numpy
import numpy.linalg
import psi4
import oepdev
from .scf import SCF
from .opdm import Density

__all__ = ["DensityDecomposition"]


class DensityDecomposition(Density):
    """
 -------------------------------------------------------------------------------------------------------------

 Density-Based Decomposition Scheme of Mandado and Hermida-Ramon
 with partitioning of polarization deformation density into induction,
 dispersion and charge-transfer contributions.

                           --> DDS <--
               --> Density Decomposition Scheme <--

 References: 
  * Mandado and Hermida-Ramon, J. Chem. Theory Comput. 2011, 7, 633-641. (JCTC 2011)
  * Błasiak, J. Chem. Phys. 2018 149 (16), 164115. (JCP 2018)

 -------------------------------------------------------------------------------------------------------------

 Constructor arguments and options:

  o aggregate - Psi4 molecular aggregate with at least two fragments
  o method    - QM method (hf, mp2, cc2, ccsd)
  o acbs      - use aggregate-centred basis set for calculations of wavefunctions. 
                Otherwise use monomer-centred basis sets and composite Hadamard 
                addition of AO spaces. ACBS=False result in no correction for BSSE.
  o jk_type   - type of Psi4 JK object.
  o no_cutoff - cutoff for natural occupancies threshold. All natural orbitals
                with occupancies less or equal to the threshold will be neglected.
  o xc_scale  - scaling parameter for exchange-correlation density
  o l_dds     - compute also linear DDS total interaction energy
  o kwargs    - additional Psi4-relevant options.

 -------------------------------------------------------------------------------------------------------------

 Usage example:

   solver = DensityDecomposition(aggr, method='hf', acbs=True, jk_type='direct', no_cutoff=0.000, xc_scale=1.0, l_dds=True, **kwargs) 
   solver.compute(polar_approx=False)
                                                                                                            
   dD_pauli = solver.deformation_density('pau')
   dD_pol   = solver.deformation_density('pol')
   dD       = solver.deformation_density('fqm')

   # dictionaries:
   # 1. accessing variables
   solver.vars 
   # 2. accessing aggregate data
   solver.matrix
   # 3. accessing unperturbed fragment data (expert)
   solver.data

                                                                                                            
   print(solver)
   solver.print_out()

 -------------------------------------------------------------------------------------------------------------
                                                                      Last Revision: Gundelfingen, 11 Jan 2019
"""
    def __init__(self, aggregate, method='hf', acbs=True, jk_type='direct', no_cutoff=0.000, xc_scale=1.0, l_dds=True, 
                       cc_relax=True, taylor=False, erase_dpol_offdiag=False, verbose=False, xc_scale_test=False, **kwargs):
        "Initialize all attributes"
        # molecular aggregate
        self.aggregate   = aggregate    
        # level of theory
        self.method      = method       
        # unperturbed fragment data
        self.data        = {"wfn": [],  # unperturbed wavefunction
                            "ene": [],  # unperturbed total energy
                            "odm": [],  # unperturbed one-particle density matrix
                            "non": [],  # unperturbed natural occupations
                            "noc": [],  # unperturbed natural orbital wavefunction coefficients
                            "nmo": [],  # number of molecular orbitals
                            "nbf": [],  # number of basis functions 
                            "ofm": [],  # MO offset
                            "ofb": [],  # AO offset
                            } 
        # aggregate data: matrices in aggregate AO/MO basis
        self.matrix = {"d"   : None,    # unperturbed 1-electron density
                       "c"   : None,    # unperturbed LCAO-NO coefficients
                       "n"   : None,    # unperturbed NO occupation numbers
                       "doo" : None,    # orthogonalized unpolarized 1-electron density
                       "coo" : None,    # orthogonalized unpolarized LCAO-NO coefficients
                       "noo" : None,    # orthogonalized unpolarized NO occupation numbers
                       "dqm" : None,    # orthogonalized polarized 1-electron density (full QM 1-electron density)
                       "cqm" : None,    # orthogonalized polarized LCAO-NO coefficients (full QM)
                       "nqm" : None,    # orthogonalized polarized NO occupation numbers (full QM)
                       "dpp" : None,    # orthogonalized polarized 1-electron density (approximated)
                       "cpp" : None,    # orthogonalized polarized LCAO-NO coefficients (approximated)
                       "npp" : None,    # orthogonalized polarized NO occupation numbers (approximated)
                       "sqm" : None,    # overlap matrix in AO basis
                       "chf" : None,    # unperturbed LCAO-MO coefficients: HF level
                       }
        # variables
        self.vars        = {"e_cou_1"     : None,   # coulombic energy, 1el part
                            "e_cou_2"     : None,   # coulombic energy, 2el part
                            "e_cou_t"     : None,   # coulombic energy
                            "e_rep_1"     : None,   # repulsion energy, 1el part
                            "e_rep_2"     : None,   # repulsion energy, 2el part
                            "e_rep_t"     : None,   # repulsion energy, (1+2)el part
                            "e_exc_t"     : None,   # exchange energy
                            "e_exr_t"     : None,   # exchange-repulsion energy
                            "e_pol_t"     : None,   # polarization energy (exact)
                            "e_pol_a"     : None,   # polarization energy (approximated)
                            "e_pol_1"     : None,   # polarization energy, 1el part
                            "e_pol_2"     : None,   # polarization energy, 2el part
                            "e_pol_e"     : None,   # polarization energy, (1+2)el part
                            "e_exp_t"     : None,   # polarization energy, exchange-correlation part (exact)
                            "e_exp_a"     : None,   # polarization energy, exchange-correlation part (approximated)
                            "e_pol_ind_t" : None,   # induction part of polarization energy, (1+2)el part
                            "e_pol_disp_t": None,   # dispersion part of polarization energy, (1+2)el part
                            "e_pol_ct_t"  : None,   # charge-transfer part of polarization energy, (1+2)el part
                            "e_dds_t"     : 0.0 ,   # total interaction energy (DDS)
                            "e_dds_l_t"   : 0.0 ,   # total interaction energy (L-DDS)
                            "e_fqm_t"     : None,   # total interaction energy (full QM)
                            }
        # recommended XC density scaling factors
        self.xc_recommended = {"hf"          : 1.00000,
                               "mp2/6-311G**": 0.89   ,
                               }
        # additional options
        self.kwargs      = kwargs
        # basis-set mode (ACBS: aggregate-centred basis set)
        self.acbs        = acbs
        # cutoff for NO occupancies
        self.no_cutoff   = no_cutoff
        # relax density matrix for coupled cluster
        self.cc_relax = cc_relax
        # linear DDS
        self.l_dds = l_dds
        # verbose mode
        self.verbose = verbose
        # scaling factor of exchange-correlation occupation weight
        self.xc_scale = xc_scale
        self.xc_scale_test= xc_scale_test

        self.taylor = taylor
        self.erase_dpol_offdiag = erase_dpol_offdiag

        # what is already computed
        self.monomers_computed            = False   # unperturbed monomenr wavefunctions at arbitrary level of theory 
        self.densities_computed           = False   # density matrices necessary in full basis set
        self.energy_coulomb_computed      = False   # partitioning of 1st-order electrostatic energy
        self.energy_pauli_computed        = False   # partitioning of Pauli exchange-repulsion energy
        self.energy_polar_computed        = False   # partitioning of polarization energy (full QM)
        self.energy_polar_approx_computed = False   # partitioning of polarization energy (approximated)
        self.energy_ind_computed          = False   # induction part of polarization energy                                  
        self.energy_disp_computed         = False   # dispersion part of polarization energy
        self.energy_ct_compute            = False   # charge-transfer part of polarization energy
        self.energy_full_QM_computed      = False   # full QM total interaction energy
        self.dms_ind_computed             = False   # density matrix polarization susceptibility tensors for induction
        self.dms_disp_computed            = False   # density matrix polarization susceptibility tensors for dispersion
        self.dms_ct_computed              = False   # density matrix polarization susceptibility tensors for charge-transfer

        # basis set and JK object for entire aggregate
        self.bfs = psi4.core.BasisSet.build(aggregate, "BASIS", psi4.core.get_global_option("BASIS"), **kwargs)
        self.global_jk = psi4.core.JK.build(self.bfs, jk_type=jk_type)
        self.global_jk.set_memory(int(5e8))
        self.global_jk.initialize()
        Density.__init__(self, None, self.global_jk)

        # sizing
        self.nmo_t      = None
        self.nbf_t      = self.bfs.nbf()

        # sanity checks
        assert(self.aggregate.nfragments() >= 2), "Aggregate needs to have at least two fragments!"
        assert(not (self.method.lower()=='hf' and xc_scale!=1.0)), "Scaling parameter for XC density has to be 1.0 in HF theory!"
        None


    # ---- public interface ---- #

    def compute(self, polar_approx = True): 
        """\
 Perform the full density and interaction energy decompositions.

 Options: 
  o polar_approx - in addition to exact polarization energy, compute 
                   also the approximated polarization energy using
                   NO-expansion of exchange-correlation 2-electron density
                   and exact Pauli, polarization and unperturbed 1-electron densities.
                   Default: True

 Notes:
  o Exact polarization energy is always computed as a difference between
    the full QM interaction energy and all the remaining energies (Coulombic,
    exchange and repulsion energies).
"""
        # compute monomer wavefunctions
        if not self.monomers_computed          : self.compute_monomers()
        # compute full QM
        if not self.energy_full_QM_computed    : self.compute_full_QM()
        # compute densities
        if not self.densities_computed         : self.compute_densities()
        # compute interaction energies
        if not self.energy_coulomb_computed    : self.compute_coulomb()
        if not self.energy_pauli_computed      : self.compute_pauli()
        if not self.energy_polar_computed      : self.compute_polar()
        if polar_approx:
            if not self.energy_polar_approx_computed  : self.compute_polar_approx()

    def deformation_density(self, name):
        """\
 Compute the deformation 1-particle density matrix. 
 Returns density matrix in AO basis of entire molecular aggregate.

 Possible <name> entries:
  o fqm      - full QM deformation density
  o pau      - Pauli-repulsion denformation density
  o pol      - polarization deformation density
  o ind      - induction part of polarization deformation density
  o dis      - dispersion part of polarization deformation density
  o ct       - charge-transfer part of polarization deformation density
"""
        Doo, D = self._deformation_density_pauli()

        if   name.lower().startswith("fqm"):
             dD = self.matrix["dqm"] - D
        elif name.lower().startswith("pau"): 
             dD = Doo - D
        elif name.lower().startswith("pol"):
             dD = self.matrix["dqm"] - Doo
        elif name.lower().startswith("ind"):
             raise NotImplementedError
        elif name.lower().startswith("dis"):
             raise NotImplementedError
        elif name.lower().startswith("ct"):
             raise NotImplementedError
        elif name.lower().startswith("d0"):
             dD = D
        else: 
             raise ValueError("Incorrect name of density. Only pau, pol, ind, dis, ct or fqm are available.")
        return dD


    # ---- public interface (expert) ---- #

    def compute_monomers(self):
        "Compute all isolated (unperturbed) monomer wavefunctions"
        self._compute_monomers()
        self.monomers_computed = True

    def compute_full_QM(self):
        "Compute full QM interaction energy"
        assert self.monomers_computed is True

        self.aggregate.activate_all_fragments()
        c_ene, c_wfn = psi4.properties(self.method, molecule=self.aggregate, return_wfn=True, 
                                       properties=["DIPOLE","NO_OCCUPATIONS"], **self.kwargs)
        if self.method.lower().startswith('cc') and self.cc_relax is True:
           psi4.core.set_global_option("DERTYPE", "FIRST")
           psi4.core.set_global_option("WFN", self.method.upper())
           psi4.core.ccdensity(c_wfn)

        D    = c_wfn.Da       (           ).to_array(dense=True)
        S    = c_wfn.S        (           ).to_array(dense=True)
        C_hf = c_wfn.Ca_subset("AO", "ALL").to_array(dense=True)
        N, C = self.natural_orbitals(D, S=S, C=C_hf,
                                        orthogonalize_mo = True,
                                        order='descending', no_cutoff=0.0,
                                        return_ao_orthogonal = False,
                                        ignore_large_n = False,
                                        renormalize = False)
        if self.verbose is True: print("Sum of natural orbital occupations for nqm= %13.6f" % N.sum())
        self.matrix["dqm"] = D
        self.matrix["nqm"] = N
        self.matrix["cqm"] = C
        self.matrix["sqm"] = S
        self.matrix["chf"] = C_hf

        e_fqm_t = c_ene
        for i in range(self.aggregate.nfragments()):
            e_fqm_t -= self.data["ene"][i]
        self.vars["e_fqm_t"] = e_fqm_t

        self.energy_full_QM_computed = True

        # XC scaling factor test
        xc_scale_test = math.sqrt(self._f_test(N))
        A = 1.00 # 9.01266
        B = 0.00 #-7.03388
        xc_scale_test = A * xc_scale_test + B
        print(" Model XC Scale = %15.6f" % xc_scale_test)
        if self.xc_scale_test: self.xc_scale = xc_scale_test
        return

    def compute_densities(self):
        "Compute all the necessary density matrices"
        assert self.monomers_computed is True

        S_mo_t    = numpy.zeros((self.nmo_t, self.nmo_t), numpy.float64)
        W_mo_t    = numpy.zeros((self.nmo_t            ), numpy.float64)
        C_ao_mo_t = numpy.zeros((self.nbf_t, self.nmo_t), numpy.float64)

        for i in range(self.aggregate.nfragments()):
            wfn_i = self.data["wfn"][i]
            non_i = self.data["non"][i]
            noc_i = self.data["noc"][i]
            nmo_i = self.data["nmo"][i]
            nbf_i = self.data["nbf"][i]
            ofm_i = self.data["ofm"][i]
            ofb_i = self.data["ofb"][i]
            S_ao_ii = wfn_i.S().to_array(dense=True)
            S_mo_ii = self.triplet(noc_i.T, S_ao_ii, noc_i)
            S_mo_t[ofm_i:ofm_i+nmo_i, ofm_i:ofm_i+nmo_i] = S_mo_ii
            W_mo_t[ofm_i:ofm_i+nmo_i] = numpy.sqrt(non_i)
            C_ao_mo_t[ofb_i:ofb_i+nbf_i,ofm_i:ofm_i+nmo_i] = noc_i

            for j in range(i):
                wfn_j = self.data["wfn"][j]
                noc_j = self.data["noc"][j]
                nmo_j = self.data["nmo"][j]
                nbf_j = self.data["nbf"][j]
                ofm_j = self.data["ofm"][j]
                ofb_j = self.data["ofb"][j]
                S_ao_ij = self._compute_overlap_ao(wfn_i, wfn_j)
                S_mo_ij = self.triplet(noc_i.T, S_ao_ij, noc_j)
                S_mo_t[ofm_i:ofm_i+nmo_i, ofm_j:ofm_j+nmo_j] = S_mo_ij 
                S_mo_t[ofm_j:ofm_j+nmo_j, ofm_i:ofm_i+nmo_i] = S_mo_ij.T

        n_mo_t = W_mo_t * W_mo_t

        WSW = self.triplet(numpy.diag(W_mo_t), S_mo_t, numpy.diag(W_mo_t))      # mo::mo
        WSWm12 = self.matrix_power(WSW, -0.5)                                   # mo::mo

        if 0:
            # test of WSW-1/2 Taylor expansion
            delta = S_mo_t.copy()
            for i in range(delta.shape[0]): delta[i,i] = 0.0
            delta2 = self.doublet(delta, delta)
            delta3 = self.doublet(delta, delta2)
            Wm12 = self.matrix_power(numpy.diag(W_mo_t), -0.5)
            Wm1  = self.matrix_power(numpy.diag(W_mo_t), -1.0)
           #print(Wm1.diagonal())
           #print(W_mo_t)
            WSWm12_a0= Wm1.copy()
            WSWm12_a1= Wm1.copy()
            WSWm12_a1-= self.triplet(Wm12, delta , Wm12) * 1.0/2.0
            WSWm12_a2 = WSWm12_a1.copy()
            #WSWm12_a2+=(self.doublet(Wm1, delta) + self.doublet(delta, Wm1)) * 1.0/4.0
            WSWm12_a2+= self.triplet(Wm12, delta2, Wm12) * 3.0/8.0
            WSWm12_a3 = WSWm12_a2.copy()
            WSWm12_a3-= self.triplet(Wm12, delta3, Wm12) * 5.0/16.0
    
            def rms(m1, m2): return math.sqrt(((m1-m2) * (m1-m2)).sum())
    
           #print()
           #print(WSWm12.diagonal())
           #print(WSWm12_a0.diagonal())
           #print(WSWm12_a1.diagonal())
           #print(WSWm12_a2.diagonal())
           #print(WSWm12_a3.diagonal())
    
            #print(WSWm12   [0,1])
            #print(WSWm12_a1[0,1])
            #print(WSWm12_a2[0,1])
    
            print(rms(WSWm12, WSWm12_a0))
            print(rms(WSWm12, WSWm12_a1))
            print(rms(WSWm12, WSWm12_a2))
            print(rms(WSWm12, WSWm12_a3))
            if self.taylor is not False:
                if   self.taylor == 0: WSWm12 = WSWm12_a0
                elif self.taylor == 1: WSWm12 = WSWm12_a1
                elif self.taylor == 2: WSWm12 = WSWm12_a2
                elif self.taylor == 3: WSWm12 = WSWm12_a3
                else: raise NotImplementedError

        K = self.triplet(WSWm12, numpy.diag(n_mo_t), WSWm12)                    # mo::mo
        CW  = self.doublet(C_ao_mo_t, numpy.diag(W_mo_t))                       # ao::mo
        Doo_ao_t = self.triplet(CW, K, CW.T)                                    # ao::ao
        D_ao_t = self.triplet(C_ao_mo_t, numpy.diag(n_mo_t), C_ao_mo_t.T)       # ao::ao


        if self.taylor is not False:
           delta = S_mo_t.copy()
           for i in range(delta.shape[0]): delta[i,i] = 0.0
           delta2 = self.doublet(delta, delta)
           delta3 = self.doublet(delta, delta2)

           K0 = numpy.diag(n_mo_t)
           K1 = K0 - self.triplet(numpy.diag(W_mo_t), delta , numpy.diag(W_mo_t))
           K2 = K1 + self.triplet(numpy.diag(W_mo_t), delta2, numpy.diag(W_mo_t))
           K3 = K2 - self.triplet(numpy.diag(W_mo_t), delta3, numpy.diag(W_mo_t))

           Kr = self.triplet(numpy.diag(W_mo_t), K, numpy.diag(W_mo_t))
           def rms(m1, m2): return math.sqrt(((m1-m2) * (m1-m2)).sum())

           print(" RMS Taylor = 0: %14.6f" % rms(K0, Kr))
           print(" RMS Taylor = 1: %14.6f" % rms(K1, Kr))
           print(" RMS Taylor = 2: %14.6f" % rms(K2, Kr))
           print(" RMS Taylor = 3: %14.6f" % rms(K3, Kr))

           if   self.taylor == 0: Doo_ao_t = self.triplet(C_ao_mo_t, K0, C_ao_mo_t.T)
           elif self.taylor == 1: Doo_ao_t = self.triplet(C_ao_mo_t, K1, C_ao_mo_t.T)
           elif self.taylor == 2: Doo_ao_t = self.triplet(C_ao_mo_t, K2, C_ao_mo_t.T)
           elif self.taylor == 3: Doo_ao_t = self.triplet(C_ao_mo_t, K3, C_ao_mo_t.T)
           else: raise NotImplementedError

        # NO analysis of Antisymmetrized wavefunction
        noo, coo = self.natural_orbitals(Doo_ao_t.copy(), S=self.matrix["sqm"], C=self.matrix["chf"],
                                         orthogonalize_mo = True,
                                         order='descending', no_cutoff=0.0,
                                         return_ao_orthogonal = False,
                                         ignore_large_n = False,
                                         renormalize = False)
        if self.verbose is True: print("Sum of natural orbital occupations for noo= %13.6f" % noo.sum())

        # save
        self.matrix["n"  ] = n_mo_t
        self.matrix["c"  ] = C_ao_mo_t
        self.matrix["d"  ] = D_ao_t
        self.matrix["noo"] = noo
        self.matrix["coo"] = coo
        self.matrix["doo"] = Doo_ao_t

        # set up state of the object
        self.densities_computed = True


        # test exchange-polarization commutability hypothesis
        if 0:
           def rms(m1, m2): return math.sqrt(((m1-m2) * (m1-m2)).sum())                                                                             
           doo, do = self._deformation_density_pauli()
           d       = self.matrix['dqm']
           delta_dpol = self.deformation_density('pol')
                                                                                                                                                    
           do_pol = do + delta_dpol
           n, c = self.natural_orbitals(do_pol.copy(), S=self.matrix["sqm"], C=self.matrix["chf"],
                                         orthogonalize_mo = True,
                                         order='descending', no_cutoff=0.0,
                                         return_ao_orthogonal = False,
                                         ignore_large_n = True,
                                         renormalize = True)
           w = numpy.sqrt(n)
           print(n)
           print("Sum of natural orbital occupations for ntest= %13.6f" % n.sum())
                                                                                                                                                    
           s = self.triplet(c.T, self.matrix["sqm"], c)
           wsw = self.triplet(numpy.diag(w), s, numpy.diag(w))
           wswm12 = self.matrix_power(wsw, -0.5)
           k = self.triplet(wswm12, numpy.diag(n), wswm12)
           cw= self.doublet(c, numpy.diag(w))
           d_test = self.triplet(cw, k, cw.T)
                                                                                                                                                    
           n_el_test = 2.0 * self.doublet(d_test, self.matrix["sqm"]).trace()
           n_el_o    = 2.0 * self.doublet(do    , self.matrix["sqm"]).trace()
           n_el_qm   = 2.0 * self.doublet(d     , self.matrix["sqm"]).trace()
           print(" N-el FQM = %13.6f" % n_el_qm)
           print(" N-el D0  = %13.6f" % n_el_o )
           print(" N-el test= %13.6f" % n_el_test )
                                                                                                                                                    
           print(" RMS dqm-do= %13.6f" % rms(d, do))
           print(" RMS dqm-dt= %13.6f" % rms(d, d_test))
           print(" RMS dqm-doo=%13.6f" % rms(d, self.matrix["doo"]))
                                                                                                                                                    
           # compute polarized 0-th order density for 1st monomer
           print(" Computing SCF of first monomer in presence of second one")
           s1 = SCF(self.data["wfn"][0])
                                                                                                                                                    
           # --- V nuc
           V1_nuc = self._V_n(1)
           # --- J
           n_1 = self._block_fragment_mo_matrix(numpy.diag(self.matrix["n"]), keep_frags=[1], off_diagonal=False)
           c_1 = self.matrix["c"]
           D_1 = self.triplet(c_1, n_1, c_1.T)
                                                                                                                                                    
           self.global_jk.C_clear()                                           
           self.global_jk.C_left_add(psi4.core.Matrix.from_array(D_1, ""))
           I = numpy.identity(D_1.shape[0], numpy.float64)
           self.global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
           self.global_jk.compute()
           J = self.global_jk.J()[0].to_array(dense=True)
           K = self.global_jk.K()[0].to_array(dense=True)
                                                                                                                                                    
           V1 = self._V_n(1) + 2.0 * J 
                                                                                                                                                    
           s1.run(V_ext=V1, guess=self.data["wfn"][0].Fa(), verbose=True, maxit=100)
           print(self.data["ene"][0], s1.E)
                                                                                                                                                    
           print(" Computing SCF of second monomer in presence of first one")
           s2 = SCF(self.data["wfn"][1])
                                                                                                                                                    
           # --- V nuc
           V2_nuc = self._V_n(0)
           # --- J
           n_2 = self._block_fragment_mo_matrix(numpy.diag(self.matrix["n"]), keep_frags=[0], off_diagonal=False)
           c_2 = self.matrix["c"]
           D_2 = self.triplet(c_2, n_2, c_2.T)
                                                                                                                                                    
           self.global_jk.C_clear()                                           
           self.global_jk.C_left_add(psi4.core.Matrix.from_array(D_2, ""))
           I = numpy.identity(D_2.shape[0], numpy.float64)
           self.global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
           self.global_jk.compute()
           J = self.global_jk.J()[0].to_array(dense=True)
           K = self.global_jk.K()[0].to_array(dense=True)
                                                                                                                                                    
           V2 = self._V_n(0) + 2.0 * J 
                                                                                                                                                    
           s2.run(V_ext=V2, guess=self.data["wfn"][1].Fa(), verbose=True, maxit=100)
           print(self.data["ene"][1], s2.E)
                                                                                                                                                    
           n1, c1 = self.natural_orbitals(s1.D.copy(), orthogonalize_first=self.matrix["sqm"], order='descending', no_cutoff=0.0)
           n2, c2 = self.natural_orbitals(s2.D.copy(), orthogonalize_first=self.matrix["sqm"], order='descending', no_cutoff=0.0)
                                                                                                                                                    
           # density matrix
           dp = numpy.zeros(self.matrix["sqm"].shape, numpy.float64)
           for i in range(n1.size): dp += n1[i] * numpy.outer(c1[:,i], c1[:,i])
           for i in range(n2.size): dp += n2[i] * numpy.outer(c2[:,i], c2[:,i])
                                                                                                                                                    
           n_el_dp = 2.0 * self.doublet(self.matrix["sqm"], dp).trace()
           print(" N-el dp  = %13.6f" % n_el_dp )
           print(" RMS dqm-dp= %13.6f" % rms(d, dp))
                                                                                                                                                    
           #n, c = self.natural_orbitals(dp.copy(), orthogonalize_first=self.matrix["sqm"], order='descending', no_cutoff=0.0)
           n = numpy.hstack((n1, n2))
           c = numpy.hstack((c1, c2))
           w = numpy.sqrt(n)
                                                                                                                                                    
           # compute test dp matrix
           s = self.triplet(c.T, self.matrix["sqm"], c)
           wsw = self.triplet(numpy.diag(w), s, numpy.diag(w))
           wswm12 = self.matrix_power(wsw, -0.5)
           k = self.triplet(wswm12, numpy.diag(n), wswm12)
           cw= self.doublet(c, numpy.diag(w))
           dp_test = self.triplet(cw, k, cw.T)
           print(" RMS dqm-dt= %13.6f" % rms(d, dp_test))
        return 

    def compute_coulomb(self):
        "Compute coulombic energy and its components"
        assert self.monomers_computed is True
        assert self.densities_computed is True
        # orthogonalized and non-orthogonalized OPDM's

        # core Hamiltonian
        mints = psi4.core.MintsHelper(self.bfs)
        H = mints.ao_potential()
        T = mints.ao_kinetic()
        H.add(T)
        del mints

        e_cou_n = e_cou_1 = e_cou_2 = e_cou_t = 0.0

        # nuclear repulsion energy
        for i in range(self.aggregate.nfragments()):
            for j in range(i):
                e_cou_n += self._nuclear_repulsion_energy(i, j)

        # electronic contribution
        for i in range(self.aggregate.nfragments()):
            V_i = self._V_n(i)
            if self.acbs is False:
                D_i = self._block_fragment_ao_matrix(self.matrix["d"], keep_frags=[i], off_diagonal=False)
            else:
               n_i = self._block_fragment_mo_matrix(numpy.diag(self.matrix["n"]), keep_frags=[i], off_diagonal=False)
               c_i = self.matrix["c"]
               D_i = self.triplet(c_i, n_i, c_i.T)
              
            for j in range(i):
               V_j = self._V_n(j)
               if self.acbs is False:
                  D_j = self._block_fragment_ao_matrix(self.matrix["d"], keep_frags=[j], off_diagonal=False)
               else:
                  n_j = self._block_fragment_mo_matrix(numpy.diag(self.matrix["n"]), keep_frags=[j], off_diagonal=False)
                  c_j = self.matrix["c"]
                  D_j = self.triplet(c_j, n_j, c_j.T)

               e_cou_1 += 2.0 * self.compute_1el_energy(D_i, V_j)
               e_cou_1 += 2.0 * self.compute_1el_energy(D_j, V_i)
               e_cou_2 += 4.0 * self.compute_2el_energy(D_i, D_j, type='j')
        e_cou_t = e_cou_n + e_cou_1 + e_cou_2

        # save
        self.vars["e_cou_n"] = e_cou_n
        self.vars["e_cou_1"] = e_cou_1
        self.vars["e_cou_2"] = e_cou_2
        self.vars["e_cou_t"] = e_cou_t
        self.vars["e_dds_t"]+= e_cou_t
        if self.l_dds:
           self.vars["e_dds_l_t"]+= e_cou_t

        self.energy_coulomb_computed = True
        return

    def compute_pauli(self):
        "Compute repulsion energy and its components"
        assert self.monomers_computed is True
        assert self.densities_computed is True
        # orthogonalized and non-orthogonalized OPDM's
        Doo, D = self._deformation_density_pauli()
        # deformation density matrix
        dD = Doo - D
        # core Hamiltonian
        mints = psi4.core.MintsHelper(self.bfs)
        H = mints.ao_potential()
        T = mints.ao_kinetic()
        H.add(T)
        del mints

        # ------- one-electron part
        e_rep_1 = 2.0 * self.compute_1el_energy(dD, H.to_array(dense=True))
        # ------- two-electron part
        e_rep_2 = 2.0 * self.compute_2el_energy(dD, dD + 2.0 * D, type='j')

        if self.l_dds:
           self.vars["i_dpau_dpau"] = 2.0 * self.compute_2el_energy(dD, dD, type='j')
           self.vars["i_dpau_d"   ] = 4.0 * self.compute_2el_energy(dD,  D, type='j')

        # ---     Exchange energy
        Dunp = self._generalized_density_matrix(self.matrix["noo"], self.matrix["coo"])

        e_exc_t =-self.compute_2el_energy(Dunp, Dunp, type='k')

        for i in range(self.aggregate.nfragments()):
            w_i = numpy.sqrt(self.matrix["n"])
            #print(" Test: Molecule %d: %16.3f" % (i+1, self._f_test(w_i)))
            c_i = self.matrix["c"]
            scale = math.sqrt(self.xc_scale)
            #scale = math.sqrt(2.0*self._f_test(w_i))
            w_i*= scale
            w_i = self._block_fragment_mo_matrix(numpy.diag(w_i), keep_frags=[i], off_diagonal=False)
            D_i = self.triplet(c_i, w_i, c_i.T)
            #D_i = self._block_fragment_ao_matrix(D, keep_frags=[i], off_diagonal=False) -> acbs=False (old)
            e_exc_t += self.compute_2el_energy(D_i, D_i, type='k')

        e_rep_t = e_rep_1 + e_rep_2
        e_exr_t = e_exc_t + e_rep_t
        # save
        self.vars["e_rep_1"] = e_rep_1
        self.vars["e_rep_2"] = e_rep_2
        self.vars["e_rep_t"] = e_rep_t
        self.vars["e_exc_t"] = e_exc_t
        self.vars["e_exr_t"] = e_exr_t
        self.vars["e_dds_t"]+= e_exr_t
        if self.l_dds:
           self.vars["e_dds_l_t"] += e_exr_t
           self.vars["e_dds_l_t"] -= self.vars["i_dpau_dpau"]

        self.energy_pauli_computed = True
        return

    def compute_polar(self):
        "Compute exact polarization energy"
        assert self.monomers_computed is True
        assert self.energy_coulomb_computed is True

        self.vars["e_pol_t"] = self.vars["e_fqm_t"] - self.vars["e_cou_t"] - self.vars["e_exr_t"]
        return

    def compute_polar_approx(self):
        "Compute approximate polarization interaction energy"
        assert self.monomers_computed is True
        assert self.densities_computed is True
        assert self.energy_full_QM_computed is True

        # orthogonalized and non-orthogonalized OPDM's
        Doo, D = self._deformation_density_pauli()
        dD_pau = self.deformation_density("pau")
        dD_pol = self.deformation_density("pol")

        if self.erase_dpol_offdiag is True: 
           for i in range(self.aggregate.nfragments()):
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               for j in range(self.aggregate.nfragments()):
                   if i!=j:
                      ofb_j = self.data["ofb"][j]
                      nbf_j = self.data["nbf"][j]
                      dD_pol[ofb_i:ofb_i+nbf_i, ofb_j:ofb_j+nbf_j].fill(0.0)
           #n = int(dD_pol.shape[0]/2)                            
           #dD_pol[:n,n:].fill(0.0)
           #dD_pol[n:,:n].fill(0.0)

        # core Hamiltonian
        mints = psi4.core.MintsHelper(self.bfs)
        H = mints.ao_potential()
        T = mints.ao_kinetic()
        H.add(T)
        del mints

        e_pol_1 = 2.0 * self.compute_1el_energy(dD_pol, H.to_array(dense=True))
        e_pol_2 = 2.0 * self.compute_2el_energy(dD_pol, dD_pol + 2.0 * (dD_pau + D), type='j')

        if self.l_dds:
           self.vars["i_dpol_dpol"] = 2.0 * self.compute_2el_energy(dD_pol, dD_pol, type='j')
           self.vars["i_dpol_dpau"] = 4.0 * self.compute_2el_energy(dD_pol, dD_pau, type='j')
           self.vars["i_dpol_d"   ] = 4.0 * self.compute_2el_energy(dD_pol,  D    , type='j')
           self.vars["i_dpol_hcor"] = e_pol_1

        # exchange-polarization energy
        Dpol = self._generalized_density_matrix(self.matrix["nqm"], self.matrix["cqm"])
        Dunp = self._generalized_density_matrix(self.matrix["noo"], self.matrix["coo"])

        e_exp_a = self.compute_2el_energy(Dunp, Dunp, type='k')
        e_exp_a-= self.compute_2el_energy(Dpol, Dpol, type='k')

        e_pol_e = e_pol_1 + e_pol_2
        e_pol_a = e_pol_e + e_exp_a

        self.vars["e_pol_e"   ] = e_pol_e 
        self.vars["e_pol_1"   ] = e_pol_1
        self.vars["e_pol_2"   ] = e_pol_2
        self.vars["e_exp_a"   ] = e_exp_a
        self.vars["e_pol_a"   ] = e_pol_a

        self.vars["e_exp_t"   ] = self.vars["e_pol_t"] - e_pol_1 - e_pol_2
        self.vars["e_dds_t"   ]+= e_pol_a
        if self.l_dds:
           self.vars["e_dds_l_t"] += e_pol_a
           self.vars["e_dds_l_t"] -= self.vars["i_dpol_dpau"] + self.vars["i_dpol_dpol"]

        self.energy_polar_approx_computed = True
        return


    # ---- printers ---- #

    def __repr__(self):
        "Print final results"
        log = "\n\n"
        log+= " ===> Density-Based interaction energy decomposition <=== \n\n"

        if self.monomers_computed:
           log += "   %s/%s level of theory.\n"                  %(self.method.upper(), self.bfs.name())
           log += "   Total number of basis functions  = %d\n"   % self.nbf_t
           log += "   Total number of natural orbitals = %d\n"   % self.nmo_t
           log += "   Natural orbital threshold        = %.8f\n" % self.no_cutoff
           log += "\n\n"

        if self.energy_coulomb_computed:
           log += "   DDS Results\n"                                                                              
           log += " ---------------------------------------------------------------------------------------\n"
           log += "                                  [A.U.]                [kcal/mole]          [kJ/mole]  \n"
           log += " ---------------------------------------------------------------------------------------\n"
           log += "   Electrostatics         " + self._print_line(self.vars["e_cou_t"]) + "\n"
          #log += "     E-coul(nuclear)      " + self._print_line(self.vars["e_cou_n"]) + "\n"
          #log += "     E-coul(1)            " + self._print_line(self.vars["e_cou_1"]) + "\n"
          #log += "     E-coul(2)            " + self._print_line(self.vars["e_cou_2"]) + "\n"
           log += "\n"

        if self.energy_pauli_computed:
           log += "   Exchange-Repulsion     " + self._print_line(self.vars["e_exr_t"]) + "\n"
           log += "\n"
           log += "     E-exchange           " + self._print_line(self.vars["e_exc_t"]) + "\n"
           log += "     E-repul              " + self._print_line(self.vars["e_rep_t"]) + "\n"
           log += "       E-repul(1)         " + self._print_line(self.vars["e_rep_1"]) + "\n"
           log += "       E-repul(2)         " + self._print_line(self.vars["e_rep_2"]) + "\n"
           if self.l_dds:
            log+= "         E-dpau-d         " + self._print_line(self.vars["i_dpau_d"   ]) + "\n"
            log+= "         E-dpau-dpau      " + self._print_line(self.vars["i_dpau_dpau"]) + "\n"
           log += "\n"


        if self.energy_full_QM_computed:
           log += "   Polarization           " + self._print_line(self.vars["e_pol_t"]) + "\n"
           log += "\n"
        if self.energy_polar_approx_computed:
           log += "     E-polar(el)          " + self._print_line(self.vars["e_pol_e"]) + "\n"
           log += "       E-polar(1)         " + self._print_line(self.vars["e_pol_1"]) + "\n"
           log += "       E-polar(2)         " + self._print_line(self.vars["e_pol_2"]) + "\n"
           if self.l_dds:
            log+= "         E-dpol-hcore     " + self._print_line(self.vars["i_dpol_hcor"]) + "\n"
            log+= "         E-dpol-d         " + self._print_line(self.vars["i_dpol_d"   ]) + "\n"
            log+= "         E-dpol-dpol      " + self._print_line(self.vars["i_dpol_dpol"]) + "\n"
            log+= "         E-dpol-dpau      " + self._print_line(self.vars["i_dpol_dpau"]) + "\n"
           log += "     E-ex-pol             " + self._print_line(self.vars["e_exp_t"]) + "\n"
           log += "       E-ex-pol(no)       " + self._print_line(self.vars["e_exp_a"]) + "\n"
           log += "\n"
           log += "   Polarization (approx)  " + self._print_line(self.vars["e_pol_a"]) + "\n"
        if self.energy_full_QM_computed:
           log += "\n"

        if self.energy_full_QM_computed:
           log += "   Supramolecular Energy  " + self._print_line(self.vars["e_fqm_t"]) + "\n"
           log += "   DDS Energy             " + self._print_line(self.vars["e_dds_t"]) + "\n" 
           if self.l_dds:
            log +="   L-DDS Energy           " + self._print_line(self.vars["e_dds_l_t"]) + "\n" 
           log += "\n"
           err  = self.vars["e_dds_t"] - self.vars["e_fqm_t"]
           log += "     DDS Error            " + self._print_line(err) + "\n"
           if self.l_dds:
            err  = self.vars["e_dds_l_t"] - self.vars["e_fqm_t"]
            log +="     L-DDS Error          " + self._print_line(err) + "\n"
           log += " ---------------------------------------------------------------------------------------\n"
           log += "\n"

        return str(log)

    def print_out(self):
        "Print out to the Psi4 output file"
        psi4.core.print_out(str(self))



    # ---- public utilities ---- # 

    def doublet(self, A, B):
        "Compute double matrix product: result = AB"
        return numpy.dot(A, B)

    def triplet(self, A, B, C):
        "Compute triple matrix product: result = ABC"
        return numpy.dot(A, numpy.dot(B, C))

    def matrix_power(self, M, x, eps=1.0e-6):
        "Computes the well-behaved matrix power. All eigenvalues below or equal eps are ignored"
        E, U = numpy.linalg.eigh(M)
        Ex = E.copy(); Ex.fill(0.0)
        for i in range(len(E)):
            if (E[i]>0.0+eps) : Ex[i] = numpy.power(E[i], x)
            else: Ex[i] = 0.0
        Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
        return Mx

    def rms(self, m1, m2): 
        "RMS between two matrices"
        return math.sqrt(((m1-m2) * (m1-m2)).sum())                                                                             

    # --- protected interface --- #

    def _compute_monomers(self):
        "Compute monomer wavefunctions in aggregate-centred basis set"
        self.nmo_t = 0
        ofm_n = 0
        ofb_n = 0
        for n in range(self.aggregate.nfragments()):
            id_m = [n+1]
            if self.acbs is True: 
               id_g = list(range(1, self.aggregate.nfragments()+1))
               id_g.pop(n)
               current_molecule = self.aggregate.extract_subsets(id_m, id_g)
            else      : 
               current_molecule = self.aggregate.extract_subsets(id_m)

            c_ene, c_wfn = psi4.properties(self.method, molecule=current_molecule, return_wfn=True, 
                                           properties=['DIPOLE', 'NO_OCCUPATIONS'], **self.kwargs)

            if self.method.lower().startswith('cc') and self.cc_relax is True:
               psi4.core.set_global_option("DERTYPE", "FIRST")
               psi4.core.set_global_option("WFN", self.method.upper())
               psi4.core.ccdensity(c_wfn)

            D = c_wfn.Da().to_array(dense=True)
            N, C = self.natural_orbitals(D, S=c_wfn.S().to_array(dense=True), C=c_wfn.Ca_subset("AO","ALL").to_array(dense=True),
                                         orthogonalize_mo = True,
                                         order='descending', no_cutoff=0.0,
                                         return_ao_orthogonal = False,
                                         ignore_large_n = False,
                                         renormalize = False)
            if self.verbose is True: print("Sum of natural orbital occupations for n(%d)= %13.6f" % (n, N.sum()))

            self.data['wfn'].append(c_wfn)
            self.data['ene'].append(c_ene)
            self.data['odm'].append(D)
            self.data['non'].append(N)
            self.data['noc'].append(C)

            self.data['nmo'].append(C.shape[1])
            self.nmo_t += self.data['nmo'][n]
            self.data['nbf'].append(c_wfn.basisset().nbf()) 
            self.data['ofm'].append(ofm_n)
            self.data['ofb'].append(ofb_n)
            ofm_n += self.data['nmo'][n]
            if self.acbs is False: 
               ofb_n += self.data['nbf'][n]

            psi4.core.clean()
        return

    def _deformation_density_pauli(self):
        Doo = self.matrix["doo"]
        D   = self.matrix["d"]
        return Doo, D

    def _V(self, mol, bfs):
      "Potential matrix of mol's nuclei in bfs AO basis"
      ep = psi4.core.ExternalPotential()
      BohrToAngstrom = 0.5291772086
      for a in range(mol.natom()):
          q = mol.Z(a)
          x = mol.x(a) * BohrToAngstrom
          y = mol.y(a) * BohrToAngstrom
          z = mol.z(a) * BohrToAngstrom
          ep.addCharge(q, x, y, z)
      V = ep.computePotentialMatrix(bfs).to_array(dense=True)
      return V

    def _V_n(self, n):
        "Potential matrix due to n-th fragment nuclei in ACBS"
        id_m = [n+1]
        id_g = list(range(1, self.aggregate.nfragments()+1))
        id_g.pop(n)
        frag = self.aggregate.extract_subsets(id_m, id_g)
        V_n = self._V(frag, self.bfs)
        return V_n

    def _nuclear_repulsion_energy(self, i, j):
        "Compute nuclear repulsion energy between molecules i and j"
        mol_i = self.aggregate.extract_subsets(i+1)
        mol_j = self.aggregate.extract_subsets(j+1)
        energy = 0.0
        for ni in range(mol_i.natom()):
            xi = mol_i.x(ni)
            yi = mol_i.y(ni)
            zi = mol_i.z(ni)
            Zi = mol_i.Z(ni)
            for nj in range(mol_j.natom()):
                xj = mol_j.x(nj)
                yj = mol_j.y(nj)
                zj = mol_j.z(nj)
                Zj = mol_j.Z(nj)

                rij = math.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                energy += Zi*Zj/rij
        return energy





    # ---- protected utilities ---- # 

    def _generalized_density_matrix(self, n, c):
        "Compute occupation-weighted 1-electron density matrix in AO basis"
        ns = numpy.sqrt(n.real); N = c.shape[0]
        D  = numpy.zeros((N, N), numpy.float64)
        #scale = self.xc_scale ** (1.0 - ns)
        #print(" Test: %16.3f" % (self._f_test(n)))
        scale = math.sqrt(self.xc_scale)
        #print("Scale", scale)
        #scale = math.sqrt(self._f_test(n))
        #print("Next Scale", scale)
        ns *= scale
        for i in range(n.size):
            D += ns[i] * numpy.outer(c.real[:,i], c.real[:,i]) 
        return D

    def _compute_overlap_ao(self, wfn_i, wfn_j):
        "Compute AO overlap matrix between fragments"
        mints = psi4.core.MintsHelper(wfn_i.basisset())
        S_ij = mints.ao_overlap(wfn_i.basisset(), wfn_j.basisset()).to_array(dense=True)
        return S_ij

    def _block_fragment_ao_matrix(self, M, keep_frags=[], off_diagonal=False):
        "Returns the block-fragment form of matrix M (AO-basis) with fragments as each block. Keep only indicated fragments."
        M_frag = M.copy(); M_frag.fill(0.0)
        for i in range(self.aggregate.nfragments()):
            if i in keep_frags:
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               M_frag[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i] = numpy.array(M[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i], numpy.float64)
               # keep offdiagonals
               if off_diagonal is True:
                  raise NotImplementedError
        return M_frag

    def _block_fragment_mo_matrix(self, M, keep_frags=[], off_diagonal=False):
        "Returns the block-fragment form of matrix M (MO-basis) with fragments as each block. Keep only indicated fragments."
        M_frag = M.copy(); M_frag.fill(0.0)
        for i in range(self.aggregate.nfragments()):
            if i in keep_frags:
               ofm_i = self.data["ofm"][i]
               nmo_i = self.data["nmo"][i]
               M_frag[ofm_i:ofm_i+nmo_i, ofm_i:ofm_i+nmo_i] = numpy.array(M[ofm_i:ofm_i+nmo_i, ofm_i:ofm_i+nmo_i], numpy.float64)
               # keep offdiagonals
               if off_diagonal is True:
                  raise NotImplementedError
        return M_frag

    def _block_diagonal_ao_matrix(self, M, zero_frags=[]):
        "Returns the block-diagonal form of matrix M (AO-basis) with fragments as each block. Zero-out indicated fragments"
        M_diag = M.copy(); M_diag.fill(0.0)
        for i in range(self.aggregate.nfragments()):
            if not i in zero_frags:
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               M_diag[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i] = numpy.array(M[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i], numpy.float64)
        return M_diag

    def _print_line(self, value, type='e'):
        "Print one line of report in multiple units"
        line = ''
        # energy output
        # ---> [A.U.]  [kcal/mole]  [kJ/mole]
        if type == 'e':
           c1 = 627.5096080306   # A.U. to kcal/mole
           c2 = 2625.5002        # A.U. to kJ/mole
        # other output
        else: raise NotImplementedError

        v1 = value * c1
        v2 = value * c2
        line = "%18.8f   %18.8f   %18.8f" % (value, v1, v2)
        return line
    def _f_test(self, n):
        N = n.sum()
        _f1, _f2 = 0.0, 0.0
        for i in range(len(n)):
            for j in range(len(n)):
                _f1 += math.sqrt(n[i] * n[j])
                #_f2 += math.sqrt(n[i] * n[j])
        _f1 = N*N / _f1
        #_f2 = N / _f2
        return _f1


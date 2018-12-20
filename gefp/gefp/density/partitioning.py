#!/usr/bin/python3
#*-* coding: utf-8 *-*

import sys
import psi4
import oepdev
import numpy
import numpy.linalg

__all__ = ["DensityDecomposition"]

class DensityDecomposition:
    """
 Density-Based Decomposition Scheme from Mandado and Hermida-Ramon
 JCTC 2011
"""
    def __init__(self, aggregate, method='hf', ACBS=True, jk_type='direct', no_cutoff=0.000, **kwargs):
        "Initialize all attributes"
        # molecular aggregate
        self.aggregate   = aggregate    
        # level of theory
        self.method      = method       
        # wavefunction data of each unperturbed fragment
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
        # variables
        self.vars        = {"e_cou_1"     : None,   # coulombic energy, 1el part
                            "e_cou_2"     : None,   # coulombic energy, 2el part
                            "e_cou_t"     : None,   # coulombic energy
                            "e_rep_1"     : None,   # repulsion energy, 1el part
                            "e_rep_2"     : None,   # repulsion energy, 2el part
                            "e_rep_t"     : None,   # repulsion energy, (1+2)el part
                            "e_exc_t"     : None,   # exchange energy
                            "e_exr_t"     : None,   # exchange-repulsion energy
                            "e_pol_t"     : None,   # polarization energy
                            "e_pol_ind_t" : None,   # induction part of polarization energy, (1+2)el part
                            "e_pol_disp_t": None,   # dispersion part of polarization energy, (1+2)el part
                            "e_pol_ct_t"  : None,   # charge-transfer part of polarization energy, (1+2)el part
                            "e_t"         : None,   # total interaction energy
                            }
        # additional options
        self.kwargs      = kwargs
        # basis-set mode (ACBS: aggregate-centred basis set)
        self.acbs        = ACBS
        # cutoff for NO occupancies
        self.no_cutoff   = no_cutoff

        # what is already computed
        self.monomers_computed      = False   # unperturbed monomenr wavefunctions at arbitrary level of theory
        self.coulomb_computed       = False   # partitioning of 1st-order electrostatic energy
        self.pauli_computed         = False   # partitioning of Pauli exchange-repulsion energy
        self.polarization_computed  = False   # partitioning of polarization energy
        self.in_induction_computed  = False   # induction part of polarization energy
        self.in_dispersion_computed = False   # dispersion part of polarization energy
        self.in_ct_computed         = False   # charge-transfer part of polarization energy
        self.dms_ind_computed       = False   # density matrix polarization susceptibility tensors for induction
        self.dms_disp_computed      = False   # density matrix polarization susceptibility tensors for dispersion
        self.dms_ct_computed        = False   # density matrix polarization susceptibility tensors for charge-transfer

        # basis set and JK object for entire aggregate
        self.bfs = psi4.core.BasisSet.build(aggregate, "BASIS", psi4.core.get_global_option("BASIS"), puream=-1)
        self.global_jk = psi4.core.JK.build(self.bfs, jk_type=jk_type)
        self.global_jk.set_memory(int(5e8))
        self.global_jk.initialize()

        # sizing
        self.nmo_t      = None
        self.nbf_t      = self.bfs.nbf()

        # sanity checks
        #if self.acbs: 
        #   print (" Only monomer-centred AO basis set computations are implemented so far.")
        #   sys.exit(1)

    def compute_monomers(self):
        "Compute all isolated (unperturbed) monomer wavefunctions"
        if self.acbs: self._compute_monomers(do_acbs=True )  # aggregate-centred AO basis 
        else:         self._compute_monomers(do_acbs=False)  # monomer-centred AO basis
        self.monomers_computed = True

    def _compute_monomers(self, do_acbs):
        "Compute monomer wavefunctions in aggregate-centred basis set"
        self.nmo_t = 0
        ofm_n = 0
        ofb_n = 0
        for n in range(self.aggregate.nfragments()):
            id_m = [n+1]
            if do_acbs: 
               id_g = list(range(1, self.aggregate.nfragments()+1))
               id_g.pop(n)
               current_molecule = self.aggregate.extract_subsets(id_m, id_g)
            else      : 
               current_molecule = self.aggregate.extract_subsets(id_m)

            c_ene, c_wfn = psi4.properties(self.method, molecule=current_molecule, return_wfn=True, 
                                           properties=['DIPOLE', 'NO_OCCUPATIONS'])
            D = numpy.array(c_wfn.Da(), numpy.float64)
            N, C = self.natural_orbitals(D, orthogonalize_first=c_wfn.S(), order='descending')

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
            if not do_acbs: 
               ofb_n += self.data['nbf'][n]

            psi4.core.clean()
        return

    def compute_densities(self):
        "Compute all the density matrices necessary"
        raise NotImplementedError
        return

    def compute(self): 
        "Perform the full density decomposition"
        # compute monomer wavefunctions
        if not self.monomers_computed: self.compute_monomers()
        raise NotImplementedError

    def _deformation_density_pauli(self):
        Doo, D = self._deformation_density_pauli()
        return Doo, D
    #def _deformation_density_pauli_acbs(self): --> deprecate
    #    raise NotImplementedError
    def _deformation_density_pauli(self):
        "Compute Pauli deformation density"
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
            S_ao_ii = wfn_i.S()
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
        K = self.triplet(WSWm12, numpy.diag(n_mo_t), WSWm12)                    # mo::mo
        CW  = self.doublet(C_ao_mo_t, numpy.diag(W_mo_t))                       # ao::mo

        Doo_ao_t = self.triplet(CW, K, CW.T)                                    # ao::ao
        D_ao_t = self.triplet(C_ao_mo_t, numpy.diag(n_mo_t), C_ao_mo_t.T)       # ao::ao

        return Doo_ao_t, D_ao_t

    def compute_repulsion(self):
        "Compute repulsion energy and its components"
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

        # ---     Coulombic energy
        #e_coul_1 = 
        # ---     Repulsion energy
        # ------- one-electron part
        e_rep_1 = 2.0 * self.compute_1el_energy(dD, numpy.array(H))
        e_rep_2 = 2.0 * self.compute_2el_energy(dD, dD + 2.0 * D, type='j')
        e_rep_t  = e_rep_1 + e_rep_2
        # ---     Exchange energy
        e_exc_t =-1.0 * ( self.compute_2el_energy(Doo, Doo, type='k', k_diag=False) - \
                          self.compute_2el_energy(D  , D  , type='k', k_diag=True, d_diag=self.data["odm"])   )
        e_exr_t = e_exc_t + e_rep_t
        # save
        self.vars["e_rep_1"] = e_rep_1
        self.vars["e_rep_2"] = e_rep_2
        self.vars["e_rep_t"] = e_rep_t
        self.vars["e_exc_t"] = e_exc_t
        self.vars["e_exr_t"] = e_exr_t

        self.pauli_computed = True
        return

    def __repr__(self):
        log = " Density-Based interaction energy decomposition\n\n"
        if self.monomers_computed:
           log += " Monomers computed at %s/%s level of theory.\n" % (self.method.upper(), self.bfs.name())
        #if self.
        if self.pauli_computed:
           log += " @Sec: Pauli repulsion energy:\n"
           log += "       E_REP_1=%12.6f E_REP_2=%12.6f E_REP_T=%12.6f\n" % (self.vars["e_rep_1"],
                                                                             self.vars["e_rep_2"],
                                                                             self.vars["e_rep_t"])
           log += "       E_EXC  =%12.6f E_EXREP=%12.6f\n"                % (self.vars["e_exc_t"],
                                                                             self.vars["e_exr_t"])

        return str(log)

    def natural_orbitals(self, D, orthogonalize_first=None, order='descending', original_ao_mo=True):
        "Compute the Natural Orbitals from a given ODPM"
        if orthogonalize_first is not None:
           S = orthogonalize_first
           D_ = self._orthogonalize_OPDM(D, S)
        else:
           D_ = D
        n, U = numpy.linalg.eigh(D_)
        n[numpy.where(n<0.0)] = 0.0
        if original_ao_mo:
           assert orthogonalize_first is not None
           U = numpy.dot(self._orthogonalizer(orthogonalize_first), U)
        if self.no_cutoff != 0.0:
           ids = numpy.where(n>=self.no_cutoff)
           n = n[ids]
           U =(U.T[ids]).T
        if order=='ascending': 
           pass
        elif order=='descending':
           n = n[  ::-1]
           U = U[:,::-1]
        else: raise ValueError("Incorrect order of NO orbitals. Possible only ascending or descending.")
        return n, U

    def deformation_density(self, name):
        "Compute the deformation density matrix"
        if name.lower() == "pauli": 
           Doo, D = self._deformation_density_pauli()
        dD = Doo - D
        return dD

    def compute_1el_energy(self, D, Hcore):
        "Compute generalized 1-electron energy"
        energy = numpy.dot(D, Hcore).trace()
        return energy

    def compute_2el_energy(self, D_left, D_right, type='j', k_diag=False, d_diag=None):
        "Compute generalized 2-electron energy"
        assert self.global_jk is not None
        if type.lower() == 'k' and k_diag is True:
           assert d_diag is not None
           if self.acbs is False: JorK = self._diagonal_K(d_diag, None)
           else                 : JorK = self._diagonal_K(None, D_left)
        else:
           self.global_jk.C_clear()                                           
           self.global_jk.C_left_add(psi4.core.Matrix.from_array(D_left, ""))
           I = numpy.identity(D_left.shape[0], numpy.float64)
           self.global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
           self.global_jk.compute()
           if   type.lower() == 'j': JorK = numpy.array(self.global_jk.J()[0])
           elif type.lower() == 'k': JorK = numpy.array(self.global_jk.K()[0])
           else: raise ValueError("Incorrect type of JK matrix. Only J ro K allowed.")
        energy = numpy.dot(JorK, D_right).trace()
        return energy


    # ---- public utilities ---- # 

    def doublet(self, A, B):
        "Compute double matrix product: result = AB"
        return numpy.dot(A, B)

    def triplet(self, A, B, C):
        "Compute triple matrix product: result = ABC"
        return numpy.dot(A, numpy.dot(B, C))

    def matrix_power(self, M, x):
        "Computes the well-behaved matrix power"
        E, U = numpy.linalg.eigh(M)
        Ex = E.copy(); Ex.fill(0.0)
        for i in range(len(E)):
            if (E[i]>0.0) : Ex[i] = numpy.power(E[i], x)
            else: Ex[i] = 1.0
        Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
        return Mx


    # ---- protected utilities ---- # 

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

    def _compute_overlap_ao(self, wfn_i, wfn_j):
        "Compute AO overlap matrix between fragments"
        mints = psi4.core.MintsHelper(wfn_i.basisset())
        S = numpy.array(mints.ao_overlap(wfn_i.basisset(),
                                         wfn_j.basisset()))
        return S

    def _diagonal_K(self, D_list, D):
        "Compute K matrix from only fragment-diagonal contributions of density matrix in D_list"
        #assert self.acbs is False
        if self.acbs is False:
           assert D is None
           K = numpy.zeros((self.nbf_t, self.nbf_t), numpy.float64)
           for i in range(self.aggregate.nfragments()):                                  
               jk = psi4.core.JK.build(self.data["wfn"][i].basisset(), jk_type='direct')
               jk.set_memory(int(5e8))
               jk.initialize()
               jk.C_clear()
                                                                                         
               D_i = numpy.array(D_list[i], numpy.float64)
               jk.C_left_add(psi4.core.Matrix.from_array(D_i, ""))
               I = numpy.identity(D_i.shape[0], numpy.float64)
               jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
               jk.compute()
               K_i = numpy.array(jk.K()[0])
                                                                                         
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               K[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i] = K_i
        else:
           assert D_list is None
           K = numpy.zeros((self.nbf_t, self.nbf_t), numpy.float64)
           for i in range(self.aggregate.nfragments()):                                  
               self.global_jk.C_clear()
                                                                                         
               D_i = self._block_fragment_ao_matrix(D, keep_frags=[i], off_diagonal=False)
               self.global_jk.C_left_add(psi4.core.Matrix.from_array(D_i, ""))
               I = numpy.identity(D_i.shape[0], numpy.float64)
               self.global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
               self.global_jk.compute()
               K_i = numpy.array(self.global_jk.K()[0])
                                                                                         
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               K += K_i

        return K

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

    def _block_diagonal_ao_matrix(self, M, zero_frags=[]):
        "Returns the block-diagonal form of matrix M (AO-basis) with fragments as each block. Zero-out indicated fragments"
        M_diag = M.copy(); M_diag.fill(0.0)
        for i in range(self.aggregate.nfragments()):
            if not i in zero_frags:
               ofb_i = self.data["ofb"][i]
               nbf_i = self.data["nbf"][i]
               M_diag[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i] = numpy.array(M[ofb_i:ofb_i+nbf_i, ofb_i:ofb_i+nbf_i], numpy.float64)
        return M_diag


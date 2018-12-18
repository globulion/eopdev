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
    def __init__(self, aggregate, method='hf', ACBS=True, jk_type='direct', **kwargs):
        "Initialize all attributes"
        self.aggregate   = aggregate
        self.method      = method
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
        self.kwargs      = kwargs
        self.acbs        = ACBS
        self.monomers_computed = False

        self.bfs = psi4.core.BasisSet.build(aggregate, "BASIS", psi4.core.get_global_option("BASIS"), puream=-1)
        self.global_jk = psi4.core.JK.build(self.bfs, jk_type=jk_type)
        self.global_jk.set_memory(int(5e8))
        self.global_jk.initialize()

        # sanity checks
        if self.acbs: 
           print (" Only monomer-centred AO basis set computations are implemented so far.")
           sys.exit(1)

    def compute_monomers(self):
        "Compute all isolated (unperturbed) monomer wavefunctions"
        if self.acbs: self._compute_monomers(do_acbs=True )  # aggregate-centred AO basis 
        else:         self._compute_monomers(do_acbs=False)  # monomer-centred AO basis
        self.monomers_computed = True

    def _compute_monomers(self, do_acbs):
        "Compute monomer wavefunctions in aggregate-centred basis set"
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
            if not do_acbs: 
               self.data['nmo'].append(c_wfn.basisset().nbf())
               self.data['nbf'].append(c_wfn.basisset().nbf()) 
               self.data['ofm'].append(ofm_n)
               self.data['ofb'].append(ofb_n)
               ofm_n += self.data['nmo'][n]
               ofb_n += self.data['nbf'][n]

            psi4.core.clean()
        return

    #def _compute_monomers_mcbs(self):  ### -> deprecate
    #    "Compute monomer wavefunctions in monomer-centred basis set"
    #    raise NotImplementedError

    def compute(self): 
        "Perform the full density decomposition"
        # compute monomer wavefunctions
        if not self.monomers_computed: self.compute_monomers()
        raise NotImplementedError
    def compute_repulsion(self):
        "Compute repulsion energy and its components"
        # orthogonalized and non-orthogonalized OPDM's
        Doo, D = self._deformation_density_pauli()
        # deformation density matrix
        dD = Doo - D
        # core Hamiltonian
        mints = psi4.core.MintsHelper(self.bfs)
        V = mints.ao_potential()
        T = mints.ao_kinetic()
        #V.add(T)

        # ---     Repulsion energy
        # ------- one-electron part
        #e_rep_1 = 2.0 * self.compute_1el_energy(dD, numpy.array(V))
        #e_rep_2 = 2.0 * self.compute_2el_energy(dD, dD + 2.0 * dD, type='j')
        #e_rep   = e_rep_1 + e_rep_2
        #print(e_rep_1)
        #print(e_rep_2)
        #print(e_rep  )
        e_rep_1_k = 2.0 * self.compute_1el_energy(dD, numpy.array(T))
        e_rep_1_p = 2.0 * self.compute_1el_energy(dD, numpy.array(V))
        e_rep_1   = e_rep_1_k + e_rep_1_p
        e_rep_2_e = 4.0 * self.compute_2el_energy(dD, D, 'j')
        e_rep_2_p = 2.0 * self.compute_2el_energy(dD, dD, 'j')
        e_rep_2   = e_rep_2_e + e_rep_2_p
        e_rep     = e_rep_1 + e_rep_2
        print(e_rep_1_k, e_rep_1_p, e_rep_1)
        print(e_rep_2_e, e_rep_2_p, e_rep_2)
        print(e_rep)
        return

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
        if order=='ascending': return n, U
        else:                  return n[::-1], U[:,::-1]

    def matrix_power(self, M, x):
        "Computes the well-behaved matrix power"
        E, U = numpy.linalg.eigh(M)
        Ex = E.copy(); Ex.fill(0.0)
        for i in range(len(E)):
            if (E[i]>0.0) : Ex[i] = E[i]**x
            else: Ex[i] = 1.0
        Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
        return Mx
    def compute_1el_energy(self, D, Hcore):
        "Compute generalized 1-electron energy"
        energy = numpy.dot(D, Hcore).trace()
        return energy
    def compute_2el_energy(self, D_left, D_right, type='j'):
        "Compute generalized 2-electron energy"
        assert self.global_jk is not None
        self.global_jk.C_clear()
        self.global_jk.C_left_add(psi4.core.Matrix.from_array(D_left, ""))
        I = numpy.identity(D_left.shape[0], numpy.float64)
        self.global_jk.C_right_add(psi4.core.Matrix.from_array(I, ""))
        self.global_jk.compute()
        if   type.lower() == 'j': JorK = numpy.array(self.global_jk.J()[0])
        elif type.lower() == 'k': JorK = numpy.array(self.global_jk.K()[0])
        energy = numpy.dot(JorK, D_right).trace()
        return energy
    def deformation_density(self, name):
        "Compute the deformation density matrix"
        if name.lower() == "pauli": 
           Doo, D = self._deformation_density_pauli()
        dD = Doo - D
        return dD
    def _deformation_density_pauli(self):
        if self.acbs: Doo, D = self._deformation_density_pauli_acbs()
        else        : Doo, D = self._deformation_density_pauli_mcbs()
        return Doo, D
    def _deformation_density_pauli_acbs(self):
        raise NotImplementedError
    def _deformation_density_pauli_mcbs(self):
        "Compute Pauli deformation density"
        nmo_t = 0; nbf_t = 0
        for i in range(self.aggregate.nfragments()):
            nmo_t += self.data["nmo"][i]
            nbf_t += self.data["nbf"][i]
        S_mo_t    = numpy.zeros((nmo_t, nmo_t), numpy.float64)
        W_mo_t    = numpy.zeros((nmo_t       ), numpy.float64)
        C_ao_mo_t = numpy.zeros((nbf_t, nmo_t), numpy.float64)

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

        WSW = self.triplet(numpy.diag(W_mo_t), S_mo_t, numpy.diag(W_mo_t))
        CW  = self.doublet(C_ao_mo_t, numpy.diag(W_mo_t))
        WSWm12 = self.matrix_power(WSW, -0.5)
        K = self.triplet(WSWm12, numpy.diag(n_mo_t), WSWm12)

        Doo_ao_t = self.triplet(CW, K, CW.T)
        D_ao_t = self.triplet(C_ao_mo_t, numpy.diag(n_mo_t), C_ao_mo_t.T)
        #dD_ao_pauli = Doo_ao_t - D_ao_t

        return Doo_ao_t, D_ao_t

    # ---- public utilities ---- # 

    def doublet(self, A, B):
        "Compute double matrix product: result = AB"
        return numpy.dot(A, B)

    def triplet(self, A, B, C):
        "Compute triple matrix product: result = ABC"
        return numpy.dot(A, numpy.dot(B, C))

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
       

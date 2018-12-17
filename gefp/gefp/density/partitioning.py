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
    def __init__(self, aggregate, method='hf', ACBS=True, **kwargs):
        "Initialize all attributes"
        self.aggregate   = aggregate
        self.method      = method
        self.data        = {"wfn": [],
                            "ene": [],
                            "odm": [],
                            "non": [],
                            "noc": [],
                            "nmo": [],
                            "off": []}
        self.kwargs      = kwargs
        self.acbs        = ACBS
        self.monomers_computed = False
        if not self.acbs: 
           print (" Only aggregate-centred AO basis set computations are implemented so far.")
           sys.exit(1)

    def compute_monomers(self):
        "Compute all isolated (unperturbed) monomer wavefunctions in aggregate-centred AO basis set"
        if self.acbs: self._compute_monomers_acbs()
        else:         self._compute_monomers_mcbs()
        self.monomers_computed = True

    def _compute_monomers_acbs(self):
        "Compute monomer wavefunctions in aggregate-centred basis set"
        off_n = 0
        for n in range(self.aggregate.nfragments()):
            id_m = [n+1]
            id_g = list(range(1, self.aggregate.nfragments()+1))
            id_g.pop(n)
            current_molecule = self.aggregate.extract_subsets(id_m, id_g)
            c_ene, c_wfn = psi4.properties(self.method, molecule=current_molecule, return_wfn=True, 
                                           properties=['DIPOLE', 'NO_OCCUPATIONS'])
            D = numpy.array(c_wfn.Da(), numpy.float64)
            n, C = self.natural_orbitals(D, orthogonalize_first=None, order='descending')

            # unperturbed wavefunction
            self.data['wfn'].append(c_wfn)
            # unperturbed total energy
            self.data['ene'].append(c_ene)
            # unperturbed one-particle density matrix
            self.data['odm'].append(D)
            # unperturbed natural occupations
            self.data['non'].append(n)
            # unperturbed natural orbital wavefunction coefficients
            self.data['noc'].append(C)
            # number of molecular orbitals
            self.data['nmo'].append(wfn.basisset().nbf())
            # MO offset
            self.data['off'].append(off_n)

            off_n = self.data['nmo'][n]
            psi4.core.clean()
        return

    def _compute_monomers_mcbs(self):
        "Compute monomer wavefunctions in monomer-centred basis set"
        raise NotImplementedError

    def compute(self): 
        "Perform the full density decomposition"
        if not self.monomers_computed: self.compute_monomers()
        raise NotImplementedError

    def natural_orbitals(self, D, orthogonalize_first=None, order='descending'):
        """Compute the Natural Orbitals from a given ODPM"""
        if orthogonalize_first is not None:
           S = orthogonalize_first
           D_ = self._orthogonalize_OPDM(D, S)
        else:
           D_ = D
        n, U = numpy.linalg.eigh(D_)
        if order=='ascending': return n, U
        else:                  return n[::-1], U[:,::-1]

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
    def matrix_power(self, M, x):
        "Computes the well-behaved matrix power"
        E, U = numpy.linalg.eigh(m)           
        Ex = E.copy(); Ex.fill(0.0)
        for i in range(len(E)):
            if (E[i]>0.0) : Ex[i] = E[i]**x
            else: Ex[i] = 1.0
        Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
        return Mx
    def deformation_density(self, name):
        "Compute the deformation density matrix"
        if name.lower() == "pauli": dD = self._deformation_density_pauli()
        return dD
    def _deformation_density_pauli(self):
        "Compute Pauli deformation density"
        nmo_t = 0.0
        for i in range(self.aggregates.nfragments()):
            nmo_t += self.data["nmo"][i]
        S_mo_t = numpy.zeros((nmo_t, nmo_t), numpy.float64)
        W_mo_t = numpy.zeros((nmo_t), numpy.float64)

        for i in range(self.aggregates.nfragments()):
            wfn_i = self.data["wfn"][i]
            noc_i = self.data["noc"][i]
            nmo_i = self.data["nmo"][i]
            off_i = self.data["off"][i]
            S_ao_ii = wfn_i.S()
            S_mo_ii = self.triplet(noc_i.T, S_ao_ii, noc_i)
            S_mo_t[off_i:off_i+nmo_i, off_i:off_i+nmo_i] = S_mo_ii

            for j in range(i):
                wfn_j = self.data["wfn"][j]
                noc_j = self.data["noc"][j]
                nmo_j = self.data["nmo"][j]
                off_j = self.data["off"][j]
                S_ao_ij = self._compute_overlap_ao(wfn_i, wfn_j)
                S_mo_ij = self.triplet(noc_i.T, S_ao_ij, noc_j)
                S_mo_t[off_i:off_i+nmo_i, off_j:off_j+nmo_j] = S_mo_ij 
                S_mo_t[off_j:off_j+nmo_j, off_i:off_i+nmo_i] = S_mo_ij.T

        n_mo_t = W_mo_t * W_mo_t

        WSW = self.triplet(numpy.diag(W_mo_t), S_mo_t, numpy.diag(W_mo_t))
        CW  = self.doublet(C_t, numpy.diag(W_mo_t))
        WSWm12 = self.matrix_power(WSW, -0.5)
        K = self.triplet(WSWm12, numpy.diag(n_mo_t), WSWm12)

        Doo_t = self.triplet(CW, K, CW.T)
        D_t = self.triplet(C_t, numpy.diag(n_mo_t), C_t.T)
        dD_pauli = Doo_t - D_t

        return dD_pauli
    def _compute_overlap_ao(self, wfn_i, wfn_j):
        "Compute AO overlap matrix between fragments"
        raise NotImplementedError
        return S
    def doublet(self, A, B):
        "Compute double matrix product: result = AB"
        return numpy.dot(A, B)
    def triplet(self, A, B, C):
        "Compute triple matrix product: result = ABC"
        return numpy.dot(A, numpy.dot(B, C))

        

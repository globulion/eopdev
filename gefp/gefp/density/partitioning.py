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
                            "ene": []}
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
        for n in range(self.aggregate.nfragments()):
            id_m = [n+1]
            id_g = list(range(1, self.aggregate.nfragments()+1))
            id_g.pop(n)
            current_molecule = self.aggregate.extract_subsets(id_m, id_g)
            c_ene, c_wfn = psi4.properties(self.method, molecule=current_molecule, return_wfn=True, 
                                           properties=['DIPOLE', 'NO_OCCUPATIONS'])
            self.data['wfn'].append(c_wfn)
            self.data['ene'].append(c_ene)
            psi4.core.clean()
        return

    def _compute_monomers_mcbs(self):
        "Compute monomer wavefunctions in monomer-centred basis set"
        raise NotImplementedError

    def compute(self): 
        "Perform the full density decomposition"
        if not self.monomers_computed: self.compute_monomers()
        raise NotImplementedError

    def compute_no(self, D, orthogonalize_first=None, order='descending'):
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

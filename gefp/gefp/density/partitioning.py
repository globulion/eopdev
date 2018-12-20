#!/usr/bin/python3
#*-* coding: utf-8 *-*

import sys
import psi4
import oepdev
import math
import numpy
import numpy.linalg

__all__ = ["DensityDecomposition"]

class DensityDecomposition:
    """
 Density-Based Decomposition Scheme from Mandado and Hermida-Ramon
 JCTC 2011
 Usage:
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
        # matrices in aggregate AO/MO basis
        self.matrix = {"d"   : None,    # unperturbed 1-electron density
                       "c"   : None,    # unperturbed LCAO-NO coefficients
                       "n"   : None,    # unperturbed occupation numbers
                       "doo" : None,    # orthogonalized unpolarized 1-electron density
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
                            "e_t"         : 0.0 ,   # total interaction energy
                            "e_fqm_t"     : None,   # total interaction energy (full QM)
                            }
        # additional options
        self.kwargs      = kwargs
        # basis-set mode (ACBS: aggregate-centred basis set)
        self.acbs        = ACBS
        # cutoff for NO occupancies
        self.no_cutoff   = no_cutoff

        # what is already computed
        self.monomers_computed      = False   # unperturbed monomenr wavefunctions at arbitrary level of theory
        self.densities_computed     = False   # density matrices necessary in full basis set
        self.energy_coulomb_computed= False   # partitioning of 1st-order electrostatic energy
        self.energy_pauli_computed  = False   # partitioning of Pauli exchange-repulsion energy
        self.energy_polar_computed  = False   # partitioning of polarization energy
        self.energy_ind_computed    = False   # induction part of polarization energy
        self.energy_disp_computed   = False   # dispersion part of polarization energy
        self.energy_ct_compute      = False   # charge-transfer part of polarization energy
        self.energy_full_QM_computed= False   # full QM total interaction energy
        self.dms_ind_computed       = False   # density matrix polarization susceptibility tensors for induction
        self.dms_disp_computed      = False   # density matrix polarization susceptibility tensors for dispersion
        self.dms_ct_computed        = False   # density matrix polarization susceptibility tensors for charge-transfer

        # basis set and JK object for entire aggregate
        self.bfs = psi4.core.BasisSet.build(aggregate, "BASIS", psi4.core.get_global_option("BASIS"), puream=-1, **kwargs)
        self.global_jk = psi4.core.JK.build(self.bfs, jk_type=jk_type)
        self.global_jk.set_memory(int(5e8))
        self.global_jk.initialize()

        # sizing
        self.nmo_t      = None
        self.nbf_t      = self.bfs.nbf()

        # sanity checks
        None


    # ---- public interface ---- #

    def compute(self): 
        "Perform the full density decomposition"
        # compute monomer wavefunctions
        if not self.monomers_computed          : self.compute_monomers()
        # compute densities
        if not self.densities_computed         : self.compute_densities()
        # compute interaction energies
        if not self.energy_coulomb_computed    : self.compute_coulomb()
        if not self.energy_pauli_computed      : self.compute_pauli()
        if not self.energy_full_QM_computed    : self.compute_full_QM()

    def compute_monomers(self):
        "Compute all isolated (unperturbed) monomer wavefunctions"
        self._compute_monomers()
        self.monomers_computed = True

    def compute_densities(self):
        "Compute all the necessary density matrices"
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

        # save
        self.matrix["n"] = n_mo_t
        self.matrix["c"] = C_ao_mo_t
        self.matrix["d"] = D_ao_t
        self.matrix["doo"] = Doo_ao_t

        # set up state of the object
        self.densities_computed = True
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
               e_cou_2 += 4.0 * self.compute_2el_energy(D_i  , D_j  , type='j')
        e_cou_t = e_cou_n + e_cou_1 + e_cou_2

        # save
        self.vars["e_cou_n"] = e_cou_n
        self.vars["e_cou_1"] = e_cou_1
        self.vars["e_cou_2"] = e_cou_2
        self.vars["e_cou_t"] = e_cou_t
        self.vars["e_t"    ]+= e_cou_t

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
        e_rep_1 = 2.0 * self.compute_1el_energy(dD, numpy.array(H))
        e_rep_2 = 2.0 * self.compute_2el_energy(dD, dD + 2.0 * D, type='j')
        e_rep_t  = e_rep_1 + e_rep_2
        # ---     Exchange energy
        e_exc_t =-self.compute_2el_energy(Doo, Doo, type='k')
        for i in range(self.aggregate.nfragments()):
            if self.acbs is False:
               D_i = self._block_fragment_ao_matrix(D, keep_frags=[i], off_diagonal=False)
            else: 
               n_i = self._block_fragment_mo_matrix(numpy.diag(self.matrix["n"]), keep_frags=[i], off_diagonal=False)
               c_i = self.matrix["c"]
               D_i = self.triplet(c_i, n_i, c_i.T)
            e_exc_t += self.compute_2el_energy(D_i  , D_i  , type='k')
        e_exr_t = e_exc_t + e_rep_t
        # save
        self.vars["e_rep_1"] = e_rep_1
        self.vars["e_rep_2"] = e_rep_2
        self.vars["e_rep_t"] = e_rep_t
        self.vars["e_exc_t"] = e_exc_t
        self.vars["e_exr_t"] = e_exr_t
        self.vars["e_t"    ]+= e_exr_t

        self.energy_pauli_computed = True
        return

    def compute_full_QM(self):
        "Compute full QM interaction energy"
        assert self.monomers_computed is True

        self.aggregate.activate_all_fragments()
        c_ene, c_wfn = psi4.energy(self.method, molecule=self.aggregate, return_wfn=True, **self.kwargs)

        e_fqm_t = c_ene
        for i in range(self.aggregate.nfragments()):
            e_fqm_t -= self.data["ene"][i]
        self.vars["e_fqm_t"] = e_fqm_t
        self.vars["e_pol_t"] = e_fqm_t - self.vars["e_cou_t"] - self.vars["e_exr_t"]

        self.energy_full_QM_computed = True
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
        if   name.lower().startswith("pau"): 
             Doo, D = self._deformation_density_pauli()
             dD = Doo - D
        elif name.lower().startswith("pol"):
             raise NotImplementedError
        elif name.lower().startswith("ind"):
             raise NotImplementedError
        elif name.lower().startswith("dis"):
             raise NotImplementedError
        else: 
             raise ValueError("Incorrect name of density. Only pau, pol, ind or dis are available.")
        return dD

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
        else: raise ValueError("Incorrect type of JK matrix. Only J ro K allowed.")
        energy = numpy.dot(JorK, D_right).trace()
        return energy



    # ---- printers ---- #

    def __repr__(self):
        log = "\n\n"
        log+= " ===> Density-Based interaction energy decomposition <=== \n\n"
        if self.monomers_computed:
           log += " Monomers computed at %s/%s level of theory.\n" % (self.method.upper(), self.bfs.name())
           log += " Total number of basis functions  = %d\n" % self.nbf_t
           log += " Total number of natural orbitals = %d\n" % self.nmo_t
           log += "\n"
        if self.energy_coulomb_computed:
           log += " @Sec: Coulombic energy:\n"
           log += "       E_COU_1=%12.6f E_COU_2=%12.6f E_NUC_T=%12.6f\n" % (self.vars["e_cou_1"],
                                                                             self.vars["e_cou_2"],
                                                                             self.vars["e_cou_n"])
           log += "       E_COU_T=%12.6f\n"                               %  self.vars["e_cou_t"] 
        if self.energy_pauli_computed:
           log += " @Sec: Pauli repulsion energy:\n"
           log += "       E_REP_1=%12.6f E_REP_2=%12.6f E_REP_T=%12.6f\n" % (self.vars["e_rep_1"],
                                                                             self.vars["e_rep_2"],
                                                                             self.vars["e_rep_t"])
           log += "       E_EXC  =%12.6f E_EXREP=%12.6f\n"                % (self.vars["e_exc_t"],
                                                                             self.vars["e_exr_t"])
        if self.energy_full_QM_computed:
           log += " @Sec: Full QM energy:\n"
           log += "       E_POL_T=%12.6f\n"                                % self.vars["e_pol_t"]
           log += "       E_FQM_T=%12.6f\n"                                % self.vars["e_fqm_t"]
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

    def matrix_power(self, M, x):
        "Computes the well-behaved matrix power"
        E, U = numpy.linalg.eigh(M)
        Ex = E.copy(); Ex.fill(0.0)
        for i in range(len(E)):
            if (E[i]>0.0) : Ex[i] = numpy.power(E[i], x)
            else: Ex[i] = 1.0
        Mx = numpy.dot(U, numpy.dot(numpy.diag(Ex), U.T))
        return Mx


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
      V = numpy.array(ep.computePotentialMatrix(bfs), numpy.float64)
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




#!/usr/bin/python3
import psi4
import numpy

# helper methods
def oei(wfn_i, wfn_j, t='s'):
    "One-electron integrals"
    mints = psi4.core.MintsHelper(wfn_i.basisset())
    if   t == 's': O_ij = mints.ao_overlap(wfn_i.basisset(), wfn_j.basisset()).to_array(dense=True)
    elif t == 'v': O_ij = mints.ao_potential(wfn_i.basisset(), wfn_j.basisset()).to_array(dense=True)
    elif t == 't': O_ij = mints.ao_kinetic(wfn_i.basisset(), wfn_j.basisset()).to_array(dense=True)
    return O_ij
def transform(A, B, C):
    return numpy.linalg.multi_dot([A.T, B, C])
def get_V(mol, bfs):
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


def ct_energy(wfn_1, wfn_2):
    "CT energy from EFP2 method. Takes two RHF wavefunctions with canonical orbitals."
    # LCAO matrices
    C_occ_A = wfn_1.Ca_subset("AO","OCC").to_array(dense=True)
    C_occ_B = wfn_2.Ca_subset("AO","OCC").to_array(dense=True)
    C_vir_A = wfn_1.Ca_subset("AO","VIR").to_array(dense=True)
    C_vir_B = wfn_2.Ca_subset("AO","VIR").to_array(dense=True)
    # AO matrices
    F_ao_11 = wfn_1.Fa().to_array(dense=True)
    F_ao_22 = wfn_2.Fa().to_array(dense=True)

    S_ao_12 = oei(wfn_1, wfn_2, 's')
    T_ao_12 = oei(wfn_1, wfn_2, 't')
    VA_ao_12 = oei(wfn_1, wfn_2, 'v')
    VB_ao_12 = oei(wfn_2, wfn_1, 'v').T

    T_ao_11 = oei(wfn_1, wfn_1, 't')
    T_ao_22 = oei(wfn_2, wfn_2, 't')
    VA_ao_22 = get_V(wfn_1.molecule(), wfn_2.basisset())
    VB_ao_11 = get_V(wfn_2.molecule(), wfn_1.basisset())
    # MO matrices
    F_mo_11 = transform(C_occ_A, F_ao_11, C_occ_A)
    T_mo_11 = transform(C_vir_A, T_ao_11, C_vir_A)
    F_mo_22 = transform(C_occ_B, F_ao_22, C_occ_B)
    T_mo_22 = transform(C_vir_B, T_ao_22, C_vir_B)

    VB_mo_12 = transform(C_occ_A, VB_ao_12, C_vir_B)
    VB_mo_11o= transform(C_occ_A, VB_ao_11, C_occ_A)
    VB_mo_11v= transform(C_occ_A, VB_ao_11, C_vir_A)

    VA_mo_12 = transform(C_vir_A, VA_ao_12, C_occ_B)
    VA_mo_22o= transform(C_occ_B, VA_ao_22, C_occ_B)
    VA_mo_22v= transform(C_occ_B, VA_ao_22, C_vir_B)

    S_mo_12  = transform(C_occ_A, S_ao_12, C_occ_B)
    S_mo_1o2v= transform(C_occ_A, S_ao_12, C_vir_B)
    S_mo_1v2v= transform(C_vir_A, S_ao_12, C_vir_B)

    S_mo_2o1v = transform(C_occ_B, S_ao_12.T, C_vir_A)
    S_mo_2v1v = S_mo_1v2v.T

    T_mo_2v2o = transform(C_vir_B, T_ao_22, C_occ_B)
    T_mo_1o2o = transform(C_occ_A, T_ao_12, C_occ_B)
    T_mo_1v2o = transform(C_vir_A, T_ao_12, C_occ_B)

    T_mo_1v1o = transform(C_vir_A, T_ao_11, C_occ_A)
    T_mo_2o1o = transform(C_occ_B, T_ao_12.T, C_occ_A)
    T_mo_2v1o = transform(C_vir_B, T_ao_12.T, C_occ_A)


    # 2-electron contribution to potential matrices
    mints = psi4.core.MintsHelper(wfn_1)
    eri_1122 = mints.ao_eri(wfn_1.basisset(), wfn_1.basisset(), wfn_2.basisset(), wfn_2.basisset());
    eri_2211 = mints.ao_eri(wfn_2.basisset(), wfn_2.basisset(), wfn_1.basisset(), wfn_1.basisset());
    eri_1222 = mints.ao_eri(wfn_1.basisset(), wfn_2.basisset(), wfn_2.basisset(), wfn_2.basisset());
    eri_2111 = mints.ao_eri(wfn_2.basisset(), wfn_1.basisset(), wfn_1.basisset(), wfn_1.basisset());
    #print(wfn_1.basisset().nbf(), wfn_2.basisset().nbf())
    #print(eri_1122.shape, eri_1222.shape, eri_2111.shape)
    eri_O1_V2_O2_O2 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_1222, C_occ_A, C_vir_B, C_occ_B, C_occ_B)
    eri_O1_O1_O2_O2 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_1122, C_occ_A, C_occ_A, C_occ_B, C_occ_B)
    eri_O1_V1_O2_O2 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_1122, C_occ_A, C_vir_A, C_occ_B, C_occ_B)
    #
    eri_O2_V1_O1_O1 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_2111, C_occ_B, C_vir_A, C_occ_A, C_occ_A)
    eri_O2_O2_O1_O1 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_2211, C_occ_B, C_occ_B, C_occ_A, C_occ_A)
    eri_O2_V2_O1_O1 = numpy.einsum("ijkl,ia,jb,kc,ld->abcd", eri_2211, C_occ_B, C_vir_B, C_occ_A, C_occ_A)

    n_occ_1 = C_occ_A.shape[1]
    n_occ_2 = C_occ_B.shape[1]
    n_vir_1 = C_vir_A.shape[1]
    n_vir_2 = C_vir_B.shape[1]

    # VB_mo_12
    for i in range(n_occ_1):
        for n in range(n_vir_2):
            v = VB_mo_12[i, n]
            for j in range(n_occ_2): v += 2.0 * eri_O1_V2_O2_O2[i,n,j,j]
            VB_mo_12[i, n] = v
    # VB_mo_11o
    for i in range(n_occ_1):
        for k in range(n_occ_1):
            v = VB_mo_11o[i, k]
            for j in range(n_occ_2): v += 2.0 * eri_O1_O1_O2_O2[i,k,j,j]
            VB_mo_11o[i, k] = v
    # VB_mo_11v
    for i in range(n_occ_1):
        for m in range(n_vir_1):
            v = VB_mo_11v[i, m]
            for j in range(n_occ_2): v += 2.0 * eri_O1_V1_O2_O2[i,m,j,j]
            VB_mo_11v[i, m] = v

    # VA_mo_12
    for m in range(n_vir_1):
        for j in range(n_occ_2):
            v = VA_mo_12[m, j]
            for k in range(n_occ_1): v += 2.0 * eri_O2_V1_O1_O1[j,m,k,k]
            VA_mo_12[m, j] = v
    # VA_mo_22o
    for j in range(n_occ_2):
        for l in range(n_occ_2):
            v = VA_mo_22o[j, l]
            for k in range(n_occ_1): v += 2.0 * eri_O2_O2_O1_O1[j, l, k, k]
            VA_mo_22o[j, l] = v
    # VA_mo_22v
    for j in range(n_occ_2):
        for n in range(n_vir_2):
            v = VA_mo_22v[j, n]
            for k in range(n_occ_1): v += 2.0 * eri_O2_V2_O1_O1[j, n, k, k]
            VA_mo_22v[j, n] = v


    # Term A--->B

    # VB_I term
    VB_I_1o2v = VB_mo_12 - numpy.dot(VB_mo_11o, S_mo_1o2v) - numpy.dot(VB_mo_11v, S_mo_1v2v)

    # VB_2 term
    T_II_2v2o = T_mo_2v2o - numpy.dot(S_mo_1o2v.T, T_mo_1o2o) - numpy.dot(S_mo_1v2v.T, T_mo_1v2o)
    VB_II_1o2v = VB_I_1o2v + numpy.dot(T_II_2v2o, S_mo_12.T).T

    # normalization
    norm_2v = 1.0 / (1.0 - (S_mo_1o2v**2).sum(axis=0) - (S_mo_1v2v**2).sum(axis=0) )
    # denominator
    d_1o2v = 1.0 / ( F_mo_11.diagonal()[:,numpy.newaxis] - T_mo_22.diagonal()[numpy.newaxis,:] )

    E_AB = d_1o2v * VB_I_1o2v * VB_II_1o2v * norm_2v[numpy.newaxis, :]
    e_AB = 2.0 * E_AB.sum()


    # Term B--->A

    # VA_I term
    VA_I_2o1v = VA_mo_12.T - numpy.dot(VA_mo_22o, S_mo_2o1v) - numpy.dot(VA_mo_22v, S_mo_1v2v.T)

    # VB_2 term
    T_II_1v1o = T_mo_1v1o - numpy.dot(S_mo_2o1v.T, T_mo_2o1o) - numpy.dot(S_mo_1v2v, T_mo_2v1o)
    VA_II_2o1v = VA_I_2o1v + numpy.dot(S_mo_12.T, T_II_1v1o.T)

    # normalization
    norm_1v = 1.0 / (1.0 - (S_mo_2o1v**2).sum(axis=0) - (S_mo_2v1v**2).sum(axis=0) )
    # denominator
    d_2o1v = 1.0 / ( F_mo_22.diagonal()[:,numpy.newaxis] - T_mo_11.diagonal()[numpy.newaxis,:] )

    E_BA = d_2o1v * VA_I_2o1v * VA_II_2o1v * norm_1v[numpy.newaxis, :]
    e_BA = 2.0 * E_BA.sum()



    print(" A---->B energy: %14.6f" % e_AB)
    print(" B---->A energy: %14.6f" % e_BA)
    return e_AB + e_BA

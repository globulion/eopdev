#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 EFP2 Interaction Energy Solver Module.
 Bartosz Błasiak, Gundelfingen, 5 Nov 2019
"""

import os
import sys
import math
import numpy
import scipy.linalg
import scipy.optimize
import psi4
import oepdev
from abc import ABC, abstractmethod

__all__ = ["induction"]

def t_2(ra, rb):
    "Compute dipole-dipole interaction tensor"
    rab = ra - rb; r = numpy.linalg.norm(rab)
    T = numpy.outer(rab, rab) * 3.0 / r**5
    T-= numpy.identity(3)/r**3
    return T

def electric_field(cphf, r):
    "Compute electric field from CPHF object at position r"
    f = numpy.zeros(3)

    # prepare
    D = cphf.wfn().Da().to_array(dense=True)
    mints = psi4.core.MintsHelper(cphf.wfn().basisset())
    mol = cphf.wfn().molecule()

    ints = mints.electric_field(origin=r)

    # nuclear contribution
    rx, ry, rz = r
    for i in range(mol.natom()):
        x = mol.x(i)
        y = mol.y(i)
        z = mol.z(i)
        Z = numpy.float64(mol.Z(i))

        rai = numpy.array([rx-x, ry-y, rz-z])
        rai_norm = numpy.linalg.norm(rai)
        f += Z * rai / rai_norm**3

    # add electronic contribution
    f[0] += (D @ ints[0].to_array(dense=True)).trace() 
    f[1] += (D @ ints[1].to_array(dense=True)).trace() 
    f[2] += (D @ ints[2].to_array(dense=True)).trace() 

    return f



def induction(*cphfs):
    "Compute EFP2 induction energy from list of CPHF objects."
    n_frag = len(cphfs)
    # extract LMO centroids, distributed polarizabilities and other data
    lmoc = []
    dpol = []
    nmol = [x.ndocc() for x in cphfs]

    for cphf in cphfs:
        lmoc_cphf = [cphf.lmo_centroid(x) for x in range(cphf.ndocc())]
        lmoc.append(lmoc_cphf)
        dpol_cphf = [cphf.polarizability(x) for x in range(cphf.ndocc())]
        dpol.append(dpol_cphf)

    # allocate D and F matrices
    DIM = sum(DIM)
    D = numpy.zeros((3*DIM, 3*DIM))
    F = numpy.zeros( 3*DIM)

    # compute D and F matrices
    i_pol = 0
    for n in range(n_frag):
        N   = nmol[n]  # number of LMOc's
        An  = dpol[n]  # list of distributed polarizabilities
        Ln  = lmoc[n]  # list of LMOc's

        i_pol += N
        for a in range(N):  
            ix3 = 3*(i_pol - N) + 3*(a - 1)
            D[ix3:ix3+3,ix3:ix3+3] = numpy.linalg.inv(An[a].to_array(dense=True))

            ra = La[a].to_array(dense=True)

            j_pol = 0
            for m in range(n_frag):
              if (n!=m):
                M = nmol[m]
                LM = lmoc[m]

                rb = Lb[b].to_array(dense=True)

                j_pol += M
                for b in range(M):
                    jx3 = 3*(j_pol - M) + 3*(b - 1)
                    Tab = t_2(ra, rb)

                    D[ix3:ix3+3,jx3:jx3+3] =-Tab

                    #
                    F[ix3:ix3+3] += electric_field(cphfs[b], ra)
                
            
    # compute induction energy
    P = F.T @ numpy.linalg.inv(D)
    e =-0.5 * P.dot(F)

    return e

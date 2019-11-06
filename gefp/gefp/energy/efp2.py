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
    D = cphf.wfn().Da().to_array(dense=True) * 2.0
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



def induction(*cphfs, return_all=False):
    """
 Compute EFP2 induction energy from list of CPHF objects. 

 Returns:
   (return_all = True) induction energy, induced dipole moments vector, electric field vector
   (return_all = False) induction energy
"""

    n_frag = len(cphfs)

    # extract LMO centroids, distributed polarizabilities and other data
    lmoc = []
    dpol = []
    nmos = [x.nocc() for x in cphfs]

    for cphf in cphfs:
        lmoc_cphf = [cphf.lmo_centroid(x) for x in range(cphf.nocc())]
        lmoc.append(lmoc_cphf)
        dpol_cphf = [cphf.polarizability(x) for x in range(cphf.nocc())]
        dpol.append(dpol_cphf)

    # allocate D and F matrices
    DIM = sum(nmos)
    D = numpy.zeros((3*DIM, 3*DIM))
    F = numpy.zeros( 3*DIM)

    # compute D and F matrices
    i_pol = 0
    for n in range(n_frag):
        N   = nmos[n]  # number of LMOc's
        An  = dpol[n]  # list of distributed polarizabilities
        Ln  = lmoc[n]  # list of LMOc's

        i_pol += N
        for a in range(N):  
            ix3 = 3*(i_pol - N) + 3*a
            #print(n,a,ix3, ix3+3)

            D[ix3:ix3+3,ix3:ix3+3] = numpy.linalg.inv(An[a].to_array(dense=True))  # diagonal elements of D matrix

            ra = Ln[a].to_array(dense=True)

            j_pol = 0
            for m in range(n_frag):
              M = nmos[m]
              j_pol += M

              # consider other fragments than n-th fragment
              if (n!=m):
                Lm = lmoc[m]

                F[ix3:ix3+3] += electric_field(cphfs[m], ra)

                for b in range(M):
                    rb = Lm[b].to_array(dense=True)
                    jx3 = 3*(j_pol - M) + 3*b
                    Tab = t_2(ra, rb)

                    D[ix3:ix3+3,jx3:jx3+3] =-Tab
                    #
                
    # compute induction energy
    P = numpy.linalg.inv(D) @ F
    e =-0.5 * F @ P

    if return_all: return e, P, F
    else: return e

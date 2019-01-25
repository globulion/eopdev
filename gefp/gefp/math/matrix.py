#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Matrix module.
 Bartosz Błasiak, Gundelfingen, Jan 2019
"""

__all__ = ["Superimposer"     ,
           "rotation_matrix"  ,
           "rotate_ao_matrix" ]

import sys
import math
import numpy
import psi4

### SVDSuperimposer from BIOPYTHON PACKAGE
# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# [!] zmiana w oryginalnym kodzie:
#        w run() zmieniono na:     
#           av1=sum(coords,axis=0)/self.n  
#           av2=sum(reference_coords,axis=0)/self.n 
#        z powodu złej wersji metody sum()
class Superimposer:
    """\
 SVDSuperimposer finds the best rotation and translation to put
 two point sets on top of each other (minimizing the RMSD). This is 
 eg. useful to superimpose crystal structures.  

 SVD stands for Singular Value Decomposition, which is used to calculate
 the superposition.

 Reference:

 Matrix computations, 2nd ed. Golub, G. & Van Loan, CF., The Johns 
 Hopkins University Press, Baltimore, 1989
 """
    def __init__(self):
        self._clear()

    # Private methods

    def _clear(self):
        self.reference_coords=None
        self.coords=None
        self.transformed_coords=None
        self.rot=None
        self.tran=None
        self.rms=None
        self.init_rms=None

    def _rms(self, coords1, coords2):
        "Return rms deviations between coords1 and coords2."
        diff=coords1-coords2
        l=coords1.shape[0]
        return numpy.sqrt(numpy.sum(numpy.sum(diff*diff))/l)

    # Public methods
    
    def set(self, reference_coords, coords):
        """
        Set the coordinates to be superimposed.
        coords will be put on top of reference_coords.

        o reference_coords: an NxDIM array
        o coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        # clear everything from previous runs
        self._clear()
        # store cordinates
        self.reference_coords=reference_coords
        self.coords=coords
        n=reference_coords.shape
        m=coords.shape
        if n!=m or not(n[1]==m[1]==3):
            raise Exception("Coordinate number/dimension mismatch.")
        self.n=n[0]

    def run(self):
        "Superimpose the coordinate sets."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        coords=self.coords
        reference_coords=self.reference_coords
        # center on centroid
        av1=numpy.sum(coords,axis=0)/self.n  
        av2=numpy.sum(reference_coords,axis=0)/self.n    
        coords=coords-av1
        reference_coords=reference_coords-av2
        # correlation matrix
        a=numpy.dot(numpy.transpose(coords), reference_coords)
        u, d, vt=numpy.linalg.svd(a)
        self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        # check if we have found a reflection
        if numpy.linalg.det(self.rot)<0:
            vt[2]=-vt[2]
            self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        self.tran=av2-numpy.dot(av1, self.rot)

    def get_transformed(self):
        "Get the transformed coordinate set."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords=numpy.dot(self.coords, self.rot)+self.tran
        return self.transformed_coords

    def get_rotran(self):
        "Right multiplying rotation matrix and translation."
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        "Root mean square deviation of untransformed coordinates."
        if self.coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rms is None:
            self.init_rms=self._rms(self.coords, self.reference_coords)
        return self.init_rms

    def get_rms(self):
        "Root mean square deviation of superimposed coordinates."
        if self.rms is None:
            transformed_coords=self.get_transformed()
            self.rms=self._rms(transformed_coords, self.reference_coords)
        return self.rms

def rotation_matrix(initial=None,final=None):
    """\
 Returns rotation matrix and rms from SVD superposition of two structures (Kabsch algorithm).
 The initial structure is rotated into final one. The transformation is defined as follows:

 final = numpy.dot(initial, rot) + transl

 Returns: rot, rms
"""
    sup = Superimposer()
    sup.set(final,initial)
    sup.run()
    rms = sup.get_rms()
    rot, transl = sup.get_rotran()
    return rot, rms

def rotate_ao_matrix(M, rot_3d, bfs):
    """\
 Rotation of matrices in AO basis due to 3D rotation of basis set.

 M      - matrix in AO basis to be rotated
 rot_3d - 3 x 3 rotation matrix (Cartesian)
 bfs    - psi4.core.BasisSet in which M is represented

 Returns: rotated matrix M
"""
    # number of basis functions per shell type
    if bfs.has_puream():  nam = {0: 1, 1: 3, 2: 5, 3: 7, 4: 9}
    else:                 nam = {0: 1, 1: 3, 2: 6, 3:10, 4:15}

    max_am = bfs.max_am()
    assert max_am <= 4, "Agnular momenta larger than 4 are not supported!"
    assert bfs.has_puream() is False, "Sorry. Only Cartesian basis sets (puream = False) are supported at present."

    # indices per angular momentum
    idx = {x: [] for x in range(max_am+1)}
    for i in range(bfs.nbf()):
        i_shell = bfs.ao_to_shell(i)
        am      = bfs.shell(i_shell).am
        idx[am].append(i)
    for am in range(max_am+1):
        assert len(idx[am]) % nam[am] == 0

    # build rotation matrix
    R = numpy.identity(M.shape[0], numpy.float64)

    #if max_am > 1: r2= numpy.matrix_power(rot_3d, 2).transpose([0,2,1,3])
    #if max_am > 2: r3= numpy.matrix_power(rot_3d, 3).transpose([1-1,4-1,2-1,5-1,3-1,6-1])

    def populate_R(am, r):
        n_p_groups = int(len(idx[am]) / nam[am])
        g_c = 0
        for group in range(n_p_groups):
            g_n = g_c + nam[am]
            idx_g = idx[am][g_c:g_n]
            for ir,i in enumerate(idx_g):
                for jr,j in enumerate(idx_g):
                    R[i,j] = r[ir,jr]
            g_c+= nam[am]

    # --- s block
    None

    # --- p block
    if max_am > 0: populate_R(1, rot_3d.copy())

    # --- d block
    if max_am > 1: raise NotImplementedError
    
    # rotate
    M_rot = numpy.linalg.multi_dot([R.T, numpy.array(M), R])

    # return
    return M_rot

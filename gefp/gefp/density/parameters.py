#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Parameters Module.
 Bartosz Błasiak, Wrocław, Apr 2019
"""

from abc import ABC, abstractmethod
from .partitioning import Density
import numpy

__all__ = ["Guess", "rearrange_eigenvectors"]

class Guess(ABC):
    def __init__(self, n=None, c=None, matrix=None):
        super(Guess, self).__init__()
        self._n = None
        self._c = None
        self._matrix = None
        if n is not None: self._n = n.copy()
        if c is not None: self._c = c.copy()
        if matrix is not None: self._matrix = matrix.copy()
        self._type = None

    @classmethod
    def create(cls, n=None, c=None, matrix=None, t='matrix'):
        if matrix is None:
           if   t.lower() == 'matrix': return Matrix_Guess(n, c)
           elif t.lower() == 'nc'    : return     NC_Guess(n, c)
           else: raise ValueError("Not recognized guess type. Only MATRIX and NC are possible.")
        elif t.lower() == 'nc'       : return     NC_Guess(n, c)
        else                         : return Matrix_Guess(matrix=matrix, n=n, c=c)

    def update(self, S=None, C=None):
        if S is not None and C is not None:
         self._n, self._c = Density.natural_orbitals(self.matrix(), S, C, 
                              orthogonalize_mo=True,
                              order='descending', 
                              no_cutoff=0.0, 
                              renormalize=False, 
                              return_ao_orthogonal=False,
                              ignore_large_n=True)
        else:
         self._n, self._c = Density.natural_orbitals(self.matrix(), None, None,
                              orthogonalize_mo=False,
                              order='descending',
                              no_cutoff=0.0,
                              renormalize=False,
                              return_ao_orthogonal=False,
                              ignore_large_n=True)

    def matrix(self): 
        if self._type == 'nc': 
            return Density.generalized_density(self._n, self._c)
        else: 
            if self._matrix is None: 
                assert self._n is not None and self._c is not None
                return Density.generalized_density(self._n, self._c)
            else: 
                return self._matrix.copy()

    def copy(self): 
        return Guess.create(n=self._n, c=self._c, matrix=self.matrix(), t=self._type)

    @abstractmethod
    def pack(self): pass
    def unpack(self): return self._n, self._c

    def __add__(self, other):
        assert self._type == other._type, "Cannot add Guess objects of different type!"
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            if self._type == 'nc':
               n1 = self._n; n2 = other._n 
               c1 = self._c; c2 = other._c
               n  = n1 + n2
               c  = c1 + c2
               m  = None
            else:
               n  = None
               c  = None
               m1 = self.matrix(); m2 = other.matrix()
               if m1 is not None and m2 is not None:  m = m1 + m2
               else:                                  m = None
            new = Guess.create(n, c, m, self._type)
        else: raise TypeError("Addition of only 'Guess' type objects is defined.")
        return new

    def __sub__(self, other):
        assert self._type == other._type, "Cannot subtract Guess objects of different type!"
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            if self._type == 'nc':
               n1 = self._n; n2 = other._n 
               c1 = self._c; c2 = other._c
               n  = n1 - n2
               c  = c1 - c2
               m  = None
            else:
               n  = None
               c  = None
               m1 = self.matrix(); m2 = other.matrix()
               if m1 is not None and m2 is not None:  m = m1 - m2
               else:                                  m = None
            new = Guess.create(n, c, m, self._type)
        else: raise TypeError("Subtraction of only 'Guess' type objects is defined.")
        return new

    def __rmul__(self, other):
        if (not isinstance(other, Guess) and isinstance(self, Guess)):
            if self._type == 'nc':
               n  = self._n * other                 
               c  = self._c * other
               m  = None
            else:
               n  = None
               c  = None
               m  = self.matrix() * other
            new = Guess.create(n, c, m, self._type)
        else: raise TypeError("Right multiplication works only for 'X'*'Guess' where 'X' is instance of other class than 'Guess'.")
        return new


class NC_Guess(Guess):
    def __init__(self, n=None, c=None, matrix=None):
        super(NC_Guess, self).__init__(n, c, matrix)
        self._type = 'nc'

    def pack(self): return numpy.hstack([self._n, self._c.ravel()])

class Matrix_Guess(Guess):
    def __init__(self, n=None, c=None, matrix=None):
        super(Matrix_Guess, self).__init__(n, c, matrix)
        self._type = 'matrix'

    def pack(self): return self.matrix()


def _reorder(P,sim,axis=0):
    """Reorders the tensor according to <axis> (default is 0). 
<sim> is the list of pairs from 'order' function. 
In normal numbers (starting from 1...).
Copied from LIBBBG code."""
    P_new = numpy.zeros(P.shape,dtype=numpy.float64)
    if   axis==0:
         for i,j in sim:
             P_new[i-1] = P[j-1]
    elif axis==1:
         for i,j in sim:
             P_new[:,i-1] = P[:,j-1]
    elif axis==2:
         for i,j in sim:
             P_new[:,:,i-1] = P[:,:,j-1]
    return P_new

def _order(R,P,start=0,lprint=1):
    """order list: adapted from LIBBBG code"""
    new_P = P.copy()
    sim   = []
    rad =  []
    for i in range(len(R)-start):
        J = 0+start
        r = 1.0E+100
        rads = []
        for j in range(len(P)-start):
            r_ = numpy.sum(( R[i+start]-P[j+start])**2)
            r__= numpy.sum((-R[i+start]-P[j+start])**2)
            if r__<r_: r_=r__
            rads.append(r_)
            if r_<r:
               r=r_
               J = j
        sim.append((i+1,J+1))
        new_P[i+start] = P[J+start]
        rad.append(rads)
    for i in range(len(R)-start):
        s = numpy.sum(numpy.sign(new_P[i])/numpy.sign(R[i]))
        if lprint: print("%10d %f" %(i+1,s))
        r_ = sum(( R[i+start]-new_P[i+start])**2)
        r__= sum((-R[i+start]-new_P[i+start])**2)
       
        #if s < -154: 
        #   print "TUTAJ s < -154"
        #   #new_P[i]*=-1.
        if r__<r_:
          if lprint: print("    HERE r__ < r_ (sign reversal)")
          new_P[i]*=-1.
    return new_P, sim#, array(rad,dtype=float)

def rearrange_eigenpairs(n, u, u_ref):
    r,sim = _order(u_ref.T, u.T, lprint=0)
    w  = _reorder(u.T, sim)
    m  = _reorder(n  , sim)
    w  = w.T
    return m, w



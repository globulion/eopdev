#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Parameters Module.
 Bartosz Błasiak, Wrocław, Apr 2019
"""

from abc import ABC, abstractmethod
from .partitioning import Density

__all__ = ["Guess"]

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

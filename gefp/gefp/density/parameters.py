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
        self._n = None
        self._c = None
        self._matrix = None
        if n is not None: self._n = n.copy()
        if c is not None: self._c = c.copy()
        if matrix is not None: self._matrix = matrix.copy()
        self._type = None
        super(Guess, self).__init__()

    @classmethod
    def create(cls, n=None, c=None, matrix=None, t='matrix'):
        if matrix is None:
           if   t.lower() == 'matrix': return Matrix_Guess(n, c)
           elif t.lower() == 'nc'    : return     NC_Guess(n, c)
           else: raise ValueError("Not recognized guess type. Only MATRIX and NC are possible.")
        else: return Matrix_Guess(matrix=matrix)

    def matrix(self): 
        if self._type == 'nc': 
            return Density.generalized_density(self._n, self._c)
        else: 
            if self._matrix is None: 
                assert self._n is not None and self._c is not None
                return Density.generalized_density(self._n, self._c)
            else: 
                return self._matrix

    def copy(self): 
        if self._type == 'nc':
           return Guess.create(n=self._n.copy(), c=self._c.copy(), t=self._type)
        else:
           if self._matrix is None:
              assert self._n is not None and self._c is not None
              return Guess.create(matrix=self._matrix.copy())
           else:
              Guess.create(n=self._n.copy(), c=self._c.copy(), t=self._type)

    @abstractmethod
    def pack(self): pass
    def unpack(self): return self._n, self._c

    def __add__(self, other):
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            if self._matrix is None:
               n1 = self._n; n2 = other._n 
               c1 = self._c; c2 = other._c
               n  = n1 + n2
               c  = c1 + c2
               new = Guess.create(n, c, None, self._type)
            else:
               m1 = self._matrix; m2 = other._matrix
               m  = m1 + m2
               new = Guess.create(matrix=m)
        return new

    def __sub__(self, other):
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            if self._matrix is None:
               n1 = self._n; n2 = other._n          
               c1 = self._c; c2 = other._c
               n  = n1 - n2
               c  = c1 - c2
               new = Guess.create(n, c, None, self._type)
            else: 
               m1 = self._matrix; m2 = other._matrix
               m  = m1 - m2
               new = Guess.create(matrix=m)
        return new

    def __rmul__(self, other):
        if (not isinstance(other, Guess) and isinstance(self, Guess)):
            if self._matrix is None:
               n  = self._n * other                 
               c  = self._c * other
               new = Guess.create(n, c, None, self._type)
            else:
               new = Guess.create(matrix=self._matrix * other)
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

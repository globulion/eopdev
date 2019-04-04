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
    def __init__(self, n, c): 
        self._n = n
        self._c = c
        self._type = None
        super(Guess, self).__init__()

    @classmethod
    def create(cls, n, c, t='matrix'):
        if   t.lower() == 'matrix': return Matrix_Guess(n, c)
        elif t.lower() == 'nc'    : return     NC_Guess(n, c)
        else: raise ValueError("Not recognized guess type. Only MATRIX and NC are possible.")
    def matrix(self): return Density.generalized_density(self._n, self._c)
    def copy(self): return Guess.create(self._n.copy(), self._c.copy(), self._type)
    @abstractmethod
    def pack(self): pass
    def unpack(self): return self._n, self._c

    def __add__(self, other):
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            n1 = self._n; n2 = other._n
            c1 = self._c; c2 = other._c
            n  = n1 + n2
            c  = c1 + c2
            new = Guess.create(n, c, self._type)
        return new

    def __sub__(self, other):
        if (isinstance(self, Guess) and isinstance(other, Guess)):
            n1 = self._n; n2 = other._n
            c1 = self._c; c2 = other._c
            n  = n1 - n2
            c  = c1 - c2
            new = Guess.create(n, c, self._type)
        return new

    def __rmul__(self, other):
        if (not isinstance(other, Guess) and isinstance(self, Guess)):
            n  = self._n * other
            c  = self._c * other
            new = Guess.create(n, c, self._type)
        return new


class NC_Guess(Guess):
    def __init__(self, n, c):
        super(NC_Guess, self).__init__(n, c)
        self._type = 'nc'

    def pack(self): return numpy.hstack([self._n, self._c.ravel()])

class Matrix_Guess(Guess):
    def __init__(self, n, c):
        super(Matrix_Guess, self).__init__(n, c)
        self._type = 'matrix'

    def pack(self): return self.matrix()

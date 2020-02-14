#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Orthonormalization module.
 Bartosz Błasiak, Gundelfingen, Jan 2019
 Created as gefp module: Gundelfingen, 14 Feb 2020
"""

import numpy

__all__ = ["GrammSchmidt"]

class GrammSchmidt:
  def __init__(self, V):
      self.V = V
      self.n = V.shape[1]
  def normalize(self):
      for i in range(self.n):
          self.V[:,i]/= numpy.linalg.norm(self.V[:,i])
  def orthonormalize(self):
      self.orthogonalize()
      self.normalize()
  def orthogonalize(self):
      for k in range(1, self.n):
          for i in range(0, k):
              self.V[:,k] -= self.proj(self.V[:,i], self.V[:,k])
      #U = self.V.copy()
      #for k in range(1, self.n):
      #    U[:,k] = self.V[:,k]
      #    for i in range(0, k):
      #        U[:,k] -= self.proj(U[:,i], self.V[:,k])
      #self.V = U.copy()
  def orthogonalize_vector(self, d, normalize=False):
      for i in range(self.n):
          d -= self.proj(self.V[:,i], d)
      if normalize: d/= numpy.linalg.norm(d)
      return d
  def append(self, d):
      self.V = numpy.hstack((self.V, d[:,numpy.newaxis]))
      self.n += 1
  
  def proj(self, u, v):
      return v.dot(u)/ u.dot(u) * u 

if __name__ == '__main__':
   numpy.random.seed(0)
   V = numpy.random.random((234,6)) - 0.5

   s = Schmidt(V)
   s.orthonormalize()
   V_orig=s.V.copy()

   v = numpy.random.random((234))
   v_= s.orthogonalize_vector(v, normalize=True)
   s.append(v_)
   V = s.V
   print(V.T @ V)

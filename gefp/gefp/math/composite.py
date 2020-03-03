#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 Tensor Composite Indices Module.
 Bartosz Błasiak, Gundelfingen, Feb 2020
"""

__all__ = ["symmetry_matrix", "number_of_elements",
           "retrieve_matrix",
           "retrieve_tensor",
           "retrieve_supertensor",
           "form_futf",
           "form_superfutf",
           "partial_contraction",
           "partial_contraction_with_trace"]

import math
import numpy

def symmetry_matrix(n):
    "2D Cartesian symmetry numbers r_ij such that 1 if i==j, 2 otherwise."
    r = numpy.ones((n,n)) * 2.0
    for i in range(n): r[i,i] = 1.0
    return r

# number of composite elements {uw}
number_of_elements = lambda n: int(n*(n+1)/2)
    
def retrieve_matrix(g_triu):
    "For preliminary test and design to extend to tensor version"
    n =int ( (math.sqrt(1+g_triu.size*8)-1)/2 )
    Y = numpy.zeros((n,n))
    Y[numpy.triu_indices(n)] = g_triu.copy()
    Yt = Y.copy().transpose((1,0))
    u = numpy.diag( Yt.diagonal(axis1=0,axis2=1) )
    Yt-=u
    Y+=Yt
    return Y

def retrieve_tensor(g_triu):
    """
 -------------------------------------------------------------------------------------------------

 Retrieve full tensor from its flattened upper triangular form (FUTF) 
 correspoding to two first indices.

 -------------------------------------------------------------------------------------------------
 Description

 Consider a tensor G of dimention (n,n,[d1,d2,...]) which is symmetric
 with respect to the first two axes, i.e.,
 
 G[i,j] = G[j,i]

 The flattened upper triangular form (FUTF) is defined as

 g_triu = G[numpy.triu_indices(n)]
 
 There inverse transformation is then defined to be:
 
 G = retrieve_tensor(g_triu)

 This works for all tensors of rank r >= 2. 

 -------------------------------------------------------------------------------------------------
                                                         Bartosz Błasiak, 27 Feb 2020 Gundelfingen
"""
    n =int ( (math.sqrt(1+g_triu.shape[0]*8)-1)/2 )

    s = [n,n] + list(g_triu.shape[1:])
    s = tuple(s)
 
    Y = numpy.zeros(s)

    Y[numpy.triu_indices(n)] = g_triu.copy()

    dim = len(s)
    s = numpy.arange(dim); s[0] = 1; s[1] = 0 
    s = tuple(s)

    Yt = Y.copy().transpose(s)
    if dim == 2:
       for i in range(n): Yt[i,i] = 0.0
    else:
       for i in range(n): Yt[i,i].fill(0.0)

    Y+=Yt
    return Y

def retrieve_supertensor(H_triu, m1=1):
    """
 -------------------------------------------------------------------------------------------------

 Retrieve full super-tensor from its doubly flattened upper triangular form (FUTF).

 -------------------------------------------------------------------------------------------------
 Description

 Consider a super-tensor H of dimention (n,n,[d1,d2,...,dk] | m,m,[d1,d2,...,dl]) which 
 is symmetric with respect to the first two axes, i.e.,
 
 H[i,j] = H[j,i]

 and also symmetric with respect to m1 and m1+1 axes, i.e.,

 H[:,m1,m1+1,:] = H[:,m1+1,m1,:]

 where m1 is the first axis of dimension m.
 The flattened upper triangular form (FUTF) is defined as

 H_triu = form_superfutf(H, m1)
 
 There inverse transformation is then defined to be:
 
 H = retrieve_super_tensor(H_triu)

 This works for all tensors of rank r >= 4. 

 -------------------------------------------------------------------------------------------------
 Example

 Take an array H of the following structure:

               (n,n,[dk1,dk2,dk3] | m,m,[dl1])
                | |  |   |   |      | |  |
 axes:          0 1  2   3   4      5 6  7
                | |  |   |   |      | |  |
 For example,   | |  |   |   |      | |  |
                | |  |   |   |      | |  |
 H.shape =     (4,4, 6,  7,  8,     2,2, 4)

 and symmetric with respect to axes 0,1 as well as 5,6.
 Here, one has

  k = 3
  l = 1
  m1 = 5 

 The FUTF of H will be of shape

 H_triu.shape = (10,2,3,4,3,4)
                 -      -

 and the FUTF axes are underlined. Note that the dimension of the FUTF axis corresponding
 to an upper triangle of matrix dimension (n,n) is given by

 n_ = n*(n+1) / 2

 -------------------------------------------------------------------------------------------------
                                                         Bartosz Błasiak, 27 Feb 2020 Gundelfingen
"""
    n1 = 0
    # H_triu                          #  n_   (d1,...,dk) |  m_   (d1,...,dl)
    Y = retrieve_tensor(H_triu)       # [n,n] (d1,...,dk) |  m_   (d1,...,dl)
    
    # determine s
    o = H_triu.shape
    o1= o[:m1]
    o2= o[m1:]
    K = len(o1)
    L = len(o2)

    s1= m1 + numpy.arange(L) + 1
    s2=      numpy.arange(K+1)
    s = tuple(numpy.hstack((s1,s2)))
    Y = Y.transpose(s)                #  m_   (d1,...,dl) | [n,n] (d1,...,dk)
    Y = retrieve_tensor(Y)            # [m,m] (d1,...,dl) | [n,n] (d1,...,dk)
    #
    N1= 1 + L
    s1= N1 + numpy.arange(K+1)
    s2=      numpy.arange(L+1)

    s = tuple(numpy.hstack((s1,s2)))
    Y = Y.transpose(s)                # [n,n] (d1,...,dk) | [m,m] (d1,...,dl)
    return Y    


def form_futf(G): 
    """
 -------------------------------------------------------------------------------------------------

 Create a flattened upper triangular form (FUTF) of a tensor.

 -------------------------------------------------------------------------------------------------
 Description

 Consider a tensor G of dimention (n,n,[d1,d2,...]) which is symmetric
 with respect to the first two axes, i.e.,
 
 G[i,j] = G[j,i]

 The flattened upper triangular form (FUTF) is defined as

 g_triu = G[numpy.triu_indices(n)]
 
 There inverse transformation is then defined to be:
 
 G = retrieve_tensor(g_triu)

 This works for all tensors of rank r >= 2. 

 -------------------------------------------------------------------------------------------------
                                                         Bartosz Błasiak, 27 Feb 2020 Gundelfingen
"""
    n = G.shape[0]
    return G[numpy.triu_indices(n)]


def form_superfutf(H, m1=2):
    """
 -------------------------------------------------------------------------------------------------

 Create a doubly flattened upper triangular form (FUTF) of a super-tensor.

 -------------------------------------------------------------------------------------------------
 Description

 See description of retrieve_supertensor.

 -------------------------------------------------------------------------------------------------
                                                         Bartosz Błasiak, 27 Feb 2020 Gundelfingen
"""

    n1 = 0
    N = H.shape[n1]
    M = H.shape[m1]
    # H                               # [n,n] (d1,...,dk) | [m,m] (d1,...,dl)
    H_1 = H[numpy.triu_indices(N)]    #  n_   (d1,...,dk) | [m,m] (d1,...,dl)

    # determine s
    o = H.shape
    o1= o[:m1]
    o2= o[m1:]
    K = len(o1) - 1
    L = len(o2) - 1

    s1= m1 + numpy.arange(L + 1) - 1
    s2=      numpy.arange(K) 
    s = tuple(numpy.hstack((s1,s2)))

    H_1 = H_1.transpose(s)            # [m,m] (d1,...,dl) |  n_   (d1,...,dk)
    H_1 = H_1[numpy.triu_indices(M)]  #  m_   (d1,...,dl) |  n_   (d1,...,dk)

    N1= 1 + L
    s1= N1 + numpy.arange(K) - 1
    s2=      numpy.arange(L) 
    s = tuple(numpy.hstack((s1,s2)))
  
    H_1 = H_1.transpose(s)            #  n_   (d1,...,dk) |  m_   (d1,...,dl)
    return H_1


def partial_contraction(W, X):
    """
 -----------------------------------------------------------------------

 Partial contraction of two 2-rank tensors Y = W · X

 -----------------------------------------------------------------------

 From two matrices W and X generate a 3-rank tensor Y such that:

 Y_{a|bc} [W,X] = \sum_{p ≥ a}  W_{pb} X_{pc}

 -----------------------------------------------------------------------
                                Bartosz Błasiak, Gundelfingen 2 Mar 2020
"""
    n = W.shape[0]
    Y = numpy.zeros((n,n,n))

    v = numpy.zeros((n,n))
    for a in range(n):
        i = n-1-a
        v+= numpy.outer(W[i,:], X[i,:])
        Y[i] = v.copy()

    return Y

def partial_contraction_with_trace(W, X): 
    """
 -----------------------------------------------------------------------

 Partial contraction of two 2-rank tensors with trace Y = Tr[ W · X ]

 -----------------------------------------------------------------------

 From two matrices W and X generate a 2-rank tensor Y such that:

 Y_{a|b} [W,X] = \sum_{p ≥ a}  W_{pb} X_{pa}

 -----------------------------------------------------------------------
                                Bartosz Błasiak, Gundelfingen 2 Mar 2020
"""
    Y = W.T @ numpy.tril(X)
    return Y.T

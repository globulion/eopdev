#!/usr/bin/python

from sympy import KroneckerDelta, Symbol, cos, sin, Matrix
import sympy

x = Symbol('x')
# transformation matrix
X = lambda i,j,I,J: KroneckerDelta(i,j) * (1-KroneckerDelta(i,I))*(1-KroneckerDelta(j,J)) \
                  + cos(x) * KroneckerDelta(i,j)*(KroneckerDelta(i,I)+KroneckerDelta(j,J)) \
                  + sin(x) * (1-KroneckerDelta(i,j)) * (KroneckerDelta(i,I)*KroneckerDelta(j,J)-KroneckerDelta(i,J)*KroneckerDelta(j,I))
Xd= lambda i,j,I,J: -sin(x) * KroneckerDelta(i,j)*(KroneckerDelta(i,I)+KroneckerDelta(j,J)) \
                    +cos(x) * (1-KroneckerDelta(i,j)) * (KroneckerDelta(i,I)*KroneckerDelta(j,J)-KroneckerDelta(i,J)*KroneckerDelta(j,I))

# symbol space
sym = {0:'i', 1:'j', 2:'k', 3:'l', 4:'m', 5:'n'}
#sym = {0:'1', 1:'2', 2:'3', 3:'4', 4:'5', 5:'6'}

# R and P symbols
make_Rsymbol = lambda i,j,k,l,m,n: Symbol('R%c%c%c%c%c%c' % (sym[i],sym[j],sym[k],sym[l],sym[m],sym[n]))
make_Psymbol = lambda i,j,k      : Symbol('P%c%c%c'       % (sym[i],sym[j],sym[k]                     ))


i, j, k, l, m, n, I, J = sympy.symbols("i j k l m n I J")
dZ_R = Xd(k,i,I,J) * X (l,j,I,J) * X (m,i,I,J) * X (n,j,I,J) \
      +X (k,i,I,J) * Xd(l,j,I,J) * X (m,i,I,J) * X (n,j,I,J) \
      +X (k,i,I,J) * X (l,j,I,J) * Xd(m,i,I,J) * X (n,j,I,J) \
      +X (k,i,I,J) * X (l,j,I,J) * X (m,i,I,J) * Xd(n,j,I,J) 

dZ_P = Xd(j,i,I,J) * X(k,i,I,J) + X(j,i,I,J) * Xd(k,i,I,J)

#dZ_P_s = sympy.simplify(dZ_P)
#dZ_R_s = sympy.simplify(dZ_R)

#print dZ_P_s
print
print dZ_R


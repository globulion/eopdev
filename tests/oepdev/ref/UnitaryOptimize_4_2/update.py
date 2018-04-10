#!/usr/bin/python
"""
 Update 
"""
from numpy import *
def update_P(Pref, X):
    P = Pref.copy(); P.fill(0.0)
    N = Pref.shape[0]
    for i in range(N):
        for J in range(N): 
            for K in range(N):
                v = 0.0
                for j in range(N):
                    for k in range(N):
                        v += X[j,J] * X[k,K] * Pref[i,j,k] #* d(j,J) * d(k,K)
                P[i,J,K] = v
    return P
def update_P(Pref, X):
    P1 = tensordot(Pref, X, (1,0))
    P1 = P1.transpose(0,2,1)
    P1 = tensordot(P1, X, (2,0))
    return P1   
def update_R(Rref, X):
    N = Rref.shape[0]
    R = Rref.copy(); R.fill(0.0)
    for i in range(N):
        for j in range(N):
            for K in range(N):
                for L in range(N):
                    for M in range(N):
                        for NN in range(N):
                            v = 0.0
                            for k in range(N):
                                for l in range(N):
                                    for m in range(N):
                                        for n in range(N):
                                            v += X[k,K] * X[l,L] * X[m,M] * X[n,NN] * Rref[i,j,k,l,m,n] 
                            R[i,j,K,L,M,NN] = v
    return R
def update_R(Rref, X):
    R1 = tensordot(Rref, X, (2,0))
    R1 = R1.transpose(0,1,5,3,4,2)
    R1 = tensordot(R1  , X, (3,0))
    R1 = R1.transpose(0,1,2,5,4,3)
    R1 = tensordot(R1  , X, (4,0))
    R1 = R1.transpose(0,1,2,3,5,4)
    R1 = tensordot(R1  , X, (5,0))
    return R1

def eZ(P,X):
    Z = 0; N=len(X)
    for i in range(N):
      for j in range(N):
        for k in range(N):
          Z+=P[i,j,k] * X[j,i] * X[k,i]
    return Z
def eZI(P):
    Z = 0; N=len(P)
    for i in range(N): Z+=P[i,i,i]
    return Z
def eZR(R,X):
    Z = 0
    N = R.shape[0]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                #Z += P[i,j,k] * X[j,i] * X[k,i]
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            Z += R[i,j,k,l,m,n] * X[k,i] * X[l,j] * X[m,i] * X[n,j]
    return Z
N = 3; nX  = 7
X1 = random.random((N,N))
X2 = random.random((N,N))
X  = dot(X1,X2)
P = random.random((N,N,N))
R = random.random((N,N,N,N,N,N))
I = identity(N)

# form an ensemble of X matrices
Xlist = [ random.random((N,N)) for i in range(nX)]
# construct the final transformation matrix
Xacc = I.copy()
for x in Xlist: Xacc = dot(Xacc, x)
#Xacc = dot(dot(Xlist[0], Xlist[1]), Xlist[2])
    
# compute the cumulative Z
zr_1 = eZR(R, Xacc)
zp_1 = eZ(P,Xacc)
# compute the step-wise Z
Rt = R.copy(); Pt = P.copy()
for x in Xlist:
    Pt = update_P(Pt, x)
    Rt = update_R(Rt, x)
zr_2 = eZR(Rt, I)
zp_2 = eZ(Pt, I)
# compare and exit
print """
 Z values after %d transformations
 --------------------------------------
 Cumulative         : %10.4f %10.4f
 Step-Wise updated  : %10.4f %10.4f
""" % (nX, zr_1, zp_1, zr_2, zp_2)
exit()

# first transform
P1 = tensordot(P, X1, (1,0))
P1 = P1.transpose(0,2,1)
P1 = tensordot(P1, X1, (2,0))
# second transfo
P2 = tensordot(P1, X2, (1,0))
P2 = P2.transpose(0,2,1)
P2 = tensordot(P2, X2, (2,0))
# test of update
P1_t = update_P(P,X1)
R1  = update_R(R,X1)
R2 = update_R(R1,X2)
R_t= update_R(R, X)

zr_1 = eZR(R,X)
zr_2 = eZR(R2,I)
print zr_1, zr_2
exit()

print ((R_t-R2)**2).sum()

# full transform
Pf = tensordot(P, X, (1,0))
Pf = Pf.transpose(0,2,1)
Pf = tensordot(Pf, X, (2,0))

c = ((Pf-P2)**2).sum()
print c

print eZ(P,X)
print eZ(Pf,I)


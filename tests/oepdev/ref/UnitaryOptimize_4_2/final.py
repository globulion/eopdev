#!/usr/bin/python
"""
 Z = Z_P + Z_R
 Z_P = sum_ijk P_ijk X_ji X_ki
 Z_R = sum_ijklmn R_ijklmn X_ki X_lj X_mi X_nj
"""
import numpy

# Kronecker delta
d = lambda i,j: 0.0 if (i!=j) else 1.0
def make_X(x,n,I,J):
    "Form unitary Jacobi n-dimensional transform based on I,J Jacobi indices"
    assert I<J, "I must be smaller than J!"
    X = numpy.identity(n)
    X[I,I] = numpy.cos(x)
    X[J,J] = numpy.cos(x)
    X[I,J] = numpy.sin(x)
    X[J,I] =-numpy.sin(x)
    return X
def Z_P(P,x,I,J):
    "Compute Z_P function"
    n = len(P)
    X = make_X(x,n,I,J)
    Z = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                Z += P[i,j,k] * X[j,i] * X[k,i]
    return Z
def Z_R(R,x,I,J):
    "Compute Z_R function"
    N = len(R)
    X = make_X(x,N,I,J)
    Z = 0.0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            Z += R[i,j,k,l,m,n] * X[k,i] * X[l,j] * X[m,i] * X[n,j]
    return Z
def dZ_P_dx(P,x,I,J):
    "Compute derivative of Z_P wrt x (analytically)"
    assert I<J, "I must be smaller than J!"
    n = len(P)
    A = lambda i,j,k,III,JJJ: (d(III,j)*d(JJJ,i) - d(III,i)*d(JJJ,j)) * (1-d(i,j))
    B = lambda i,j,k,III,JJJ: -d(i,j)*(d(III,j)+d(JJJ,i))
    C = lambda i,j,k,III,JJJ:  d(i,k)*(1-d(III,k))*(1-d(JJJ,i))
    D = lambda i,j,k,III,JJJ:  d(i,k)*(d(III,k)+d(JJJ,i))
    E = lambda i,j,k,III,JJJ: -d(III,i)*d(JJJ,k)*(1-d(i,k))
    F = lambda i,j,k,III,JJJ:  d(III,k)*d(JJJ,i)*(1-d(i,k))
    G = lambda i,j,k,III,JJJ:  (1-d(i,k))*(d(III,k)*d(JJJ,i) - d(III,i)*d(JJJ,k))
    H = lambda i,j,k,III,JJJ: -d(i,k)*(d(III,k)+d(JJJ,i))
    II= lambda i,j,k,III,JJJ:  d(i,j)*(1-d(III,j))*(1-d(JJJ,i))
    JJ= lambda i,j,k,III,JJJ:  d(i,j)*(d(III,j)+d(JJJ,i))
    K = lambda i,j,k,III,JJJ: -d(III,i)*d(JJJ,j)*(1-d(i,j))
    L = lambda i,j,k,III,JJJ:  d(III,j)*d(JJJ,i)*(1-d(i,j))
    dZ = 0.0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a0 = 0.5 * (A(i,j,k,I,J)*D(i,j,k,I,J) + G(i,j,k,I,J)*JJ(i,j,k,I,J) + B(i,j,k,I,J)*(E(i,j,k,I,J)+F(i,j,k,I,J)) + H(i,j,k,I,J)*(K(i,j,k,I,J)+L(i,j,k,I,J)))                                                 
                a1 = A(i,j,k,I,J)*C(i,j,k,I,J) + G(i,j,k,I,J)*II(i,j,k,I,J)
                a2 = 0.5 * (A(i,j,k,I,J)*D(i,j,k,I,J) + G(i,j,k,I,J)*JJ(i,j,k,I,J) - B(i,j,k,I,J)*(E(i,j,k,I,J)+F(i,j,k,I,J)) - H(i,j,k,I,J)*(K(i,j,k,I,J)+L(i,j,k,I,J)))                                                 
                b1 = B(i,j,k,I,J)*C(i,j,k,I,J) + H(i,j,k,I,J)*II(i,j,k,I,J)
                b2 = 0.5 * (G(i,j,k,I,J)*(K(i,j,k,I,J)+L(i,j,k,I,J)) + H(i,j,k,I,J)*JJ(i,j,k,I,J) + B(i,j,k,I,J)*D(i,j,k,I,J) + A(i,j,k,I,J)*(E(i,j,k,I,J)+F(i,j,k,I,J)))
                t  = a0 + a1 * numpy.cos(x) + b1 * numpy.sin(x) + a2 * numpy.cos(2*x) + b2 * numpy.sin(2*x)
                dZ+= t * P[i,j,k]
    return dZ
def dZ_R_dx(R,x,I,J):
    "Compute derivative of Z_R wrt x (analytically)"
    assert I<J, "I must be smaller than J!"
    a_ = lambda i, k: (-d(I,i)*d(J,k)+d(I,k)*d(J,i))*(1.0-d(i,k))
    b_ = lambda i, k: d(i,k) * (d(I,k)+d(J,i))
    c_ = lambda i, k: d(i,k) * (1.0-d(I,k))*(1.0-d(J,i))
    N = len(R)
    dZ = 0.0
    a0 = 0
    a1 = 0 ; b1 = 0
    a2 = 0 ; b2 = 0
    a3 = 0 ; b3 = 0
    a4 = 0 ; b4 = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for m in range(N):
                        for n in range(N):
                            r = R[i,j,k,l,m,n] 
                            # First R batch
                            A = a_(i,k)
                            B =-b_(i,k)
                            C = a_(i,m)
                            D = c_(i,m)
                            E = b_(i,m)
                            F = a_(j,l)
                            G = c_(j,l)
                            H = b_(j,l)
                            II= a_(j,n)
                            JJ= c_(j,n)
                            K = b_(j,n)
                            a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                      +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                            a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                      +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                            a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                            a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                       -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                            a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                            b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                       +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                      +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                            b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                      -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                            # Second R batch
                            A = a_(i,m)
                            B =-b_(i,m)
                            C = a_(i,k)
                            D = c_(i,k)
                            E = b_(i,k)
                            F = a_(j,l)
                            G = c_(j,l)
                            H = b_(j,l)
                            II= a_(j,n)
                            JJ= c_(j,n)
                            K = b_(j,n)
                            a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                      +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                            a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                      +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                            a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                            a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                       -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                            a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                            b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                       +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                      +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                            b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                      -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                            # Third R batch
                            A = a_(j,l)
                            B =-b_(j,l)
                            C = a_(i,k)
                            D = c_(i,k)
                            E = b_(i,k)
                            F = a_(i,m)
                            G = c_(i,m)
                            H = b_(i,m)
                            II= a_(j,n)
                            JJ= c_(j,n)
                            K = b_(j,n)
                            a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                      +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                            a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                      +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                            a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                            a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                       -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                            a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                            b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                       +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                      +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                            b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                      -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0
                            # Fourth R batch
                            A = a_(j,n)
                            B =-b_(j,n)
                            C = a_(i,k)
                            D = c_(i,k)
                            E = b_(i,k)
                            F = a_(i,m)
                            G = c_(i,m)
                            H = b_(i,m)
                            II= a_(j,l)
                            JJ= c_(j,l)
                            K = b_(j,l)
                            a0 += r * (A*C*F*K+A*C*H*II+4*A*D*G*K+4*A*D*H*JJ+A*E*F*II+4*A*E*G*JJ+3*A*E*H*K\
                                      +3*B*C*F*II+4*B*C*G*JJ+B*C*H*K+4*B*D*F*JJ+4*B*D*G*II+B*E*F*K+B*E*H*II) / 8.0
                            a1 += r * (A*C*F*JJ+A*C*G*II+A*D*F*II+4*A*D*G*JJ+3*A*D*H*K+3*A*E*G*K+3*A*E*H*JJ\
                                      +B*C*G*K+B*C*H*JJ+B*D*F*K+B*D*H*II+B*E*F*JJ+B*E*G*II) / 4.0
                            a2 += r * (A*D*G*K+A*D*H*JJ+A*E*G*JJ+A*E*H*K-B*C*F*II-B*C*G*JJ-B*D*F*JJ-B*D*G*II) / 2.0
                            a3 += r * (-A*C*F*JJ-A*C*G*II-A*D*F*II+A*D*H*K+A*E*G*K+A*E*H*JJ-B*C*G*K-B*C*H*JJ\
                                       -B*D*F*K-B*D*H*II-B*E*F*JJ-B*E*G*II) / 4.0
                            a4 += r * (-A*C*F*K-A*C*H*II-A*E*F*II+A*E*H*K+B*C*F*II-B*C*H*K-B*E*F*K-B*E*H*II) / 8.0
                            b1 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II+3*B*C*F*JJ+3*B*C*G*II\
                                       +3*B*D*F*II+4*B*D*G*JJ+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b2 += r * (A*C*F*II+2*A*C*G*JJ+A*C*H*K+2*A*D*F*JJ+2*A*D*G*II+A*E*F*K+A*E*H*II+B*C*F*K\
                                      +B*C*H*II+2*B*D*G*K+2*B*D*H*JJ+B*E*F*II+2*B*E*G*JJ+B*E*H*K) / 4.0
                            b3 += r * (A*C*G*K+A*C*H*JJ+A*D*F*K+A*D*H*II+A*E*F*JJ+A*E*G*II-B*C*F*JJ-B*C*G*II\
                                      -B*D*F*II+B*D*H*K+B*E*G*K+B*E*H*JJ) / 4.0
                            b4 += r * (-A*C*F*II+A*C*H*K+A*E*F*K+A*E*H*II-B*C*F*K-B*C*H*II-B*E*F*II+B*E*H*K) / 8.0

    dZ = a0 + a1*numpy.cos(x) + a2*numpy.cos(2*x) + a3*numpy.cos(3*x) + a4*numpy.cos(4*x) 
    dZ+=      b1*numpy.sin(x) + b2*numpy.sin(2*x) + b3*numpy.sin(3*x) + b4*numpy.sin(4*x) 
    return dZ
def dZ_P_dx_numeric(P,x,I,J, step=0.01):
    "Compute derivative of Z_P wrt x (numerically)"
    Zp = Z_P(P,x+step,I,J)
    Zm = Z_P(P,x-step,I,J)
    dZ = (Zp - Zm)/(2*step)
    return dZ
def dZ_R_dx_numeric(R,x,I,J, step=0.01):
    "Compute derivative of Z_R wrt x (numerically)"
    Zp = Z_R(R,x+step,I,J)
    Zm = Z_R(R,x-step,I,J)
    dZ = (Zp - Zm)/(2*step)
    return dZ

N = 4; x = 0.453
I = 1; J = 3
P = numpy.random.random((N,N,N))
R = numpy.random.random((N,N,N,N,N,N))
zp = Z_P(P,x,I,J)
zr = Z_R(R,x,I,J)
dz_p_1 = dZ_P_dx        (P,x,I,J)
dz_p_2 = dZ_P_dx_numeric(P,x,I,J)
dz_r_1 = dZ_R_dx        (R,x,I,J)
dz_r_2 = dZ_R_dx_numeric(R,x,I,J)

log = """
 ----------------------------------------------------------------------------------
 Test Suite to compute analytical derivatives of 4-th order objective function of X:
 %s
 where:
  - X is unitary (N x N) Jacobi transformation matrix with parameter x.
  - R is general real (N x N x N x N x N x N) tensor
  - P is general real (N x N x N) tensor
 ----------------------------------------------------------------------------------
 
""" % __doc__
log += """
 Dimensionality (N) = %d
 I   = %d
 J   = %d
 x   = %14.3f
 Z_P = %14.3f
 Z_R = %14.3f
 
 Derivatives:     Analytical     Numerical
   dZ_P/dx  %14.4f %14.4f
   dZ_R/dx  %14.4f %14.4f
""" % (N, I, J, x, zp, zr, 
       dz_p_1, dz_p_2,
       dz_r_1, dz_r_2)
print log

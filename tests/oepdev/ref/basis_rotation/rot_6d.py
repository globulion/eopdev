#!/usr/bin/python3
from numpy import *

def make_r2(r):
    "Create 6 by 6 transformation matrix for 6D-type vector elements"
    import numpy
    R6 = numpy.zeros((6,6))
    s = 1.0
    # for 0 - XX
    R6[0,0] = r[0,0] * r[0,0]                            # XX XX
    R6[0,1] = r[0,0] * r[1,0] * 2.0 * s                  # XX XY
    R6[0,2] = r[0,0] * r[2,0] * 2.0 * s                  # XX XZ
    R6[0,3] = r[1,0] * r[1,0]                            # XX YY
    R6[0,4] = r[1,0] * r[2,0] * 2.0 * s                  # XX YZ
    R6[0,5] = r[2,0] * r[2,0]                            # XX ZZ
    # for 1 - XY 
    R6[1,0] = r[0,0] * r[0,1]                            # XY XX
    R6[1,1] =(r[0,0] * r[1,1] + r[1,0] * r[0,1]) * s     # XY XY
    R6[1,2] =(r[0,0] * r[2,1] + r[2,0] * r[0,1]) * s     # XY XZ
    R6[1,3] = r[1,0] * r[1,1]                            # XY YY
    R6[1,4] =(r[1,0] * r[2,1] + r[2,0] * r[1,1]) * s     # XY YZ
    R6[1,5] = r[2,0] * r[2,1]                            # XY ZZ
    # for 2 - XZ
    R6[2,0] = r[0,0] * r[0,2]                            # XZ XX
    R6[2,1] =(r[0,0] * r[1,2] + r[1,0] * r[0,2]) * s     # XZ XY
    R6[2,2] =(r[0,0] * r[2,2] + r[2,0] * r[0,2]) * s     # XZ XZ
    R6[2,3] = r[1,0] * r[1,2]                            # XZ YY
    R6[2,4] =(r[1,0] * r[2,2] + r[2,0] * r[1,2]) * s     # XZ YZ
    R6[2,5] = r[2,0] * r[2,2]                            # XZ ZZ
    # for 3 - YY
    R6[3,0] = r[0,1] * r[0,1]                            # YY XX
    R6[3,1] = r[0,1] * r[1,1] * 2.0 * s                  # YY XY
    R6[3,2] = r[0,1] * r[2,1] * 2.0 * s                  # YY XZ
    R6[3,3] = r[1,1] * r[1,1]                            # YY YY
    R6[3,4] = r[1,1] * r[2,1] * 2.0 * s                  # YY YZ
    R6[3,5] = r[2,1] * r[2,1]                            # YY ZZ
    # for 4 - YZ
    R6[4,0] = r[0,1] * r[0,2]                            # YZ XX
    R6[4,1] =(r[0,1] * r[1,2] + r[1,1] * r[0,2]) * s     # YZ XY
    R6[4,2] =(r[0,1] * r[2,2] + r[2,1] * r[0,2]) * s     # YZ XZ
    R6[4,3] = r[1,1] * r[1,2]                            # YZ YY
    R6[4,4] =(r[1,1] * r[2,2] + r[2,1] * r[1,2]) * s     # YZ YZ
    R6[4,5] = r[2,1] * r[2,2]                            # YZ ZZ
    # for 5 - ZZ
    R6[5,0] = r[0,2] * r[0,2]                            # ZZ XX
    R6[5,1] = r[0,2] * r[1,2] * 2.0 * s                  # ZZ XY
    R6[5,2] = r[0,2] * r[2,2] * 2.0 * s                  # ZZ XZ
    R6[5,3] = r[1,2] * r[1,2]                            # ZZ YY
    R6[5,4] = r[1,2] * r[2,2] * 2.0 * s                  # ZZ YZ
    R6[5,5] = r[2,2] * r[2,2]                            # ZZ ZZ
    return R6.T

def delta(a,b): return 1.0 if a==b else 0.0
idx = {(0,0):0, (0,1):1, (0,2):2, (1,1):3, (1,2):4, (2,2):5}
def make_r2_test(r):
    R = identity(6)
    for ap in range(3):
     for bp in range(3):
      if ap<=bp:
         apbp = idx[(ap,bp)]
         for a in range(3):
          for b in range(3):
           if a<=b:
              ab = idx[(a,b)]
              R[apbp,ab] = r[ap,a]*r[bp,b] + r[bp,a]*r[ap,b]*(1.0-delta(ap,bp))
               
    return R

def make_a(A):
    a = zeros(6)
    a[0] = A[0,0]; a[1] = A[0,1]; a[2] = A[0,2]; a[3] = A[1,1]; a[4] = A[1,2]; a[5] = A[2,2]
    return a

def make_rotation(x,y,z):
    rot_x = array([1.0,0.0,0.0,0.0,cos(x),-sin(x),0.0,sin(x),cos(x)]).reshape(3,3)
    rot_y = array([cos(y),0.0,sin(y),0.0,1.0,0.0,-sin(y),0.0,cos(y)]).reshape(3,3)
    rot_z = array([cos(z),-sin(z),0.0,sin(z),cos(z),0.0,0.0,0.0,1.0]).reshape(3,3)
    return rot_z @ rot_y @ rot_x

random.seed(0)
A = random.random((3,3)); A+=A.T
a = make_a(A)
r = make_rotation(0.4,-4.0,1.3)

R1 = make_r2(r)
R2 = make_r2_test(r)

print(R1-R2)

# test
A_rot = make_a(einsum("ij,ia,jb->ab",A,r,r))
a_rot = dot(R1.T, a)
print(A_rot)
print()
print(a_rot)


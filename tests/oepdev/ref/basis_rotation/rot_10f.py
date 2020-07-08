#!/usr/bin/python3
from numpy import *

def delta(a,b): return 1.0 if a==b else 0.0
def dif(a,b): return 1.0 - delta(a,b)
idx = {}
I = 0
for i in range(3):
    for j in range(3):
        if i<=j:
           for k in range(3):
               if j<=k:
                  idx[(i,j,k)] = I
                  I+=1

class S:
  language = None
  def __init__(self,a,b,c,r):
      self.__a = a
      self.__b = b
      self.__c = c
      self.__r = r
  def __call__(self,ap,bp,cp):
    return self.__r[ap,self.__a]*self.__r[bp,self.__b]*r[cp,self.__c]
  def get_term(self, ap, bp, cp):
      if S.language == 'py': 
         log = "r[%d,%d]*r[%d,%d]*r[%d,%d]" % (ap,self.__a,bp,self.__b,cp,self.__c)
      else:
         log = "r%d%d*r%d%d*r%d%d" % (ap,self.__a,bp,self.__b,cp,self.__c)
      return log

def make_r3_code(language='py'):
    log = ""
    data = {}
    S.language = language

    for ap in range(3):           
     for bp in range(3):
      if ap<=bp:
       for cp in range(3):
        if bp<=cp:
         apbpcp = idx[(ap,bp,cp)]
         for a in range(3):
          for b in range(3):
           if a<=b:
            for c in range(3):
             if b<=c:
              abc = idx[(a,b,c)]
              pair = (apbpcp,abc)
              data[pair] = []
              s = S(a,b,c,None)

              # A group
              if (bp==cp):
                  l = s.get_term(ap,bp,bp)
                  #log += "R[%d,%d] += %s\n" % (apbpcp, abc, l) 
                  data[pair].append(l)
                  if ap!=bp:
                     data[pair].append( s.get_term(bp,ap,bp) )
                     data[pair].append( s.get_term(bp,bp,ap) )

              # B group
              if (ap==bp and bp!=cp):
                     data[pair].append( s.get_term(cp,ap,ap) )
                     data[pair].append( s.get_term(ap,cp,ap) )
                     data[pair].append( s.get_term(ap,ap,cp) )

              # C group
              if (ap!=bp and bp!=cp):
                     data[pair].append( s.get_term(ap,bp,cp) )
                     data[pair].append( s.get_term(ap,cp,bp) )
                     data[pair].append( s.get_term(bp,ap,cp) )
                     data[pair].append( s.get_term(bp,cp,ap) )
                     data[pair].append( s.get_term(cp,ap,bp) )
                     data[pair].append( s.get_term(cp,bp,ap) )

                    
                   
  
    for pair,terms in data.items():
        if language == 'py': line = " R[%d,%d] = " % pair
        else               : line = " R10[%d][%d] = " % pair
        line+= " + ".join(terms)
        if language == 'c++': line+= ";"

        log+= line + "\n"
    print(log)

def make_r3_generated_from_code_generator(r):
    R = identity(10)
    R[0,0] = r[0,0]*r[0,0]*r[0,0]
    R[0,1] = r[0,0]*r[0,0]*r[0,1]
    R[0,2] = r[0,0]*r[0,0]*r[0,2]
    R[0,3] = r[0,0]*r[0,1]*r[0,1]
    R[0,4] = r[0,0]*r[0,1]*r[0,2]
    R[0,5] = r[0,0]*r[0,2]*r[0,2]
    R[0,6] = r[0,1]*r[0,1]*r[0,1]
    R[0,7] = r[0,1]*r[0,1]*r[0,2]
    R[0,8] = r[0,1]*r[0,2]*r[0,2]
    R[0,9] = r[0,2]*r[0,2]*r[0,2]
    R[1,0] = r[1,0]*r[0,0]*r[0,0] + r[0,0]*r[1,0]*r[0,0] + r[0,0]*r[0,0]*r[1,0]
    R[1,1] = r[1,0]*r[0,0]*r[0,1] + r[0,0]*r[1,0]*r[0,1] + r[0,0]*r[0,0]*r[1,1]
    R[1,2] = r[1,0]*r[0,0]*r[0,2] + r[0,0]*r[1,0]*r[0,2] + r[0,0]*r[0,0]*r[1,2]
    R[1,3] = r[1,0]*r[0,1]*r[0,1] + r[0,0]*r[1,1]*r[0,1] + r[0,0]*r[0,1]*r[1,1]
    R[1,4] = r[1,0]*r[0,1]*r[0,2] + r[0,0]*r[1,1]*r[0,2] + r[0,0]*r[0,1]*r[1,2]
    R[1,5] = r[1,0]*r[0,2]*r[0,2] + r[0,0]*r[1,2]*r[0,2] + r[0,0]*r[0,2]*r[1,2]
    R[1,6] = r[1,1]*r[0,1]*r[0,1] + r[0,1]*r[1,1]*r[0,1] + r[0,1]*r[0,1]*r[1,1]
    R[1,7] = r[1,1]*r[0,1]*r[0,2] + r[0,1]*r[1,1]*r[0,2] + r[0,1]*r[0,1]*r[1,2]
    R[1,8] = r[1,1]*r[0,2]*r[0,2] + r[0,1]*r[1,2]*r[0,2] + r[0,1]*r[0,2]*r[1,2]
    R[1,9] = r[1,2]*r[0,2]*r[0,2] + r[0,2]*r[1,2]*r[0,2] + r[0,2]*r[0,2]*r[1,2]
    R[2,0] = r[2,0]*r[0,0]*r[0,0] + r[0,0]*r[2,0]*r[0,0] + r[0,0]*r[0,0]*r[2,0]
    R[2,1] = r[2,0]*r[0,0]*r[0,1] + r[0,0]*r[2,0]*r[0,1] + r[0,0]*r[0,0]*r[2,1]
    R[2,2] = r[2,0]*r[0,0]*r[0,2] + r[0,0]*r[2,0]*r[0,2] + r[0,0]*r[0,0]*r[2,2]
    R[2,3] = r[2,0]*r[0,1]*r[0,1] + r[0,0]*r[2,1]*r[0,1] + r[0,0]*r[0,1]*r[2,1]
    R[2,4] = r[2,0]*r[0,1]*r[0,2] + r[0,0]*r[2,1]*r[0,2] + r[0,0]*r[0,1]*r[2,2]
    R[2,5] = r[2,0]*r[0,2]*r[0,2] + r[0,0]*r[2,2]*r[0,2] + r[0,0]*r[0,2]*r[2,2]
    R[2,6] = r[2,1]*r[0,1]*r[0,1] + r[0,1]*r[2,1]*r[0,1] + r[0,1]*r[0,1]*r[2,1]
    R[2,7] = r[2,1]*r[0,1]*r[0,2] + r[0,1]*r[2,1]*r[0,2] + r[0,1]*r[0,1]*r[2,2]
    R[2,8] = r[2,1]*r[0,2]*r[0,2] + r[0,1]*r[2,2]*r[0,2] + r[0,1]*r[0,2]*r[2,2]
    R[2,9] = r[2,2]*r[0,2]*r[0,2] + r[0,2]*r[2,2]*r[0,2] + r[0,2]*r[0,2]*r[2,2]
    R[3,0] = r[0,0]*r[1,0]*r[1,0] + r[1,0]*r[0,0]*r[1,0] + r[1,0]*r[1,0]*r[0,0]
    R[3,1] = r[0,0]*r[1,0]*r[1,1] + r[1,0]*r[0,0]*r[1,1] + r[1,0]*r[1,0]*r[0,1]
    R[3,2] = r[0,0]*r[1,0]*r[1,2] + r[1,0]*r[0,0]*r[1,2] + r[1,0]*r[1,0]*r[0,2]
    R[3,3] = r[0,0]*r[1,1]*r[1,1] + r[1,0]*r[0,1]*r[1,1] + r[1,0]*r[1,1]*r[0,1]
    R[3,4] = r[0,0]*r[1,1]*r[1,2] + r[1,0]*r[0,1]*r[1,2] + r[1,0]*r[1,1]*r[0,2]
    R[3,5] = r[0,0]*r[1,2]*r[1,2] + r[1,0]*r[0,2]*r[1,2] + r[1,0]*r[1,2]*r[0,2]
    R[3,6] = r[0,1]*r[1,1]*r[1,1] + r[1,1]*r[0,1]*r[1,1] + r[1,1]*r[1,1]*r[0,1]
    R[3,7] = r[0,1]*r[1,1]*r[1,2] + r[1,1]*r[0,1]*r[1,2] + r[1,1]*r[1,1]*r[0,2]
    R[3,8] = r[0,1]*r[1,2]*r[1,2] + r[1,1]*r[0,2]*r[1,2] + r[1,1]*r[1,2]*r[0,2]
    R[3,9] = r[0,2]*r[1,2]*r[1,2] + r[1,2]*r[0,2]*r[1,2] + r[1,2]*r[1,2]*r[0,2]
    R[4,0] = r[0,0]*r[1,0]*r[2,0] + r[0,0]*r[2,0]*r[1,0] + r[1,0]*r[0,0]*r[2,0] + r[1,0]*r[2,0]*r[0,0] + r[2,0]*r[0,0]*r[1,0] + r[2,0]*r[1,0]*r[0,0]
    R[4,1] = r[0,0]*r[1,0]*r[2,1] + r[0,0]*r[2,0]*r[1,1] + r[1,0]*r[0,0]*r[2,1] + r[1,0]*r[2,0]*r[0,1] + r[2,0]*r[0,0]*r[1,1] + r[2,0]*r[1,0]*r[0,1]
    R[4,2] = r[0,0]*r[1,0]*r[2,2] + r[0,0]*r[2,0]*r[1,2] + r[1,0]*r[0,0]*r[2,2] + r[1,0]*r[2,0]*r[0,2] + r[2,0]*r[0,0]*r[1,2] + r[2,0]*r[1,0]*r[0,2]
    R[4,3] = r[0,0]*r[1,1]*r[2,1] + r[0,0]*r[2,1]*r[1,1] + r[1,0]*r[0,1]*r[2,1] + r[1,0]*r[2,1]*r[0,1] + r[2,0]*r[0,1]*r[1,1] + r[2,0]*r[1,1]*r[0,1]
    R[4,4] = r[0,0]*r[1,1]*r[2,2] + r[0,0]*r[2,1]*r[1,2] + r[1,0]*r[0,1]*r[2,2] + r[1,0]*r[2,1]*r[0,2] + r[2,0]*r[0,1]*r[1,2] + r[2,0]*r[1,1]*r[0,2]
    R[4,5] = r[0,0]*r[1,2]*r[2,2] + r[0,0]*r[2,2]*r[1,2] + r[1,0]*r[0,2]*r[2,2] + r[1,0]*r[2,2]*r[0,2] + r[2,0]*r[0,2]*r[1,2] + r[2,0]*r[1,2]*r[0,2]
    R[4,6] = r[0,1]*r[1,1]*r[2,1] + r[0,1]*r[2,1]*r[1,1] + r[1,1]*r[0,1]*r[2,1] + r[1,1]*r[2,1]*r[0,1] + r[2,1]*r[0,1]*r[1,1] + r[2,1]*r[1,1]*r[0,1]
    R[4,7] = r[0,1]*r[1,1]*r[2,2] + r[0,1]*r[2,1]*r[1,2] + r[1,1]*r[0,1]*r[2,2] + r[1,1]*r[2,1]*r[0,2] + r[2,1]*r[0,1]*r[1,2] + r[2,1]*r[1,1]*r[0,2]
    R[4,8] = r[0,1]*r[1,2]*r[2,2] + r[0,1]*r[2,2]*r[1,2] + r[1,1]*r[0,2]*r[2,2] + r[1,1]*r[2,2]*r[0,2] + r[2,1]*r[0,2]*r[1,2] + r[2,1]*r[1,2]*r[0,2]
    R[4,9] = r[0,2]*r[1,2]*r[2,2] + r[0,2]*r[2,2]*r[1,2] + r[1,2]*r[0,2]*r[2,2] + r[1,2]*r[2,2]*r[0,2] + r[2,2]*r[0,2]*r[1,2] + r[2,2]*r[1,2]*r[0,2]
    R[5,0] = r[0,0]*r[2,0]*r[2,0] + r[2,0]*r[0,0]*r[2,0] + r[2,0]*r[2,0]*r[0,0]
    R[5,1] = r[0,0]*r[2,0]*r[2,1] + r[2,0]*r[0,0]*r[2,1] + r[2,0]*r[2,0]*r[0,1]
    R[5,2] = r[0,0]*r[2,0]*r[2,2] + r[2,0]*r[0,0]*r[2,2] + r[2,0]*r[2,0]*r[0,2]
    R[5,3] = r[0,0]*r[2,1]*r[2,1] + r[2,0]*r[0,1]*r[2,1] + r[2,0]*r[2,1]*r[0,1]
    R[5,4] = r[0,0]*r[2,1]*r[2,2] + r[2,0]*r[0,1]*r[2,2] + r[2,0]*r[2,1]*r[0,2]
    R[5,5] = r[0,0]*r[2,2]*r[2,2] + r[2,0]*r[0,2]*r[2,2] + r[2,0]*r[2,2]*r[0,2]
    R[5,6] = r[0,1]*r[2,1]*r[2,1] + r[2,1]*r[0,1]*r[2,1] + r[2,1]*r[2,1]*r[0,1]
    R[5,7] = r[0,1]*r[2,1]*r[2,2] + r[2,1]*r[0,1]*r[2,2] + r[2,1]*r[2,1]*r[0,2]
    R[5,8] = r[0,1]*r[2,2]*r[2,2] + r[2,1]*r[0,2]*r[2,2] + r[2,1]*r[2,2]*r[0,2]
    R[5,9] = r[0,2]*r[2,2]*r[2,2] + r[2,2]*r[0,2]*r[2,2] + r[2,2]*r[2,2]*r[0,2]
    R[6,0] = r[1,0]*r[1,0]*r[1,0]
    R[6,1] = r[1,0]*r[1,0]*r[1,1]
    R[6,2] = r[1,0]*r[1,0]*r[1,2]
    R[6,3] = r[1,0]*r[1,1]*r[1,1]
    R[6,4] = r[1,0]*r[1,1]*r[1,2]
    R[6,5] = r[1,0]*r[1,2]*r[1,2]
    R[6,6] = r[1,1]*r[1,1]*r[1,1]
    R[6,7] = r[1,1]*r[1,1]*r[1,2]
    R[6,8] = r[1,1]*r[1,2]*r[1,2]
    R[6,9] = r[1,2]*r[1,2]*r[1,2]
    R[7,0] = r[2,0]*r[1,0]*r[1,0] + r[1,0]*r[2,0]*r[1,0] + r[1,0]*r[1,0]*r[2,0]
    R[7,1] = r[2,0]*r[1,0]*r[1,1] + r[1,0]*r[2,0]*r[1,1] + r[1,0]*r[1,0]*r[2,1]
    R[7,2] = r[2,0]*r[1,0]*r[1,2] + r[1,0]*r[2,0]*r[1,2] + r[1,0]*r[1,0]*r[2,2]
    R[7,3] = r[2,0]*r[1,1]*r[1,1] + r[1,0]*r[2,1]*r[1,1] + r[1,0]*r[1,1]*r[2,1]
    R[7,4] = r[2,0]*r[1,1]*r[1,2] + r[1,0]*r[2,1]*r[1,2] + r[1,0]*r[1,1]*r[2,2]
    R[7,5] = r[2,0]*r[1,2]*r[1,2] + r[1,0]*r[2,2]*r[1,2] + r[1,0]*r[1,2]*r[2,2]
    R[7,6] = r[2,1]*r[1,1]*r[1,1] + r[1,1]*r[2,1]*r[1,1] + r[1,1]*r[1,1]*r[2,1]
    R[7,7] = r[2,1]*r[1,1]*r[1,2] + r[1,1]*r[2,1]*r[1,2] + r[1,1]*r[1,1]*r[2,2]
    R[7,8] = r[2,1]*r[1,2]*r[1,2] + r[1,1]*r[2,2]*r[1,2] + r[1,1]*r[1,2]*r[2,2]
    R[7,9] = r[2,2]*r[1,2]*r[1,2] + r[1,2]*r[2,2]*r[1,2] + r[1,2]*r[1,2]*r[2,2]
    R[8,0] = r[1,0]*r[2,0]*r[2,0] + r[2,0]*r[1,0]*r[2,0] + r[2,0]*r[2,0]*r[1,0]
    R[8,1] = r[1,0]*r[2,0]*r[2,1] + r[2,0]*r[1,0]*r[2,1] + r[2,0]*r[2,0]*r[1,1]
    R[8,2] = r[1,0]*r[2,0]*r[2,2] + r[2,0]*r[1,0]*r[2,2] + r[2,0]*r[2,0]*r[1,2]
    R[8,3] = r[1,0]*r[2,1]*r[2,1] + r[2,0]*r[1,1]*r[2,1] + r[2,0]*r[2,1]*r[1,1]
    R[8,4] = r[1,0]*r[2,1]*r[2,2] + r[2,0]*r[1,1]*r[2,2] + r[2,0]*r[2,1]*r[1,2]
    R[8,5] = r[1,0]*r[2,2]*r[2,2] + r[2,0]*r[1,2]*r[2,2] + r[2,0]*r[2,2]*r[1,2]
    R[8,6] = r[1,1]*r[2,1]*r[2,1] + r[2,1]*r[1,1]*r[2,1] + r[2,1]*r[2,1]*r[1,1]
    R[8,7] = r[1,1]*r[2,1]*r[2,2] + r[2,1]*r[1,1]*r[2,2] + r[2,1]*r[2,1]*r[1,2]
    R[8,8] = r[1,1]*r[2,2]*r[2,2] + r[2,1]*r[1,2]*r[2,2] + r[2,1]*r[2,2]*r[1,2]
    R[8,9] = r[1,2]*r[2,2]*r[2,2] + r[2,2]*r[1,2]*r[2,2] + r[2,2]*r[2,2]*r[1,2]
    R[9,0] = r[2,0]*r[2,0]*r[2,0]
    R[9,1] = r[2,0]*r[2,0]*r[2,1]
    R[9,2] = r[2,0]*r[2,0]*r[2,2]
    R[9,3] = r[2,0]*r[2,1]*r[2,1]
    R[9,4] = r[2,0]*r[2,1]*r[2,2]
    R[9,5] = r[2,0]*r[2,2]*r[2,2]
    R[9,6] = r[2,1]*r[2,1]*r[2,1]
    R[9,7] = r[2,1]*r[2,1]*r[2,2]
    R[9,8] = r[2,1]*r[2,2]*r[2,2]
    R[9,9] = r[2,2]*r[2,2]*r[2,2]
    return R


def make_r3_test(r):
    R = identity(10)
    for ap in range(3):
     for bp in range(3):
      if ap<=bp:
       for cp in range(3):
        if bp<=cp:
         apbpcp = idx[(ap,bp,cp)]
         for a in range(3):
          for b in range(3):
           if a<=b:
            for c in range(3):
             if b<=c:
              abc = idx[(a,b,c)]
              s = S(a,b,c,r)
              #R[apbpcp,abc] = r[ap,a]*r[bp,b]*r[cp,c] \
              #               +r[bp,a]*r[ap,b]*r[cp,c] * dif(a,b) \
              #               +r[ap,a]*r[cp,b]*r[bp,c] * dif(b,c) \
              #               +r[cp,a]*r[bp,b]*r[ap,c] * dif(a,c) \
              #               +r[bp,a]*r[cp,b]*r[ap,c] * dif(a,b) * dif(a,c) \
              #               +r[cp,a]*r[ap,b]*r[bp,c] * dif(a,b) * dif(a,c)
              A = s(ap,bp,bp)+dif(ap,bp)*(s(bp,ap,bp)+s(bp,bp,ap)) 
             #B = dif(bp,cp)*(s(cp,ap,ap)+dif(bp,cp)*(s(ap,cp,ap)+s(ap,ap,cp)))
              B = dif(bp,cp)*(s(cp,ap,ap)+s(ap,cp,ap)+s(ap,ap,cp)) # simpler than above but equivalent
              C = s(ap,bp,cp)+s(ap,cp,bp)+s(bp,ap,cp)+s(bp,cp,ap)+s(cp,ap,bp)+s(cp,bp,ap)
              R[apbpcp,abc] = delta(bp,cp) * A + delta(ap,bp) * B + dif(ap,bp) * dif(bp,cp) * C
    return R

def make_a(A):
    a = zeros(10)
    for key in idx.keys():
        a[idx[key]] = A[key]
    return a

def make_rotation(x,y,z):
    rot_x = array([1.0,0.0,0.0,0.0,cos(x),-sin(x),0.0,sin(x),cos(x)]).reshape(3,3)
    rot_y = array([cos(y),0.0,sin(y),0.0,1.0,0.0,-sin(y),0.0,cos(y)]).reshape(3,3)
    rot_z = array([cos(z),-sin(z),0.0,sin(z),cos(z),0.0,0.0,0.0,1.0]).reshape(3,3)
    return rot_z @ rot_y @ rot_x

def make_random_A():
    A = zeros((3,3,3))
    for i in range(3):
      for j in range(3):
       if j>=i:
        for k in range(3):
         if k>=j:
          rijk = random.random()
          A[i,j,k] = rijk
          A[i,k,j] = rijk
          A[j,i,k] = rijk
          A[j,k,i] = rijk
          A[k,i,j] = rijk
          A[k,j,i] = rijk
    return A

random.seed(0)
A = make_random_A()
a = make_a(A)
r = make_rotation(0.4,-4.0,1.3)

R = make_r3_test(r)
R_= make_r3_generated_from_code_generator(r)


# test
A_rot = make_a(einsum("ijk,ia,jb,kc->abc",A,r,r,r))
a_rot = dot(R_.T, a)
print(A_rot)
print()
print(a_rot)

make_r3_code(language='py')
make_r3_code(language='c++')

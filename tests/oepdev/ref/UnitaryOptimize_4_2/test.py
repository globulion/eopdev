#!/usr/bin/python
from numpy import *
random.seed(1)

from libbbg.utilities import UnitaryOptimizer_4_2 as UO

N = 3
R = random.random((N,N,N,N,N,N)) - 1.0
P = random.random((N,N,N))-1.0

s = UO(R,P,conv=1e-8,maxiter=100,verbose=True)
s.minimize()
Xmin = s.X
s.maximize()
Xmax = s.X

def p(X, symbol):
    "C++ readible declaration of const double list"
    log = "const double %s[%d] = {" % (symbol, len(X.ravel()))
    for i in range(len(X.ravel())-1):
        log += "%13.4E," % X.ravel()[i]
        if not i%7: log+= "\n"
    log += "%13.4E };\n" % X.ravel()[-1]
    print log

p(R   , 'R')
p(P   , 'P')
p(Xmin, 'Xmin_ref')
p(Xmax, 'Xmax_ref')

print Xmin
print Xmax


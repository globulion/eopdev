#!/usr/bin/python3
"""
 Extended Density Fitting Helper Library.

 Useful routines for AO basis set optimization can be found here. 

 BB, 29.07.2020, Gundelfingen
"""

import psi4
import numpy
import scipy.optimize
import scipy.linalg
import oepdev
from . import parameters

__all__ = ["compute_v", "projection", "projected_t", "projected_o", 
           "obj_numpy", "obj_oepdev", "find_aux_ao_mini", "optimize_ao_mini"]

matrix_power = scipy.linalg.fractional_matrix_power

def compute_v(Ca, Da, prim, left_axis):
    "Compute Fock-Like OEP Matrix V(TEST, O)"
    bfs_test = left_axis
    bfs_prim = prim
    # one-electron part
    mints = psi4.core.MintsHelper(bfs_prim)
    V_nuc= mints.ao_potential(bfs_test, bfs_prim).to_array(dense=True)
    V  = numpy.dot(V_nuc, Ca)

    # two-electron part
    f_pppt = psi4.core.IntegralFactory(bfs_prim, bfs_prim, bfs_prim, bfs_test)
    nt = bfs_test.nbf()
    Ca_ = psi4.core.Matrix.from_array(Ca.copy(), "")
    Da_ = psi4.core.Matrix.from_array(Da.copy(), "")
    V_2el = oepdev.calculate_OEP_basisopt_V(nt, f_pppt, Ca_, Da_)
    V += V_2el.to_array(dense=True)
    psi4.core.clean()
    return V

def projection(c_a, s_ab, s_bb):
    "R. Polly et.al, Mol.Phys 102 (2004), 2311-2321, Eq.(35-36)"
    s_bb1 = numpy.linalg.inv(s_bb)
    A = c_a.T @ s_ab
    T = A @ s_bb1
    T = T @ A.T
    return T

def projected_t(c_a, s_ab, s_bb):
    "Projected c_b orbitals"
    T = projection(c_a, s_ab, s_bb)
    t = matrix_power(T, -0.5).real
    s_bb1 = numpy.linalg.inv(s_bb)
    A = c_a.T @ s_ab
    c_b = s_bb1 @ A.T
    c_b = c_b @ t
    return c_b

def projected_o(c_a, s_ab, s_bb):
    "Overlap after projection"
    t = projected_t(c_a, s_ab, s_bb)
    o = matrix_power(t,  0.5).real
    return o

def obj_numpy(param, t_i, bsf_i, dfbasis):
    "Objective function in pure NumPy (slower)"
    bsf_m = dfbasis.basisset(param)
    mints = psi4.core.MintsHelper(bsf_i)
    s_im = mints.ao_overlap(bsf_i, bsf_m).to_array(dense=True)
    s_mm = mints.ao_overlap(bsf_m, bsf_m).to_array(dense=True)
    t = projection(t_i, s_im, s_mm)
   #o_mi = matrix_power(t,  1.0).real
    o_mi = matrix_power(t,  0.5).real
    return -o_mi.trace()

def obj_oepdev(param, t_i, bsf_i, dfbasis, mints):
    "Objective function in pure C++ level (approx 3 times faster)"
    bsf_m = dfbasis.basisset(param)
    z = oepdev.bs_optimize_projection(t_i, mints, bsf_m, bsf_i)
    return z

def find_aux_mo_mini(G, S, I=None, eps=0.0001):
    "Compute TB matrix"
    S05 = matrix_power(S, 0.5).real
    Sm05= matrix_power(S,-0.5).real

    G_ = S05 @ G
   #G_ = G

    C = G_@ G_.T
    l, u = numpy.linalg.eigh(C)                                            # u: (aux, new)
    i = numpy.argsort(l)[::-1]
    l = l[i]
    u = u[:,i]
    print(" Largest Eigenvalues of GG Covariance in orthogonal |i> basis")
    for i in range(len(l)):
        if l[i]>= 0.00000001: print(" %4d %14.8f" % (i+1, l[i]))

    if I is None:
       I = 0;
       for i in range(l.size):
           if l[i]< eps: break
           I+=1

    print(" Optimal OEP basis size = %d out of %d all basis set"% (I, G.shape[0]) )
 
    T_= u[:,:I]                                                        # (aux, new_crop)

    T = Sm05 @ T_
    return T

#class TakeMyStep(object):
#   def __init__(self, stepsize=0.5):
#       self.stepsize = stepsize
#   def __call__(self, x):
#       x[ 0] += self._b(3.00    )   # H  1s 1 e 
#       x[ 1] += self._b(1.00    )   # H  1s 1 c
#       x[ 2] += self._b(2.00    )   # H  1s 2 e 
#       x[ 3] += self._b(1.00    )   # H  1s 2 c
#       x[ 4] += self._b(1.00    )   # H  1s 3 e 
#       x[ 5] += self._b(1.00    )   # H  1s 3 c
#       #
#       x[ 6] += self._b(300.0   )   # O  1s 1 e 
#       x[ 7] += self._b(1.00    )   # O  1s 1 c
#       x[ 8] += self._b(60.0    )   # O  1s 2 e 
#       x[ 9] += self._b(1.00    )   # O  1s 2 c
#       x[10] += self._b(10.0    )   # O  1s 3 e 
#       x[11] += self._b(1.00    )   # O  1s 3 c
#       #
#       x[12] += self._b(30.0    )   # O  2s 1 e 
#       x[13] += self._b(1.00    )   # O  2s 1 c
#       x[14] += self._b(5.0     )   # O  2s 2 e 
#       x[15] += self._b(1.00    )   # O  2s 2 c
#       x[16] += self._b(1.00    )   # O  2s 3 e 
#       x[17] += self._b(1.00    )   # O  2s 3 c
#       #
#       x[18] += self._b(30.0    )   # O  2p 1 e 
#       x[19] += self._b(1.00    )   # O  2p 1 c
#       x[20] += self._b(5.0     )   # O  2p 2 e 
#       x[21] += self._b(1.00    )   # O  2p 2 c
#       x[22] += self._b(1.00    )   # O  2p 3 e 
#       x[23] += self._b(1.00    )   # O  2p 3 c
#       return x
#   def _b(self, a):
#       return numpy.random.uniform(-a*self.stepsize, a*self.stepsize)
#
c1 = {'type':'eq', 'fun': lambda x: x[ 1]+x[ 3]+x[ 5]-1.0}
c2 = {'type':'eq', 'fun': lambda x: x[ 7]+x[ 9]+x[11]-1.0}
c3 = {'type':'eq', 'fun': lambda x: x[13]+x[15]+x[17]-1.0}
c4 = {'type':'eq', 'fun': lambda x: x[19]+x[21]+x[23]-1.0}
c = [c1, c2, c3, c4]

def optimize_ao_mini(t_i, bsf_i, dfbasis, opt_global, cpp=False):
    "Target routine for AO basis set optimization"
    param_0 = dfbasis.param
    print(" Initial Z = %14.6f" % obj_numpy(param_0, t_i, bsf_i, dfbasis))

    options = {"disp": True, "maxiter": 2000, "ftol": 1.e-9, "iprint": 8}
    if cpp:
       mints = psi4.core.MintsHelper(bsf_i)
       T_i = psi4.core.Matrix.from_array(t_i)
       ARGS = (T_i, bsf_i, dfbasis, mints)
       OBJ = obj_oepdev
    else:
       ARGS = (t_i, bsf_i, dfbasis)
       OBJ = obj_numpy

    if not opt_global:
       res = scipy.optimize.minimize(OBJ, param_0, args=ARGS, tol=1.e-9, method='slsqp', 
                   options=options, bounds=dfbasis.bounds)
    else:
       take_step = parameters.TakeMyStandardSteps(dfbasis.scales)
       res = scipy.optimize.basinhopping(OBJ, param_0, niter=10, 
                              T=0.02, stepsize=0.5,
                              minimizer_kwargs={"method": 'slsqp', "options": options,
                                  "bounds": dfbasis.bounds, "args": ARGS, "constraints": dfbasis.constraints},
                              callback=None, interval=30, disp=True, niter_success=None,
                              take_step=take_step ).lowest_optimization_result #, accept_test=accept_test)
       print(res.message)
       print(" Optimal Z = %14.6f" % res.fun)  
    if not res.success: 
       raise ValueError("Optimization is not succesfull!")
    else: 
       dfbasis.param = res.x
       bsf = dfbasis.basisset()
    return dfbasis

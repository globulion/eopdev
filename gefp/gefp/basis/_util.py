#!/usr/bin/python3
"""
 Local utilities (protected interface)
 BB, 29.07.2020, Gundelfingen
"""

def COMPARE(a, b, i=None):
    if i is None:
       A = a.copy().ravel()
       B = B.copy().ravel()
    else:
       A = a[:,i]
       B = b[:,i]
    for i in range(len(A)):
       print("%14.6f %14.6f" % (A[i], B[i]))
   #print()

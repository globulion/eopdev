#!/usr/bin/python
"""
 Multipole components. In Psi4.
"""
components = lambda lx, ly, lz: "x"*lx + "y"*ly + "z"*lz
def print_order(order):
    I = 0
    for l in range(order+1):
        for ii in range(l+1):
            lx = l - ii
            for lz in range(ii+1):
                ly = ii - lz
                print "%3d %10s" % (I, components(lx, ly, lz))
                I += 1

print_order(4)

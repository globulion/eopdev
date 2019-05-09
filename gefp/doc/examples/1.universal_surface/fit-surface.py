#!/usr/bin/python3
"""
 Fit the universal surface of DMFT functional and plot it.
 Usage: ./fit-surface.py
 Note : modify input data within the script 
"""
from surface import fit_surface, plot_surface
def func_t_non(Z): return Z
def func_t_inv(Z): return 1./Z

# Input Data: Surface File and Function of Z
surface_file = 'test.dat'
func_t       = func_t_non

# Input Data: 2D Pade Approximant Form
a_ord = [(0, 0), (0, 1), (1, 0), (1, 1)]
b_ord = [(0, 1), (1, 0), (1, 1)]

# Fit and Plot
data = fit_surface(surface_file, a_ord, b_ord, func_t=func_t)
plot_surface(data)
print(" Fitted parameters:")
print(data.par)

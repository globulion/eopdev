#!/usr/bin/python3
"""
 Standard Templates for AO Basis Set Optimization

 BB, 04.08.2020, Gundelfingen
"""

atoms_by_row = {}
atoms_by_row[1] = ("H", "He")
atoms_by_row[2] = ("Li", "Be", "B", "C", "N", "O", "F", "Ne")
atoms_by_row[3] = ("Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar")
atoms_by_row[4] = ("K", "Ca")


standard_template_by_row = {}
standard_template_by_row[1] = '''\
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

standard_template_by_row[2] = '''\
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

standard_template_by_row[3] = '''\
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

standard_template_by_row[4] = '''\
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
S   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
P   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

standard_bounds_codes_by_row = {}
b1 = " E C E C E C"
standard_bounds_codes_by_row[1] =   b1 .split()
standard_bounds_codes_by_row[2] =(3*b1).split()
standard_bounds_codes_by_row[3] =(5*b1).split()
standard_bounds_codes_by_row[4] =(7*b1).split()

standard_guess_parameters_by_atom = {}
reference_symbol_by_row = {}
s = 0.3

# EFP2-OEP (Repulsion)
reference_symbol_by_row["efp2-rep"] = {}
reference_symbol_by_row["efp2-rep"][1] = "H"
reference_symbol_by_row["efp2-rep"][2] = "O"
reference_symbol_by_row["efp2-rep"][3] = "S"
reference_symbol_by_row["efp2-rep"][4] = "K"
standard_scales_by_row = {}
standard_scales_by_row["efp2-rep"] = {}
standard_scales_by_row["efp2-rep"][1] = \
  [3.00, s, 
   2.00, s, 
   1.00, s]
standard_scales_by_row["efp2-rep"][2] = \
  [ 80.00, s, 
    10.00, s, 
     2.00, s,
   
    24.00, s,
     5.00, s,
     1.00, s,

    20.00, s,
     5.00, s,
     1.00, s,]
standard_scales_by_row["efp2-rep"][3] = []
standard_scales_by_row["efp2-rep"][4] = []
standard_guess_parameters_by_atom["efp2-rep"] = {}
standard_guess_parameters_by_atom["efp2-rep"]["H"] = [\
12.8   ,          0.43, 
 1.23  ,          0.48, 
 0.24  ,          1.00-0.43-0.48,]
standard_guess_parameters_by_atom["efp2-rep"]["O"] = [\
1266.6 ,          0.42, 
123.55 ,          0.36,
22.41  ,          1.00-0.42-0.36,
10.21  ,          0.19,
3.13   ,          0.50,
0.83   ,          1.00-0.19-0.50,
62.8   ,          0.24,
9.91   ,          0.45,
1.79   ,          1.00-0.24-0.45,]
standard_guess_parameters_by_atom["efp2-rep"]["C"] = standard_guess_parameters_by_atom["efp2-rep"]["O"]
standard_guess_parameters_by_atom["efp2-rep"]["N"] = standard_guess_parameters_by_atom["efp2-rep"]["O"]
standard_guess_parameters_by_atom["efp2-rep"]["S"] = [] #TODO

# EFP2-OEP (Charge Transfer)
reference_symbol_by_row["efp2-ct"] = {}
reference_symbol_by_row["efp2-ct"][1] = "H"
reference_symbol_by_row["efp2-ct"][2] = "O"
reference_symbol_by_row["efp2-ct"][3] = "S"
reference_symbol_by_row["efp2-ct"][4] = "K"
standard_scales_by_row["efp2-ct"] = {}
standard_scales_by_row["efp2-ct"][1] = standard_scales_by_row["efp2-rep"][1]
standard_scales_by_row["efp2-ct"][2] = standard_scales_by_row["efp2-rep"][2]
standard_scales_by_row["efp2-ct"][3] = standard_scales_by_row["efp2-rep"][3]
standard_scales_by_row["efp2-ct"][4] = standard_scales_by_row["efp2-rep"][4]
standard_guess_parameters_by_atom["efp2-ct"] = {}
standard_guess_parameters_by_atom["efp2-ct"]["H"] = standard_guess_parameters_by_atom["efp2-rep"]["H"]
standard_guess_parameters_by_atom["efp2-ct"]["O"] = standard_guess_parameters_by_atom["efp2-rep"]["O"]
standard_guess_parameters_by_atom["efp2-ct"]["C"] = standard_guess_parameters_by_atom["efp2-rep"]["O"]
standard_guess_parameters_by_atom["efp2-ct"]["N"] = standard_guess_parameters_by_atom["efp2-rep"]["O"]
standard_guess_parameters_by_atom["efp2-ct"]["S"] = standard_guess_parameters_by_atom["efp2-rep"]["S"]


# EET TODO
reference_symbol_by_row["eet"] = {}
reference_symbol_by_row["eet"][1] = "H"
reference_symbol_by_row["eet"][2] = "O"
reference_symbol_by_row["eet"][3] = "S"
reference_symbol_by_row["eet"][4] = "K"
standard_guess_parameters_by_atom["eet"] = {}

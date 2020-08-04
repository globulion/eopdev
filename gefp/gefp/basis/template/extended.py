#!/usr/bin/python3
"""
 Extended Templates for AO Basis Set Optimization

 BB, 04.08.2020, Gundelfingen
"""

extended_template_by_row = {}
extended_template_by_row[1] = '''\
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

extended_template_by_row[2] = '''\
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
D   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

extended_template_by_row[3] = '''\
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
D   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

extended_template_by_row[4] = '''\
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
D   3   1.00
%20.10f             %20.10f
%20.10f             %20.10f
%20.10f             %20.10f
****
'''

extended_bounds_codes_by_row = {}
b1 = " E C E C E C"
extended_bounds_codes_by_row[1] =(2*b1).split()
extended_bounds_codes_by_row[2] =(4*b1).split()
extended_bounds_codes_by_row[3] =(6*b1).split()
extended_bounds_codes_by_row[4] =(8*b1).split()

extended_guess_parameters_by_atom = {}
s = 0.3

# EFP2-OEP (Repulsion)
extended_scales_by_row = {}
extended_scales_by_row["efp2-rep"] = {}
extended_scales_by_row["efp2-rep"][1] = \
  [3.00, s, 
   2.00, s, 
   1.00, s,
   3.00, s,
   2.00, s,
   1.00, s]
extended_scales_by_row["efp2-rep"][2] = \
  [ 80.00, s, 
    10.00, s, 
     2.00, s,
   
    24.00, s,
     5.00, s,
     1.00, s,

    20.00, s,
     5.00, s,
     1.00, s,

    10.00, s,
     5.00, s,
     1.00, s,]
extended_scales_by_row["efp2-rep"][3] = []
extended_scales_by_row["efp2-rep"][4] = []
extended_guess_parameters_by_atom["efp2-rep"] = {}
extended_guess_parameters_by_atom["efp2-rep"]["H"] = [\
12.8   ,          0.43, 
 1.23  ,          0.48, 
 0.24  ,          1.00-0.43-0.48,
12.8   ,          0.43, 
 1.23  ,          0.48, 
 0.24  ,          1.00-0.43-0.48,]
extended_guess_parameters_by_atom["efp2-rep"]["O"] = [\
1266.6 ,          0.42, 
123.55 ,          0.36,
22.41  ,          1.00-0.42-0.36,
10.21  ,          0.19,
3.13   ,          0.50,
0.83   ,          1.00-0.19-0.50,
62.8   ,          0.24,
9.91   ,          0.45,
1.79   ,          1.00-0.24-0.45,
42.8   ,          0.24,
7.91   ,          0.45,
0.79   ,          1.00-0.24-0.45,]
extended_guess_parameters_by_atom["efp2-rep"]["C"] = extended_guess_parameters_by_atom["efp2-rep"]["O"]
extended_guess_parameters_by_atom["efp2-rep"]["N"] = extended_guess_parameters_by_atom["efp2-rep"]["O"]
extended_guess_parameters_by_atom["efp2-rep"]["S"] = []

# EFP2-OEP (Charge Transfer)
extended_scales_by_row["efp2-ct"] = {}
extended_scales_by_row["efp2-ct"][1] = extended_scales_by_row["efp2-rep"][1]
extended_scales_by_row["efp2-ct"][2] = extended_scales_by_row["efp2-rep"][2]
extended_scales_by_row["efp2-ct"][3] = extended_scales_by_row["efp2-rep"][3]
extended_scales_by_row["efp2-ct"][4] = extended_scales_by_row["efp2-rep"][4]
extended_guess_parameters_by_atom["efp2-ct"] = {}
extended_guess_parameters_by_atom["efp2-ct"]["H"] = extended_guess_parameters_by_atom["efp2-rep"]["H"]
extended_guess_parameters_by_atom["efp2-ct"]["O"] = extended_guess_parameters_by_atom["efp2-rep"]["O"]
extended_guess_parameters_by_atom["efp2-ct"]["C"] = extended_guess_parameters_by_atom["efp2-rep"]["O"]
extended_guess_parameters_by_atom["efp2-ct"]["N"] = extended_guess_parameters_by_atom["efp2-rep"]["O"]
extended_guess_parameters_by_atom["efp2-ct"]["S"] = []


# EET TODO
extended_guess_parameters_by_atom["eet"] = {}

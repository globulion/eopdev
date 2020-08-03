#!/usr/bin/python3
"""
 Parameters for Basis Set Optimization Module

 Contains:
   o standard templates for auxiliary minimal AO basis set
   o guess parameters for auxiliary minimal AO basis set optimization
   o automatized basin hopping bound and step adjustment tools

 BB, 30.07.2020, Gundelfingen
"""
import numpy

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
standard_guess_parameters_by_atom["efp2-rep"]["S"] = []

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
standard_guess_parameters_by_atom["efp2-ct"]["S"] = []

# test
standard_guess_parameters_by_atom["efp2-ct"]["H"] = [\
12.8   ,          0.43, 
1.323  ,          0.48, 
 6.24  ,          1.00-0.43-0.48,]
standard_guess_parameters_by_atom["efp2-ct"]["N"] = [\
 800.6 ,          0.42, 
 67.55 ,          0.36,
  2.1  ,          1.00-0.42-0.36,
143.1  ,          0.19,
66.3   ,          0.50,
9.83   ,          1.00-0.19-0.50,
40.8   ,          0.24,
9.91   ,          0.45,
1.79   ,          1.00-0.24-0.45,]


# EET
reference_symbol_by_row["eet"] = {}
reference_symbol_by_row["eet"][1] = "H"
reference_symbol_by_row["eet"][2] = "O"
reference_symbol_by_row["eet"][3] = "S"
reference_symbol_by_row["eet"][4] = "K"
standard_guess_parameters_by_atom["eet"] = {}

class StandardizedInput:
  def __init__(self, mol, oep_type):
      self.mol = mol
      self.oep_type = oep_type.lower()
      self.template = None
      self.parameters = None
      self.constraints = None
      self.bounds_codes = None
      self.scales = None
      self.prepare_standard_template_and_starting_parameters()

  def prepare_standard_template_and_starting_parameters(self):
      "Create a standard template and starting parameters based on the molecule input"
      line = lambda symbol: '%s   0\n' % symbol
  
      atoms = self.get_atom_symbols()
  
      template = "cartesian\n****\n"
      parameters = []
      bounds_codes = []
      scales = []

      for symbol in atoms:
          I = self.get_row(symbol)

          template += line(symbol)
          template += standard_template_by_row[I]
  
          try:
             parameters += standard_guess_parameters_by_atom[self.oep_type][symbol]
          except KeyError:
             reference_symbol = reference_symbol_by_row[self.oep_type][I]
             parameters += standard_guess_parameters_by_atom[self.oep_type][reference_symbol]

          bounds_codes += standard_bounds_codes_by_row[I]
          scales += standard_scales_by_row[self.oep_type][I]
  
      self.parameters = numpy.array(parameters)
      self.template = template
      self.bounds_codes = numpy.array(bounds_codes, dtype=str)
      self.scales = numpy.array(scales)
      self.constraints = self.get_constraints(scales)

  def get_constraints(self, scales):
      constraints = []
      #idx = list(range(1,len(scales),2))
      #n_constraints = int(len(idx)/3)
      def const(x):
          i = int(x.size/6)
         #t = (x[1::2].reshape(i,3).sum(axis=1) - 1.0).prod()
          t = (x[1::2].reshape(i,3).sum(axis=1) - 1.0)
          t = t * t
          return t.sum()
          
      #for n in range(n_constraints):
      #    i = idx[3*n]
      #    print(i)
      #    constraint = {'type':'eq', 'fun': lambda x: x[i+0] + x[i+2] + x[i+4] - 1.0}
      #    constraints.append(constraint)
      #exit()
      constraints.append({'type':'eq', 'fun': const})
      return constraints

  def get_row(self, symbol):
      I = None
      for row,atom_list in atoms_by_row.items():
          if symbol in atom_list: 
             I = row
             break
      if I is None:
         raise ValueError("Atom %s has not been found in the StandardInput library" % symbol)
      return I

  def get_atom_symbols(self):
      atoms = []
      for i in range(self.mol.natom()):
          s = self.mol.symbol(i)
          if s not in atoms: atoms.append(s)
      return atoms


    
class TakeMyStandardSteps(object):
   def __init__(self, scales, stepsize=1.0):
       self.stepsize = stepsize
       self.__scales = scales
   def __call__(self, x):
       s = self.stepsize
       for i in range(len(self.__scales)):
           x[i] += self._b(self.__scales[i])
       return x
   def _b(self, a):
       return numpy.random.uniform(-a*self.stepsize, a*self.stepsize)

#class MyBounds(object):
#    def __init__(self):
#        xmin = [0.02]*7
#        xmax = [2500.0]*7
#        self.xmax = numpy.array(xmax)
#        self.xmin = numpy.array(xmin)
#    def __call__(self, **kwargs):
#        x = kwargs["x_new"]
#        tmax = bool(numpy.all(x <= self.xmax))
#        tmin = bool(numpy.all(x >= self.xmin))
#        return tmax and tmin

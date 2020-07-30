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

# EFP2-OEP (Repulsion)
reference_symbol_by_row["efp2-rep"] = {}
reference_symbol_by_row["efp2-rep"][1] = "H"
reference_symbol_by_row["efp2-rep"][2] = "O"
reference_symbol_by_row["efp2-rep"][3] = "S"
reference_symbol_by_row["efp2-rep"][4] = "K"
standard_guess_parameters_by_atom["efp2-rep"] = {}
standard_guess_parameters_by_atom["efp2-rep"]["H"] = [\
12.8   ,          0.43, 
 1.23  ,          0.48, 
 0.24  ,          0.09,]
standard_guess_parameters_by_atom["efp2-rep"]["O"] = [\
1266.6 ,          0.42         , 
123.55 ,          0.36         ,
22.41  ,          0.22         ,
10.21  ,          0.19         ,
3.13   ,          0.50         ,
0.83   ,          0.31         ,
62.8   ,          0.24         ,
9.91   ,          0.45         ,
1.79   ,          0.31         ,]
standard_guess_parameters_by_atom["efp2-rep"]["C"] = []
standard_guess_parameters_by_atom["efp2-rep"]["N"] = []
standard_guess_parameters_by_atom["efp2-rep"]["S"] = []

# EFP2-OEP (Charge Transfer)
reference_symbol_by_row["efp2-ct"] = {}
reference_symbol_by_row["efp2-ct"][1] = "H"
reference_symbol_by_row["efp2-ct"][2] = "O"
reference_symbol_by_row["efp2-ct"][3] = "S"
reference_symbol_by_row["efp2-ct"][4] = "K"
standard_guess_parameters_by_atom["efp2-ct"] = {}
standard_guess_parameters_by_atom["efp2-ct"]["H"] = []
standard_guess_parameters_by_atom["efp2-ct"]["O"] = []
standard_guess_parameters_by_atom["efp2-ct"]["C"] = []
standard_guess_parameters_by_atom["efp2-ct"]["N"] = []
standard_guess_parameters_by_atom["efp2-ct"]["S"] = []

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
      self.prepare_standard_template_and_starting_parameters()

  def prepare_standard_template_and_starting_parameters(self):
      "Create a standard template and starting parameters based on the molecule input"
      line = lambda symbol: '%s   0\n' % symbol
  
      atoms = self.get_atom_symbols()
  
      template = "cartesian\n****\n"
      parameters = []
      bounds_codes = []

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
  
      self.parameters = numpy.array(parameters)
      self.template = template
      self.bounds_codes = numpy.array(bounds_codes, dtype=str)
      self.constraints = ()

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


    
class TakeMyStep(object):
   def __init__(self, stepsize=1.0):
       self.stepsize = stepsize
   def __call__(self, x):
       s = self.stepsize
       #x[0] += numpy.random.uniform(-2.*s, 2.*s)
       #x[1:] += numpy.random.uniform(-s, s, x[1:].shape)
       x[ 0] += self._b(1.00    )   # H  1s
       x[ 1] += self._b(4.00    )   # H  1p
       x[ 2] += self._b(200.0   )   # N  1s
       x[ 3] += self._b(30.0    )   # N  2s
       x[ 4] += self._b(10.0    )   # N  2p.1
       x[ 5] += self._b(10.0    )   # N  2p.2
       x[ 6] += self._b(10.0    )   # N  3d
       return x
   def _b(self, a):
       return numpy.random.uniform(-a*self.stepsize, a*self.stepsize)

class MyBounds(object):
    def __init__(self):
        xmin = [0.02]*7
        xmax = [2500.0]*7
        self.xmax = numpy.array(xmax)
        self.xmin = numpy.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(numpy.all(x <= self.xmax))
        tmin = bool(numpy.all(x >= self.xmin))
        return tmax and tmin


mystep = TakeMyStep()
mytest = MyBounds()



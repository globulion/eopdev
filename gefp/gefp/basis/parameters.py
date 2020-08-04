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
from .template.standard import *
from .template.extended import *

class StandardizedInput:
  def __init__(self, mol, oep_type, standard='standard'):
      self.mol = mol
      self.oep_type = oep_type.lower()
      self.standard = standard.lower()
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

      if self.standard == 'standard':
         template_by_row           = standard_template_by_row
         guess_parameters_by_atom  = standard_guess_parameters_by_atom
         bounds_codes_by_row       = standard_bounds_codes_by_row
         scales_by_row             = standard_scales_by_row
      else:
         template_by_row           = extended_template_by_row
         guess_parameters_by_atom  = extended_guess_parameters_by_atom
         bounds_codes_by_row       = extended_bounds_codes_by_row
         scales_by_row             = extended_scales_by_row

      for symbol in atoms:
          I = self.get_row(symbol)

          template += line(symbol)
          template += template_by_row[I]
  
          try:
             parameters += guess_parameters_by_atom[self.oep_type][symbol]
          except KeyError:
             reference_symbol = reference_symbol_by_row[self.oep_type][I]
             parameters += guess_parameters_by_atom[self.oep_type][reference_symbol]

          bounds_codes += bounds_codes_by_row[I]
          scales += scales_by_row[self.oep_type][I]
  
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

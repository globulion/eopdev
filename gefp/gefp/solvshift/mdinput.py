import os, re, numpy
from .fragment import Fragment

def text_to_list(text,delimiter=None,dtype=int,dt=1):
    """
Transform one text line into a list. Returns list of integers.

Syntax: 
       text_to_list(text,delimiter=None)

Usage:

1) text = '     2 3 4 5 6    4    '

text_to_list(text) = [2 3 4 5 6 4]

2) text = ' 2-4'

text_to_list(text) = [2 3 4]

3) text = ' 3, 56, 6, -1, 33  '

text_to_list(text,delimiter=',') = [3 56 6 -1 33] 

4) text = ' -0.44, 0.122, 23.4544, -335.2   '
text_to_list(text,delimiter=',',dtype=float) = [-0.44, 0.122, 23.4544, -335.2]

"""
    if (('-' in text and '.' not in text) or ('-' in text and ',' not in text)):
       start,end = text.split('-')
       rang = numpy.arange(dtype(start),dtype(end)+1,dt)
    else:
       rang = text.split(delimiter)
       rang = list(map(dtype,rang))
    return numpy.array(rang)

class MDInput:
   """
 -------------------------------------------------------------------
                   Represents SLV MD input file
 -------------------------------------------------------------------
                                               April 23, 2014
                                               
 Import:
 
   from solvshift.md import MDInput
   
 Usage:
 
   inp  = MDInput(file)
   args, idx = inp.get()
 
 The <args> list and slicing indices <idx> can be now set to 
 solvshift.slvpar.Frag.set method:

   Frag_instance.set(frame[idx],*args)

 <args> is a list and contains the following information:

   args = [ ind, nmol, bsm, supl, reord ]

 For the details about these memorials read Frag class documentation
 by typing 

   from solvshift.slvpar import Frag
   help(Frag)

 Note: argument of MDInput class can be also a string containing 
       the contents of input file (instead of providing path to
       a file)

 -------------------------------------------------------------------
                              Last Revision:  9 May 2017 @Globulion
"""
   def __init__(self,file_name):
       self.__file_name = file_name if isinstance(file_name,str) else '<<from file>>'
       self.__bsm  = list()
       self.__nmol = list()
       self.__ind  = list()
       self.__supl = list()
       self.__reord= list()
       self.__frag_idx = 0
       self.__frame_slice = list()

       self.__call__(file_name)

   def get(self):
       """
Return the options for Frag.set method and slicing list

Usage:

    args, idx = MDInputInstance.get()
"""
       args = [ self.__ind, self.__nmol, self.__bsm, self.__supl, self.__reord ]
       idx  = self.__frame_slice
       return args, idx

   # auxiliary methods
   def get_idx(self):
       """get slicing list"""
       return self.__frame_slice
   
   def get_args(self):
       """get arguments for Frag.set method"""
       args = [ self.__ind, self.__nmol, self.__bsm, self.__supl, self.__reord ]
       return args

   def get_efp(self,frame):
       """parse EFP coordinates from MD frame"""
       return frame[self.__frame_slice]

   def __call__(self,file_name):
       if os.path.exists(file_name):
          fileo = open(file_name)
          text  = fileo.read()
          fileo.close()
       elif '\n' in file_name:
          text = file_name
       else:
          raise IOError(" The input file does not exists.")
       self.__input_text = text
       templ = re.compile(r'\$Frag',re.DOTALL)
       sections = re.split(templ,text)[1:]
       for section in sections:
           self._read(section)
         
       return

   def _read(self,section):
       """read one unique (Sol)EFP fragment"""
       b_reorder = False
       b_supl    = False
       if 'reorder' in section: b_reorder = True
       if 'supl'    in section: b_supl    = True
       lines = section.split('\n')[1:]
       if lines[0].startswith('!'): return

       # read the fragment parameters
       frag = Fragment(lines[0].split()[1])

       # reorder the fragment
       reord = None
       for line in lines:
           if line.startswith('!'):continue
           if 'reorder' in line:
               reord = numpy.array( map(int,line.split()[1:]), int ) - 1
       #if b_reorder:   # deprecated! Now reordering applies to MD structure
       #   frag.reorder(reord)

       # add superimpositon list
       supl = None
       for line in lines:
           if line.startswith('!'):continue
           if 'supl' in line:
               ll = line.split('supl')[-1]
               if ',' in ll: supl = text_to_list(ll, delimiter=',') - 1
               else:         supl = text_to_list(ll) - 1

       # parse atoms for (Sol)EFP fragment
       for line in lines:
           if line.startswith('!'):continue
           if 'atoms' in line:
               atoms, n_frags = line.split()[1:]
               if ',' in atoms: atoms   = text_to_list(atoms, delimiter=',')
               else:            atoms   = text_to_list(atoms)
               n_frags = int(n_frags)
            
               if len(atoms)==1: n_atoms = 1 
               else: n_atoms = atoms[-1]-atoms[0]+1
               merror  = 'MDInputError: Invalid atomic indices or fragment numbers in fragment %i' % (self.__frag_idx + 1)
               merror += '\n -----> %s' % line
               assert n_atoms%n_frags==0, merror
               n_atoms_per_mol = int( (atoms[-1]-atoms[0]+1)/n_frags )
               for i in range(n_frags):
                   self.__ind.append(self.__frag_idx)
                   self.__nmol.append(n_atoms_per_mol)
               atoms = list(numpy.array(atoms)-1)
               self.__frame_slice += atoms

       # --- in the case of no reordering, set the appropriate list
       if reord is None:
          reord = numpy.array( [ i for i in range(n_atoms_per_mol) ], int)

       # --- append to bsm and supl memorial lists
       self.__bsm.append(frag)
       self.__supl.append(supl)
       self.__reord.append(reord)
      
       # go to the next fragment
       self.__frag_idx += 1 
       return

   def __repr__(self):
       """print the input file"""
       N = 50
       log = '-'*N+'\n'
       log+= 'SLV-MD input file: '
       log+= self.__file_name+'\n\n'+self.__input_text
       log+= '-'*N+'\n'
       return str(log)

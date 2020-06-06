#!*-* coding: utf-8 *-*
# ---------------------------------------------- #
#     DISTRIBUTED MULTIPOLE ANALYSIS MODULE      #
# ---------------------------------------------- #

__all__ = ['DMA']

import numpy, numpy.linalg, copy, units, utilities2

class SVDSuperimposer(object):
    """
    >>>>> Copied from LIBBBG.utilities <<<<<
     SVDSuperimposer adapted from BIOPYTHON PACKAGE

     Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
     This code is part of the Biopython distribution and governed by its
     license.  Please see the LICENSE file that should have been included
     as part of this package.

     ------ LICENCE ------ 

     Biopython is currently released under the "Biopython License Agreement"
     (given in full below). Unless stated otherwise in individual file headers,
     all Biopython's files are under the "Biopython License Agreement".
     
     Some files are explicitly dual licensed under your choice of the
     "Biopython License Agreement" or the "BSD 3-Clause License" (both given
     in full below). This is with the intention of later offering all of
     Biopython under this dual licensing approach.
     
     
     ------------------------------------------------------------------------------
     
     
     Biopython License Agreement
     ---------------------------
     
     Permission to use, copy, modify, and distribute this software and its
     documentation with or without modifications and for any purpose and
     without fee is hereby granted, provided that any copyright notices
     appear in all copies and that both those copyright notices and this
     permission notice appear in supporting documentation, and that the
     names of the contributors or copyright holders not be used in
     advertising or publicity pertaining to distribution of the software
     without specific prior permission.
     
     THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
     WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
     WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
     CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
     OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
     OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
     OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
     OR PERFORMANCE OF THIS SOFTWARE.
     
     
     ------------------------------------------------------------------------------
     
     
     BSD 3-Clause License
     --------------------
     
     Copyright (c) 1999-2019, The Biopython Contributors
     All rights reserved.
     
     Redistribution and use in source and binary forms, with or without
     modification, are permitted provided that the following conditions are met:
     
     1. Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
     
     2. Redistributions in binary form must reproduce the above copyright notice,
        this list of conditions and the following disclaimer in the documentation
        and/or other materials provided with the distribution.
     
     3. Neither the name of the copyright holder nor the names of its contributors
        may be used to endorse or promote products derived from this software
        without specific prior written permission.
     
     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     ------
    
     [!] zmiana w oryginalnym kodzie:
            w run() zmieniono na:     
               av1=sum(coords,axis=0)/self.n  
               av2=sum(reference_coords,axis=0)/self.n 
            z powodu złej wersji metody sum()

    ----

    SVDSuperimposer finds the best rotation and translation to put
    two point sets on top of each other (minimizing the RMSD). This is 
    eg. useful to superimpose crystal structures.  

    SVD stands for Singular Value Decomposition, which is used to calculate
    the superposition.

    Reference:

    Matrix computations, 2nd ed. Golub, G. & Van Loan, CF., The Johns 
    Hopkins University Press, Baltimore, 1989
    """
    def __init__(self):
        self._clear()

    # Private methods

    def _clear(self):
        self.reference_coords=None
        self.coords=None
        self.transformed_coords=None
        self.rot=None
        self.tran=None
        self.rms=None
        self.init_rms=None

    def _rms(self, coords1, coords2):
        "Return rms deviations between coords1 and coords2."
        diff=coords1-coords2
        l=coords1.shape[0]
        return numpy.sqrt(numpy.sum(numpy.sum(diff*diff))/l)

    # Public methods
    
    def set(self, reference_coords, coords):
        """
        Set the coordinates to be superimposed.
        coords will be put on top of reference_coords.

        o reference_coords: an NxDIM array
        o coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        # clear everything from previous runs
        self._clear()
        # store cordinates
        self.reference_coords=reference_coords
        self.coords=coords
        n=reference_coords.shape
        m=coords.shape
        if n!=m or not(n[1]==m[1]==3):
            raise Exception("Coordinate number/dimension mismatch.")
        self.n=n[0]

    def run(self):
        "Superimpose the coordinate sets."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        coords=self.coords
        reference_coords=self.reference_coords
        # center on centroid
        av1=numpy.sum(coords,axis=0)/self.n  
        av2=numpy.sum(reference_coords,axis=0)/self.n    
        coords=coords-av1
        reference_coords=reference_coords-av2
        # correlation matrix
        a=numpy.dot(numpy.transpose(coords), reference_coords)
        u, d, vt=numpy.linalg.svd(a)
        self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        # check if we have found a reflection
        if numpy.linalg.det(self.rot)<0:
            vt[2]=-vt[2]
            self.rot=numpy.transpose(numpy.dot(numpy.transpose(vt), numpy.transpose(u)))
        self.tran=av2-numpy.dot(av1, self.rot)

    def get_transformed(self):
        "Get the transformed coordinate set."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords=numpy.dot(self.coords, self.rot)+self.tran
        return self.transformed_coords

    def get_rotran(self):
        "Right multiplying rotation matrix and translation."
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        "Root mean square deviation of untransformed coordinates."
        if self.coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rms is None:
            self.init_rms=self._rms(self.coords, self.reference_coords)
        return self.init_rms

    def get_rms(self):
        "Root mean square deviation of superimposed coordinates."
        if self.rms is None:
            transformed_coords=self.get_transformed()
            self.rms = self._rms(transformed_coords, self.reference_coords)
        return self.rms


def ParseDMA(file,type='coulomb',hexadecapoles=True):
    """\
>>>>> Copied from LIBBBG.utilities <<<<<
============================================================================
Parse DMA from GAUSSIAN, GAMESS or COULOMB.py file. It returns a DMA object.
Usage:
ParseDMA(file, type='coulomb', hexadecapoles=False)
----------------------------------------------------------------------------
<type>s:
1) coulomb  or c (or just nothing - it is default)
2) gaussian or gau
3) gamess   or gms
Notes:
Gaussian file has to have pop=ChelpG printout in log
Gamess reads Stone's DMA analysis
============================================================================
"""
    if   type.lower() == 'slv':
         pass
         
    elif type.lower() == 'gamess' or type.lower() == 'gms':
         data = open(file)
         # ----------------------------------
         #querry1 = "PROPERTIES FOR THE B3LYP    DFT FUNCTIONAL (RHF  TYPE) DENSITY MATRIX"
         #querry2 = "MP2 PROPERTIES...FOR THE FIRST ORDER WAVEFUNCTION"
         line = data.readline()
         #while 1:
         #      if ((querry1 in line) or (querry2 in line)): break
         #      line = data.readline()
         # ----------------------------------
         querry = " NET CHARGES AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         ZerothMoments = []
         Structure = []
         while line.split()!=[]:
               ZerothMoments.append( numpy.float64( line.split()[2] ) )
               Structure.append( numpy.array( line.split()[-3:],dtype=numpy.float64)  )
                    
               line = data.readline()
         # ----------------------------------
         querry = " FIRST MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map( numpy.float64, line.split()[1:]))
               line = data.readline()
         # ----------------------------------
         querry = " SECOND MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( numpy.float64, line.split()[1:] ))
               line = data.readline()
         # ----------------------------------
         querry = " THIRD MOMENTS AT POINTS"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         ThirdMoments = []
         while 'CPU' not in line.split():
               A = map( numpy.float64, line.split()[1:] )
               line = data.readline()
               B = map( numpy.float64, line.split() )
               ThirdMoments.append( A+B )
               line = data.readline()

         # add up the nuclear and electronic charges
         npoints = len(FirstMoments)
         natoms  = len(ZerothMoments) - npoints
         #if
         for i in range(natoms):
             ZerothMoments[i+natoms] += ZerothMoments[i]
         for i in range(natoms):
             ZerothMoments.pop(0)
             Structure.pop(0)

         data.close()
         return     DMA( q  =numpy.array(ZerothMoments)   ,
                         m  =numpy.array(FirstMoments )   ,
                         T  =numpy.array(SecondMoments)   , 
                         O  =numpy.array(ThirdMoments )   ,
                         pos=numpy.array(Structure)       ,
                         hexadecapoles=False              )

    # -----------------------------------------------------------------------------
    elif type.lower() == 'coulomb' or type.lower() == 'c':
         # determine if the file contains hexadecapole moments
         has_hexadecapoles = False
         data = open(file)
         if " Distributed fourth-order property" in data.read(): has_hexadecapoles = True
         data.close()

         # now read the file once again line by line
         data = open(file)
         line = data.readline()

         # start reading the sections
         querry = " Distributed zeroth-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(3): line = data.readline()
         ZerothMoments = []
         Structure = []
         Origin = []
         while line.split()!=[]:
               ZerothMoments.append( numpy.float64( line.split()[0] ) )
               line = data.readline()
         nfrag = len(ZerothMoments)
         # ----------------------------------
         querry = " Distributed first-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         FirstMoments = []
         while line.split()!=[]:
               FirstMoments.append( map( numpy.float64, line.split()[:]))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed second-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(4): line = data.readline()
         SecondMoments = []
         while line.split()!=[]:
               SecondMoments.append( map( numpy.float64, line.split()[:] ))
               line = data.readline()
         # ----------------------------------
         querry = " Distributed third-order property"
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(5): line = data.readline()
         ThirdMoments = []
         #while '-----' not in line or :
         for i in range(nfrag):
               A = map( numpy.float64, line.split()[:] )
               line = data.readline()
               B = map( numpy.float64, line.split()[:] )
               ThirdMoments.append( A+B )
               line = data.readline()
         # ----------------------------------
         FourthMoments = None
         if has_hexadecapoles and hexadecapoles:
            querry = " Distributed fourth-order property" 
            while 1:
                if querry in line: break
                line = data.readline()
            FourthMoments = []
            for i in range(6): line = data.readline() 
            for i in range(nfrag):
                A = map( numpy.float64, line.split()[:] )
                line = data.readline()
                B = map( numpy.float64, line.split()[:] )
                line = data.readline()
                C = map( numpy.float64, line.split()[:] )
                line = data.readline()
                FourthMoments.append( A+B+C )
            FourthMoments = numpy.array(FourthMoments)
         # ----------------------------------
         querry = " form@"
         while 1:
             if querry in line: break
             line = data.readline()
         typ = line.split()[0]
         if   typ == 'primitive': is_traceless = False
         elif typ == 'traceless': is_traceless = True
         # ------------------------------------
         querry = " Structure"
         struct = True
         while 1:
               if querry in line: break
               line = data.readline()
               if not line: 
                  struct = False
                  break

         atoms = []
         if struct:
            for i in range(3): line = data.readline()
            while line.split()!=[]:
                  coord = line.split()
                  atoms.append(units.Atom(coord[0]))
                  Structure.append( map( numpy.float64, coord[1:] ) )
                  line = data.readline()

         Structure = numpy.array(Structure,dtype=numpy.float64)                  
         
         querry = " Origins"
         origins = True
         while 1:
               if querry in line: break
               line = data.readline()
               if not line: 
                  origins = False
                  break              

         if origins:
            for i in range(3): line = data.readline()
            while line.split()!=[]:
                  coord = line.split()
                  Origin.append( map( numpy.float64, coord[1:] ) )
                  line = data.readline()  
            
            Origin = numpy.array(Origin   ,dtype=numpy.float64)

         else:
            Origin    = Structure.copy()
         
         data.close() 
         return     DMA( q=numpy.array(ZerothMoments)   ,
                         m=numpy.array(FirstMoments )   , 
                         T=numpy.array(SecondMoments)   ,
                         O=numpy.array(ThirdMoments )   , 
                         H=FourthMoments                ,
                         atoms=atoms              ,
                         pos=Structure            ,
                         origin=Origin,
                         traceless=is_traceless )
    # -----------------------------------------------------------------------------
    elif type.lower() == 'gaussian' or type.lower() == 'gau':
         data = open(file)
         # seek for atomic positions!
         querry = "Electrostatic Properties Using The SCF Density"
         line = data.readline()
         while 1:
               if querry in line: break
               line = data.readline()
               if line =='': raise Exception('No CHELPG population found!')
         for i in range(4): line = data.readline()
         Structure = []
         while ('Atomic Center' in line or 'Ghost Center' in line):
               Structure.append( numpy.array(line.split()[-3:],dtype=numpy.float64)  )
               line = data.readline()

         # seek for charges!
         querry = " Charge="
         while 1:
               if querry in line: break
               line = data.readline()
         for i in range(2): line = data.readline()
         ZerothMoments = []
         for i in range(len(Structure)):
             ZerothMoments.append( numpy.float64( line.split()[-1] ) )    
             line = data.readline()
        
         Result =     DMA(nfrag=len(Structure), hexadecapoles=False)
         Result.pos = numpy.array(Structure) * units.UNITS.AngstromToBohr
         Result.DMA[0] = numpy.array(ZerothMoments)
 
         data.close()
         return Result
     
def interchange(T,ind):
    """\
>>>>> Copied from LIBBBG.utilities <<<<<
interchange rows according to order list

Usage: ingerchange(numpy.array,order_list)

Returns: permuted numpy.array

Example: 

Assume we have initial matrix T. We want obtain T':

T =  [[ 1  6  8]         T'=  [[ 3 -4 -4]
      [ 2  5  7]               [ 7 -1  9]
      [ 3 -4 -4]      ind      [ 1  6  8]
      [ 4  1  0]     ---->     [ 4  1  0]
      [ 5 -7 -8]               [ 2  5  7]
      [ 6  6  2]               [ 5 -7 -8]
      [ 7 -1  9]]              [ 6  6  2]]
      
This is accomplished by run:
T_prime = interchange(T,ind=[3, 7, 1, 4, 2, 5, 6])
"""
    ind = numpy.array(ind)-1
    B = T.copy()
    for i,index in enumerate(ind):
        B[i] = T[index]
    return B
 
class DMA:
    """\
==============================LIBBBG:DMA================================
Represents the DMA distribution object. Inputs are q,m,T and O that mean
a set of distributed multipoles on n centers, ie.:                      
   Q = numpy.array([Q1,Q2, ... ,Qn])                                          
where Qi is i-th distributed tensor. The DMA class stores also vectors  
of position and origin of each of distributed center in a form of a     
NumPy numpy.array. As a zero object you can specify DMA(nfrag=<n>) after which
the DMA of n distributed centers (zero moments at each) will be created.
                                                                        
========================================================================
Usage:                                                                  
------------------------------------------------------------------------
DMA(nfrag=<n>)                  return zero DMA object with <n> centers 
<object>.contract(<list>)       contract the DMA <object> using <list>  
                                that specifies the indices of origins to
                                be saved. Others will be removed.       
<object>.MakeUa(<list>,change_origin=<True>,contract=<False>)           
                                create united atoms within the DMA.     
                                unset 'change_origin' option when using 
                                distributed charges                     
<object>.ChangeOrigin(new_origin_set=<0>,zero=<False>)                 
                                change the origins of distribited moments
                                to the new origin set specified by numpy.array
                                (in bohrs) of dimension (<n>,3). If you 
                                want to translate the centers to the    
                                origin of coordinate system specify only
                                zero=True.                              
<object>.ChangeUnits                                                    
<object>.MAKE_FULL()            create full DMA format                  
<object>.MakeTraceless()        transform the DMA quadrupoles and octu- 
                                poles to traceless form. Needs full for-
                                mat created first (above)               
<object>.Rotate(<rot>)          rotate the DMA about the rotation matrix
                                Needs full format created first (above) 
<object>.copy()                 return a deepcopy of <object>           
<object>.write(<file>)          write DMA distribution to a <file> in cf
<object>.interchange(<ind>)     interchange all memorials using ind list
------------------------------------------------------------------------
"get" methods:                                                          
------------------------------------------------------------------------
<object>.get_pos()              return atomic positions ndnumpy.array (in bohr)
<object>.get_origin()           return origins ndnumpy.array (in bohr)       
<object>.get_name()             return name of distribution             
<object>.get_nfrags()           return number of distributed sites      
<object>.get_natoms()           return number of atoms                  
------------------------------------------------------------------------
"set" methods:                                                          
------------------------------------------------------------------------
<object>.set_name(<name>)       set the name of distribution            
<object>.set_moments(charges=<None>,dipoles=<None>,                     
                     quadrupoles=<None>,octupoles=<None>)               
                                set moments as numpy numpy.arrays in reduced  
                                format                                  
<object>.set_structure(pos=<None>,origin=<None>,                        
                       atoms=<None>,equal=<False>)                      
                                set positions and/or origins and/or atom
                                symbols. If equal=<True> setting only po
                                sitions will set the same origins too   
                                without explicit call for origins       
------------------------------------------------------------------------
mathematical operations:                                                
------------------------------------------------------------------------
1) addition                                                             
   <object> + <number> = adds number to all moments (new <object>)      
   <object> + <object> = adds moments of two objects (new <object>)     
   <object> - <number> = substract number to all moments (new <object>)
   <object> - <object> = substract moments from each other (new <object>)
   <object> * <number> = multiply all moments by number (new <object>)  
   <object> / <number> = divide all moments by number (new <object>)    
   - <object>          = return negative moments (new <object>)         
   len(object)         = return number of distributed sites             
------------------------------------------------------------------------
========================================================================
"""
    
    # ---- CONSTRUCTORS

    def __init__(self,
                 file=None,
                 nfrag=None,
                 name='Untitled',
                 q=None,m=None,T=None,O=None,H=None,
                 pos=numpy.zeros((1,3),dtype=numpy.float64),
                 origin=numpy.zeros((1,3),dtype=numpy.float64),
                 atoms=None,
                 traceless=False,
                 hexadecapoles=None):

        # save the flag whether this object contains hexadecapoles        
        if hexadecapoles is not None:
           self.has_hexadecapoles = hexadecapoles
        else:
           if H is None: self.has_hexadecapoles = False
           else        : self.has_hexadecapoles = True

        # determine the path of DMA object construction
        if file is not None:
           self.__call__(file)

        elif nfrag is not None:
             q=numpy.zeros((nfrag),dtype=numpy.float64)
             m=numpy.zeros((nfrag,3),dtype=numpy.float64)
             T=numpy.zeros((nfrag,6),dtype=numpy.float64)
             O=numpy.zeros((nfrag,10),dtype=numpy.float64)
             if self.has_hexadecapoles: H=numpy.zeros((nfrag,15),dtype=numpy.float64)
             pos = numpy.zeros((nfrag,3),dtype=numpy.float64)
             origin = numpy.zeros((nfrag,3),dtype=numpy.float64)

        if file is None:
           # DMA distribution in reduced format (GAMESS-like)
           if self.has_hexadecapoles: self.DMA = [q,m,T,O,H]
           else:                      self.DMA = [q,m,T,O]
           #self.pos = numpy.zeros((nfrag,3),dtype=numpy.float64)
           #self.origin = numpy.zeros((nfrag,3),dtype=numpy.float64)
           # name
           self.name = name
           # origin
           self.origin = origin
           # position numpy.array
           self.pos = pos
           # number of distributed sites (fragments)
           self.nfrag = len(self.DMA[0])
           # atomic sites if any
           if atoms is None:
              #if len(atoms)>1: self.atoms = [Atom('X')]
              #else: self.atoms = [Atom('X')]*self.nfrag
              self.atoms = [units.Atom('X')]*self.nfrag
           else: self.atoms=atoms
           # if DMA FULL formatted memorial were created. Now it is not, so 'False'.
           self.full = False
           # if traceless forms were created. Now it is ordinary (primitive) DMA format so 'False'.
           self.is_traceless = traceless
 
    def __call__(self,file):
        """open DMA file"""
        dma = ParseDMA(file,'c')
        if dma.has_hexadecapoles:
           self.DMA = [dma[x] for x in [0,1,2,3,4]]
        else:
           self.DMA = [dma[x] for x in [0,1,2,3]]
        self.origin = dma.get_origin()
        self.pos   = dma.get_pos()
        self.name  = dma.get_name()
        self.nfrag = len(dma[0])
        self.atoms = dma.atoms
        self.full  = False
        self.is_traceless = dma.is_traceless
        self.has_hexadecapoles = dma.has_hexadecapoles
        return

    # ---- GET methods

    def get_charges      (self): return self[0]
    def get_dipoles      (self): return self[1]
    def get_quadrupoles  (self): return self[2]
    def get_octupoles    (self): return self[3]
    def get_hexadecapoles(self): return self[4]
    
    def get_pos(self):
        """return positions of atoms"""
        return self.pos

    def get_positions(self): return self.get_pos()
     
    def get_origin(self):
        """return origins of distributed multipoles"""
        return self.origin

    def get_origins(self): return self.get_origin()
    
    def get_name(self):
        """return the name of distribution"""
        return self.name
    
    def get_nfrags(self):
        """return number of distributed sites"""
        return self.nfrag

    def get_natoms(self):
        """return number of atoms"""
        return len(self.pos)

    def get_magnitudes(self):
        """
Return the magnitudes of distributed multipoles. Octupoles and hexadecapoles NOT implemented yet!

Warning: The 0 value is returned for octupoles. Hexadecapole magnitude is not returned at all yet!"""
        # charges
        c1 = self[0].copy()
        # dipoles
        c2 = numpy.sqrt ( ((self[1].copy())**2).sum(axis=1) )
        # quadrupoles and octupoles
        temp = self.copy()
        if not temp.traceless:
           temp.MAKE_FULL()       
           temp.MakeTraceless()
           temp.makeDMAfromFULL()
        v3 = temp[2]; v4 = temp[3]
        c3 = numpy.zeros(self.nfrag, numpy.float64)
        c4 = numpy.zeros(self.nfrag, numpy.float64)
        for i in range(self.nfrag):
            QXX = v3[i,0] 
            QYY = v3[i,1]
            QZZ = v3[i,2]
            QXY = v3[i,3]
            QXZ = v3[i,4]
            QYZ = v3[i,5]
            v = QZZ**2 + 4./3. * (QXY**2 + QXZ**2 + QYZ**2) + 1./3. * (QXX - QYY)**2
            c3[i] = numpy.sqrt(v)
        return c1, c2, c3, c4

    def get_const(self, origin=numpy.array([0.,0.,0.]) ):
        """
 Calculate the invariants of total multipole moments
 Returns the tuple of four invariant sets for charge, dipole, quadrupole
 and octupole moments:

 1) charge invariant: 
     total charge

 2) dipole invariant:
     norm of a dipole vector

 3) quadrupole invariants:
     trace(Q)
     0.5 * (trace(Q)^2 - trace(Q^2))
     det(Q)
 
     where Q is quadrupole moment (primitive form)

 4) octupole invariants:
     F1 = E_iik * E_ppk
     F2 = E_ijj * E_iqq
     F3 = E_ijk * E_ijk
     F4 = E_ijk * E_kij
     F5 = 0.5 * (e_kp * E_ijk * E_pji + e_ri * E_ijk * E_kjr)

     where E is octupole moment (primitive form)
     and e is 2-order antisymmetric (permutation) tensor

 Notes:
  1) invariants or 3-rd order Cartesian tensors:
        from. F. Ahmad, Arch. Mech. 63, 4, pp383-392 Warszawa 2011
        formulae on p.390 for F1-F5

  2) no invariants are computed for hexadecapole moments!
     """
        eij = numpy.zeros((3, 3), numpy.float64)
        eij[0,1] = eij[1,2] = eij[2,0] =  1.0
        eij[1,0] = eij[2,1] = eij[0,2] = -1.0
        
        eijk = numpy.zeros((3, 3, 3), numpy.float64)
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] =  1.0
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1.0
        tot = self.get_mult(origin)
        tot.MAKE_FULL()
        # invariants for dipole moment
        mu  = numpy.sqrt(numpy.sum(tot.DMA_FULL[2][0]**2))
        # invariants for quadrupole moment
        tr = tot.DMA_FULL[3][0].trace()
        tr2= (numpy.dot(tot.DMA_FULL[3][0],tot.DMA_FULL[3][0])).trace()
        qad_1 = tr
        qad_2 = 0.500*(tr*tr-tr2)
        qad_3 = numpy.linalg.det(tot.DMA_FULL[3][0])
        # invariants for octupole moment
        # from. F. Ahmad, Arch. Mech. 63, 4, pp383-392 Warszawa 2011
        # formulae from p.390 for F1-F5
        oct   = tot.DMA_FULL[4][0]
        oct_1 = numpy.dot(oct.trace(axis1=0,axis2=1),oct.trace(axis1=0,axis2=1))
        oct_2 = numpy.dot(oct.trace(axis1=1,axis2=2),oct.trace(axis1=1,axis2=2))
        oct_3 = numpy.tensordot(oct,oct,((0,1,2),(0,1,2)))
        oct_4 = numpy.tensordot(oct,oct,((0,1,2),(1,2,0)))
        A_kp  = numpy.tensordot(oct,oct,((0,1),(2,1)))
        B_ri  = numpy.tensordot(oct,oct,((1,2),(1,0)))
        oct_5 = 0.5*(numpy.sum(eij*A_kp,axis=None) + numpy.sum(eij.transpose()*B_ri,axis=None))
        #oct = 1./6. * (numpy.sum(tot.DMA_FULL[4][0]*eijk)) * eijk
        qad = numpy.array([qad_1, qad_2, qad_3],numpy.float64)
        oct = numpy.array([oct_1,oct_2,oct_3,oct_4,oct_5],numpy.float64)
        return tot.DMA_FULL[1][0], mu, qad, oct

    def get_mult(self, origin=numpy.zeros(3,dtype=numpy.float64), hexadecapoles=False):
        """Calculates overall primitive multipole moments from distributed multipoles.
This can be used for a test of the correctness of multipole analysis. The exact distributed
expansions should sum up to the total molecular expansion with respect to origin privided.
"""
        if hexadecapoles and not self.has_hexadecapoles:
           print(" WARNING: This DMA object does not contain hexadecapoles but they were requested. It will be ignored by Coulomb.")
           hexadecapoles = False
        overall = DMA(nfrag=1, hexadecapoles=hexadecapoles)
        #overall.name = 'Test of reproducing multipoles from charges [units: Debyes(*Angstroms^n)]'
        overall.name = 'Test of reproducing overall multipoles from CAMM/CBAMM/LMTP/DMA [units: A.U.]'
        overall.pos = numpy.zeros((1,3),dtype=numpy.float64)
        overall.pos[0] = origin
        tmp = self.copy()
        tmp.MAKE_FULL()
        tmp.ChangeOrigin(zero=1)
        
        ### make full format of DMA
        tmp.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = numpy.zeros((3),dtype=numpy.float64)
        quad = numpy.zeros((3,3),dtype=numpy.float64)
        oct  = numpy.zeros((3,3,3),dtype=numpy.float64)
        hex  = numpy.zeros((3,3,3,3),dtype=numpy.float64)
        for atom in range(self.nfrag):
            r     = tmp.origin[atom]### zmiana origin z pos!!!
            qatom = tmp.DMA[0][atom]
            ### calculate dipole moment
            mu   += tmp.DMA[1][atom]
            ### calculate quadrupole moment
            quad += tmp.DMA_FULL[3][atom]
            ### calculate octupole moment
            oct  += tmp.DMA_FULL[4][atom]
            ### calculate hexadecapole moment
            if hexadecapoles:
               hex  += tmp.DMA_FULL[5][atom]

        ### set the molecular moments into the DMA solvent object
        overall.DMA[0][0] = sum(tmp[0])
        overall.DMA[1][0] = mu
           
        overall.DMA[2][0,0] = quad[0,0]           
        overall.DMA[2][0,1] = quad[1,1]
        overall.DMA[2][0,2] = quad[2,2]
        overall.DMA[2][0,3] = quad[0,1]
        overall.DMA[2][0,4] = quad[0,2]
        overall.DMA[2][0,5] = quad[1,2] 
           
        overall.DMA[3][0,0] = oct[0,0,0]
        overall.DMA[3][0,1] = oct[1,1,1]
        overall.DMA[3][0,2] = oct[2,2,2]
        overall.DMA[3][0,3] = oct[0,0,1]
        overall.DMA[3][0,4] = oct[0,0,2]
        overall.DMA[3][0,5] = oct[0,1,1]
        overall.DMA[3][0,6] = oct[1,1,2]
        overall.DMA[3][0,7] = oct[0,2,2]
        overall.DMA[3][0,8] = oct[1,2,2]
        overall.DMA[3][0,9] = oct[0,1,2]

        if hexadecapoles:
           overall.DMA[4][0, 0]  = hex[0,0,0,0]
           overall.DMA[4][0, 1]  = hex[1,1,1,1]
           overall.DMA[4][0, 2]  = hex[2,2,2,2]
           overall.DMA[4][0, 3]  = hex[0,0,0,1]
           overall.DMA[4][0, 4]  = hex[0,0,0,2]
           overall.DMA[4][0, 5]  = hex[1,1,1,0]
           overall.DMA[4][0, 6]  = hex[1,1,1,2]
           overall.DMA[4][0, 7]  = hex[2,2,2,0]
           overall.DMA[4][0, 8]  = hex[2,2,2,1]
           overall.DMA[4][0, 9]  = hex[0,0,1,1]
           overall.DMA[4][0,10]  = hex[0,0,2,2]
           overall.DMA[4][0,11]  = hex[1,1,2,2]
           overall.DMA[4][0,12]  = hex[0,0,1,2]
           overall.DMA[4][0,13]  = hex[1,1,0,2]
           overall.DMA[4][0,14]  = hex[2,2,0,1]


        overall.MAKE_FULL()
        overall.ChangeOrigin(new_origin_set=numpy.array([origin],numpy.float64) )
        
        #overall.DMA[1] *= UNITS.BohrElectronToDebye
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
        #overall.DMA[3] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom**2
          
        return overall
       
    def get_mult_c(self, origin=numpy.zeros(3,dtype=numpy.float64), hexadecapoles=False):
        """Calculates overall primitive moments from charges.
This tool has a purpose of testing population analysis obtained by
fitting to the molecular ab initio potential or other methods.

For an example, refer to the construction of total molecular 
solvatochromic multipoles from fitted distributed charges:

H. Lee, J.-H. Choi and M. Cho, J. Chem. Phys. 137(11), 114307 (2012)
"""
        overall = DMA(nfrag=1)
        overall.name = 'Test of reproducing multipoles from charges'
        overall.pos = numpy.zeros((1,3),dtype=numpy.float64)
        overall.pos[0] = origin
        
        ### make full format of DMA
        tmp = self.copy()
        tmp.MAKE_FULL()
        
        ### compute molecular solvatochromic moments
        mu   = numpy.zeros((3),dtype=numpy.float64)
        quad = numpy.zeros((3,3),dtype=numpy.float64)
        oct  = numpy.zeros((3,3,3),dtype=numpy.float64)
        hex  = numpy.zeros((3,3,3,3),dtype=numpy.float64)
        for atom in range(self.nfrag):
            r     = tmp.pos[atom] - origin
            qatom = tmp.DMA[0][atom]
            ### calculate dipole moment
            mu   += qatom * r 
            ### calculate quadrupole moment
            quad += qatom * numpy.outer (r,r)
            ### calculate octupole moment
            oct  += qatom * numpy.outer( r, numpy.outer (r,r) ).reshape(3,3,3)
            ### calculate hexadecapole moment
            if hexadecapoles:
               hex  += qatom * numpy.outer( r, numpy.outer (r, numpy.outer(r,r) ) ).reshape(3,3,3,3)

        ### set the molecular moments into the DMA solvent object
        overall.DMA[0][0] = sum(tmp[0])
        overall.DMA[1][0] = mu
           
        overall.DMA[2][0,0] = quad[0,0]           
        overall.DMA[2][0,1] = quad[1,1]
        overall.DMA[2][0,2] = quad[2,2]
        overall.DMA[2][0,3] = quad[0,1]
        overall.DMA[2][0,4] = quad[0,2]
        overall.DMA[2][0,5] = quad[1,2] 
           
        overall.DMA[3][0,0] = oct[0,0,0]
        overall.DMA[3][0,1] = oct[1,1,1]
        overall.DMA[3][0,2] = oct[2,2,2]
        overall.DMA[3][0,3] = oct[0,0,1]
        overall.DMA[3][0,4] = oct[0,0,2]
        overall.DMA[3][0,5] = oct[0,1,1]
        overall.DMA[3][0,6] = oct[1,1,2]
        overall.DMA[3][0,7] = oct[0,2,2]
        overall.DMA[3][0,8] = oct[1,2,2]
        overall.DMA[3][0,9] = oct[0,1,2]

        if hexadecapoles:
           overall.DMA[4][0, 0]  = hex[0,0,0,0]  
           overall.DMA[4][0, 1]  = hex[1,1,1,1]
           overall.DMA[4][0, 2]  = hex[2,2,2,2]
           overall.DMA[4][0, 3]  = hex[0,0,0,1]
           overall.DMA[4][0, 4]  = hex[0,0,0,2]
           overall.DMA[4][0, 5]  = hex[1,1,1,0]
           overall.DMA[4][0, 6]  = hex[1,1,1,2]
           overall.DMA[4][0, 7]  = hex[2,2,2,0]
           overall.DMA[4][0, 8]  = hex[2,2,2,1]
           overall.DMA[4][0, 9]  = hex[0,0,1,1]
           overall.DMA[4][0,10]  = hex[0,0,2,2]
           overall.DMA[4][0,11]  = hex[1,1,2,2]
           overall.DMA[4][0,12]  = hex[0,0,1,2]
           overall.DMA[4][0,13]  = hex[1,1,0,2]
           overall.DMA[4][0,14]  = hex[2,2,0,1]

        #overall.DMA[1] *= UNITS.BohrElectronToDebye
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
        #overall.DMA[2] *= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom**2
          
        return overall
    
    
    # ---- SET methods

    def set_name(self,name):
        """sets the name for the object""" 
        self.name = name
        return
    
    def set_structure(self, pos=None, origin=None, atoms=None, equal=False):
        """sets new positions or/and origins and atoms"""
        # update positions
        if pos is not None:
           if len(self.pos) != len(pos): # it means you have to update self.atoms
              self.atoms = [ units.Atom('X') for x in range(len(pos)) ]
           self.pos = pos.copy()
        # update origins
        if origin is not None:
           self.origin = origin.copy()
        # update atoms
        if atoms is not None:
           if isinstance(atoms[0], units.Atom): 
              self.atoms = atoms
           else:
              self.atoms = [ units.Atom(x) for x in atoms.split(',') ]
        # equal the positions and origins
        if equal:
           self.origin = pos.copy()

    def set_charges       (self, charges       ): self.set_moments(charges       = charges       )
    def set_dipoles       (self, dipoles       ): self.set_moments(dipoles       = dipoles       )
    def set_quadrupoles   (self, quadrupoles   ): self.set_moments(quadrupoles   = quadrupoles   )
    def set_octupoles     (self, octupoles     ): self.set_moments(octupoles     = octupoles     )
    def set_hexadecapoles (self, hexadecapoles ): self.set_moments(hexadecapoles = hexadecapoles )

    def set_moments(self,charges=None,dipoles=None,
                         quadrupoles=None,octupoles=None,
                         hexadecapoles=None):
        """\
set multipoles given as n-d numpy.array
where n is rank of multipole moment"""
        if charges       is not None: self.DMA[0] = numpy.array(charges       ).copy()
        if dipoles       is not None: self.DMA[1] = numpy.array(dipoles       ).copy()
        if quadrupoles   is not None: self.DMA[2] = numpy.array(quadrupoles   ).copy()
        if octupoles     is not None: self.DMA[3] = numpy.array(octupoles     ).copy()
        if hexadecapoles is not None: self.DMA[4] = numpy.array(hexadecapoles ).copy()
        return

    def set_moments_all(self, other):
        """Copy the moments from other to self"""
        o = other.copy()
        self.set_moments(charges=o.get_charges(),
                         dipoles=o.get_dipoles(),
                         quadrupoles=o.get_quadrupoles(),
                         octupoles=o.get_octupoles())
        if self.has_hexadecapoles and other.has_hexadecapoles:
           self.set_moments(hexadecapoles=o.get_hexadecapoles())
        elif self.has_hexadecapoles and not other.has_hexadecapoles:
           raise TypeError, " %s object does not have hexadecapoles so they cannot be copied to %s object!" % (other.get_name(), self.get_name())
        return


    # ---- property descriptor methods
 
    def if_traceless(self):
        """Is the dma object traceless?"""
        return self.is_traceless

    def if_hexadecapoles(self):
        """Does this object contain hexadecapole moments?"""
        return self.has_hexadecapoles

    # --- auxiliary methods
     
    def interchange(self,ind):
        """
premute all moments, positions and origins using ind list 
(ind has the same meaning as in interchange function"""
        n_range = 4
        if self.has_hexadecapoles: n_range = 5
        for i in range(n_range):
            self.DMA[i]  = interchange(self.DMA[i],ind)
        self.pos = interchange(self.pos,ind)
        self.origin = interchange(self.origin,ind)
        return
        
    def copy(self):
        """creates deepcopy of the object"""
        return copy.deepcopy(self)

    def contract(self,contrlist):
        """shrink the dimensions of DMA object by eliminating rows in attributes"""
        K = len(contrlist)
        self.nfrag = K
        self.DMA[0] = self.DMA[0][contrlist]
        self.DMA[1] = self.DMA[1][contrlist]
        self.DMA[2] = self.DMA[2][contrlist]
        self.DMA[3] = self.DMA[3][contrlist]
        if self.has_hexadecapoles:
           self.DMA[4] = self.DMA[4][contrlist]
        self.origin = self.origin[contrlist]
        return

    def concatenate(self, other):
        """add multipole sites to the self.DMA object. The sites of other.DMA will be appended at the end"""
        assert isinstance(other, DMA), " Argument has to be a DMA instance!"
        assert (self.is_traceless and other.is_traceless or not self.is_traceless and not other.is_traceless), \
                " The objects have to have the same forms!(primitive or traceless)"

        N_moments = 4
        if self.has_hexadecapoles and other.has_hexadecapoles: N_moments = 5

        for i in range(N_moments):
            self.DMA[i] = numpy.concatenate((self.DMA[i],other.DMA[i].copy()))
        self.pos        = numpy.concatenate((self.pos, other.pos.copy()))
        self.origin     = numpy.concatenate((self.origin, other.origin.copy()))
        self.nfrag      = len(self.DMA[0])
        self.atoms     += other.atoms
        return

    def replace(self, other):
        """Copy the contents from one DMA instance to another"""
        o = other.copy()
        self.set_structure(o.get_origin(), equal=True)
        self.atoms = o.atoms
        self.set_moments_all(o)
        return

    def translate(self,vector):
        """translate origins and positions by a vector"""
        self.origin += vector
        self.pos    += vector
        return

    def rotate(self, rot):
        """Rotate the DMA instance by rotation matrix """
        self.MAKE_FULL()
        self.Rotate(rot)
        return

    def traceless(self):
        """Transform the DMA instance into traceless form. 
This is an irreversible operation and primitive moments cannot be retrieved from traceless form."""
        self.MAKE_FULL()
        self.MakeTraceless()
        self.makeDMAfromFULL()
        return

    def trac(self): self.traceless()

    def sup(self, xyz, suplist=None, sup_origin=False):
        """Superimpose DMA position onto xyz. 

Usage:

  rms = DMA_instance.sup(self, xyz, suplist=None, sup_origin=False)

By default, superimposition is done based on atomic positions. 
Superimposition based on DMA origins can be set by setting sup_origin=True. 
Subset of superimposed centers can be chosen by specifying suplist=[...] 
indices (using normal numbering scheme, not in Python convention!).
For example, suplist = [1,2,5] and sup_origin=False will choose
first, second and fifth atom while sup_origin=False corresponding DMA origins will
be chosen.

Returned value is an RMS of a superimposition and is given in Bohrs.
"""
        # convert suplist to Python-convention indices
        if suplist is not None: suplist = numpy.array(suplist)-1
        # begin superimposition
        s = SVDSuperimposer()
        xyz_prev = self.pos
        if sup_origin: xyz_prev = self.origin
        if len(xyz)!=0:
           if suplist is None: s.set( xyz, xyz_prev )
           else:               s.set( xyz[suplist], xyz_prev[suplist] )
           s.run()
           rms = s.get_rms()
           rot, transl = s.get_rotran()
        else:
           rot = numpy.identity(3,numpy.float64)
           transl = xyz - xyz_prev
           rms = 0.0
        self.rotate(rot)
        self.translate(transl)
        return rms

    def recenter(self, new_origin_set=0, zero=False):
        """Recenter the DMA distribution.

 Usage:  

       DMA_instance.recenter(new_origin_set=0, zero=False)
 
 Where:
 
       new_origin_set - is the set of new origins
       zero           - if True, new_origin_set is ignored
                        and all DMA will be recentered to the
                        point [0.0, 0.0, 0.0] (global origin)
"""
        self.MAKE_FULL()
        self.ChangeOrigin(new_origin_set, zero)
        return

    def MakeUa(self,ua_list,change_origin=True,contract=False):
        """\
Transforms the object to the contracted DMA form employing united atoms. 
Usage:
MakeUa( ua_list )
        
where 'ua_list' is a list of tuples ti containing integers:
        
ua+list = [ t1, t2, t3 , ... ], 
ti = ( A, a1, a2, a3, ... )
        
For i-th tuple corresponding to i-th UA the first integer A
is the atom id of UA center atom and the remaining a1, a2 etc relate to 
the atoms supposed to be contracted within a united atom UA for A-th center.
The numbers are normal numbers (not in Python convention)."""
       
        if self.has_hexadecapoles: rs = range(5)
        else:                      rs = range(4) 
        origin = self.origin
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(zero=True)
        for t in ua_list:
            for i in t[1:]:
                for c in rs:
                    self.DMA[c][t[0]-1]+= self.DMA[c][i-1]
                    if not c: self.DMA[c][i-1] = 0
                    else:     self.DMA[c][i-1].fill(0)
        if change_origin:
           self.MAKE_FULL()
           self.ChangeOrigin(new_origin_set=origin)
      
        if contract:
           contrlist = self.__contrListFromUaList(ua_list)
           self.contract(contrlist)

    def ChangeUnits(self,charges,dipoles,quadrupoles,octupoles,hexadecapoles=1.0):
        """changes the units"""
        self.DMA[0] *= charges
        self.DMA[1] *= dipoles
        self.DMA[2] *= quadrupoles
        self.DMA[3] *= octupoles
        if self.has_hexadecapoles:
           self.DMA[4] *= hexadecapoles
        return

    def write(self, file, type='c'):
        """writes  the DMA distribution in a file or in a XYZ file (structure)"""
        newfile = open(file,'w')
        if type=='c' or type=='dma':
           log = ' SLV output file. All units in AU\n\n'
           log+= self.__repr__()
           log += ' Structure\n'
           log += ' ---------\n\n'
           for i,atom in enumerate(self.atoms):
               log+= " %-10s %14.8f %14.8f %14.8f\n" % (atom.symbol, self.pos[i,0], self.pos[i,1], self.pos[i,2])
           log+='\n'
        
           # print origins in case the where different that structure
           if (self.pos.shape != self.origin.shape or not (self.pos==self.origin).all() ):
              log += ' Origins\n'
              log += ' -------\n\n'
              for i,point in enumerate(self.origin):
                  log+= " %-10s %14.8f %14.8f %14.8f\n" % ('X', point[0], point[1], point[2])
              log+='\n'
        
        elif type=='xyz':
           log = "  %i\n\n"%len(self.pos)
           pos = self.pos * units.UNITS.BohrToAngstrom
           for i,atom in enumerate(self.atoms):
              log+= " %-10s %14.8f %14.8f %14.8f\n" % (atom.symbol, pos[i,0], pos[i,1], pos[i,2])
           log+='\n'
           
        newfile.write(log)
        newfile.close() 
        return          

    # PRIVATE METHODS

    def __contrListFromUaList(self,ua_list):
       """creates the contraction list from ua_list"""
       return []
            
    def __getitem__(self,index): 
        return numpy.array(self.DMA[index])

    def __setitem__(self,index,value):
        self.DMA[index] = numpy.array(value)
    
    def __len__(self):
        """Return the number of distributed sites"""
        n_range = 4
        if self.has_hexadecapoles: n_range = 5
        length = len(self.DMA[0])
        for i in range(n_range-1):
            assert length == len(self.DMA[i]), "The number of distributed %i-th poles is not equal to number of 0-poles (charges)!" % (i+1)
        return(length)

    ### --- mathematical operations
    ###     Warning! If between 2 DMA objects they have to have
    ###     same dimensions!
    
    def __add__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] + other
           m = self.DMA[1] + other
           T = self.DMA[2] + other
           O = self.DMA[3] + other
           if self.has_hexadecapoles: 
              H = self.DMA[4] + other
              case =-1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] + self
           m = other.DMA[1] + self
           T = other.DMA[2] + self
           O = other.DMA[3] + self
           if other.has_hexadecapoles: 
              H = other.DMA[4] + self
              case = -2
        else:
           case = 3
           q = other.DMA[0] + self.DMA[0]
           m = other.DMA[1] + self.DMA[1]
           T = other.DMA[2] + self.DMA[2]
           O = other.DMA[3] + self.DMA[3]     
           if self.has_hexadecapoles and other.has_hexadecapoles: 
              H = other.DMA[4] + self.DMA[4]
              case = -3

        # return DMA without hexadecapoles
        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())

        # return DMA with hexadecapoles
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __sub__(self,other):
        q = self.DMA[0] - other.DMA[0]
        m = self.DMA[1] - other.DMA[1]
        T = self.DMA[2] - other.DMA[2]
        O = self.DMA[3] - other.DMA[3]
        if self.has_hexadecapoles and other.has_hexadecapoles:
           H = self.DMA[4] - other.DMA[4]
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())

    
    def __div__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] / other
           m = self.DMA[1] / other
           T = self.DMA[2] / other
           O = self.DMA[3] / other
           if self.has_hexadecapoles:
              H = self.DMA[4] / other 
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           q = other.DMA[0] / self
           m = other.DMA[1] / self
           T = other.DMA[2] / self
           O = other.DMA[3] / self
           if other.has_hexadecapoles:
              H = other.DMA[4] / self
              case = -2
        else:
           q = self.DMA[0] / other.DMA[0]
           m = self.DMA[1] / other.DMA[1]
           T = self.DMA[2] / other.DMA[2]
           O = self.DMA[3] / other.DMA[3]     
           if self.has_hexadecapoles and other.has_hexadecapoles:
              H = self.DMA[4] / other.DMA[4]
              case = -3

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())


    def __mul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
           if self.has_hexadecapoles:
              H = self.DMA[4] * other
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
           if other.has_hexadecapoles:
              H = other.DMA[4] * self
              case = -2
        else:
           case = 3
           q = self.DMA[0] * other.DMA[0]
           m = self.DMA[1] * other.DMA[1]
           T = self.DMA[2] * other.DMA[2]
           O = self.DMA[3] * other.DMA[3]            
           if self.has_hexadecapoles and other.has_hexadecapoles:
              H = self.DMA[4] * other.DMA[4]
              case = -3

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __rmul__(self,other):
        if (isinstance(self,DMA) and (not isinstance(other,DMA))):
           case = 1
           q = self.DMA[0] * other
           m = self.DMA[1] * other
           T = self.DMA[2] * other
           O = self.DMA[3] * other
           if self.has_hexadecapoles:
              H = self.DMA[4] * other
              case = -1
        elif ((not isinstance(self,DMA)) and isinstance(other,DMA)):
           case = 2
           q = other.DMA[0] * self
           m = other.DMA[1] * self
           T = other.DMA[2] * self
           O = other.DMA[3] * self
           if other.has_hexadecapoles:
              H = other.DMA[4] * self
              case = -2

        if case > 0:
           return DMA(q=q,m=m,T=T,O=O,
                      pos=self.pos.copy(),origin=self.origin.copy())
        else:
           return DMA(q=q,m=m,T=T,O=O,H=H,
                      pos=self.pos.copy(),origin=self.origin.copy())

    def __neg__(self):
        q = -self.DMA[0]
        m = -self.DMA[1]
        T = -self.DMA[2]
        O = -self.DMA[3]
        if self.has_hexadecapoles:
          H = -self.DMA[4]
          return DMA(q=q,m=m,T=T,O=O,H=H,
                     pos=self.pos.copy(),origin=self.origin.copy())
        else:
          return DMA(q=q,m=m,T=T,O=O,
                     pos=self.pos.copy(),origin=self.origin.copy())

    def __repr__(self):
        """PRINTOUT FORM FOR DMA REDUCED FORMAT"""
        log = "\n"
        log+= " --- %s ---\n\n" % (self.name)
        log+= " Distributed zeroth-order property\n"
        log+= " ---------------------------------\n"
        log+= "\n"
        charges = self.DMA[0]
        for xfrag in charges:
            log+= "%16.6E\n"%xfrag
        log+= "\n"

        log+= " Distributed first-order property \n"
        log+= " ---------------------------------\n"
        log+= "\n"
        dipoles = self.DMA[1]
        log+= "%16s %16s %16s\n" % ('X'.rjust(16),'Y'.rjust(16),'Z'.rjust(16))
        for xfrag in dipoles:
            log+= "%16.6E %16.6E %16.6E\n"%(xfrag[0],xfrag[1],xfrag[2])
        log+= "\n"  

        log+= " Distributed second-order property\n"
        log+= " ---------------------------------\n"
        log+= "\n"
        quadrupoles = self.DMA[2]
        log+= "%16s %16s %16s %16s %16s %16s\n" %\
             ('XX'.rjust(16),'YY'.rjust(16),'ZZ'.rjust(16),
              'XY'.rjust(16),'XZ'.rjust(16),'YZ'.rjust(16))
        for xfrag in quadrupoles:
            log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
             (xfrag[0],xfrag[1],xfrag[2],
              xfrag[3],xfrag[4],xfrag[5])
        log+= "\n"

        log+= " Distributed third-order property \n"
        log+= " ---------------------------------\n"
        log+= "\n"
        octupoles = self.DMA[3]
        log+= "%16s %16s %16s %16s %16s %16s\n" %\
            ('XXX'.rjust(16),'YYY'.rjust(16),'ZZZ'.rjust(16),
             'XXY'.rjust(16),'XXZ'.rjust(16),'XYY'.rjust(16))
        log+= "%16s %16s %16s %16s\n" %\
            ('YYZ'.rjust(16),'XZZ'.rjust(16),'YZZ'.rjust(16),'XYZ'.rjust(16))
        for xfrag in octupoles:
            log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n" %\
                   (xfrag[0],xfrag[1],xfrag[2],
                    xfrag[3],xfrag[4],xfrag[5])
            log+= "%16.6E %16.6E %16.6E %16.6E\n" %\
                   (xfrag[6],xfrag[7],xfrag[8],xfrag[9])

        if self.has_hexadecapoles:
           log+= "\n"
           log+= " Distributed fourth-order property\n"              
           log+= " ---------------------------------\n"
           log+= "\n"
           hexadecapoles = self.DMA[4]
           log+= "%16s %16s %16s\n" % ('XXXX'.rjust(16),'YYYY'.rjust(16),'ZZZZ'.rjust(16))
           log+= "%16s %16s %16s %16s %16s %16s\n" %\
                ('XXXY'.rjust(16),'XXXZ'.rjust(16),'YYYX'.rjust(16),
                 'YYYZ'.rjust(16),'ZZZX'.rjust(16),'ZZZY'.rjust(16))
           log+= "%16s %16s %16s %16s %16s %16s\n" %\
                ('XXYY'.rjust(16),'XXZZ'.rjust(16),'YYZZ'.rjust(16),
                 'XXYZ'.rjust(16),'YYXZ'.rjust(16),'ZZXY'.rjust(16))
           for xfrag in hexadecapoles:
               log+= "%16.6E %16.6E %16.6E\n"%(xfrag[0],xfrag[1],xfrag[2])
               log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
                (xfrag[3],xfrag[4],xfrag[5],
                 xfrag[6],xfrag[7],xfrag[8])
               log+= "%16.6E %16.6E %16.6E %16.6E %16.6E %16.6E\n"%\
                (xfrag[9 ],xfrag[10],xfrag[11],
                 xfrag[12],xfrag[13],xfrag[14])

        log = self.__end_multipole_section(log)
        
        return str(log) 

    def __end_multipole_section(self, log):
       log+= " "+"-"*100+"\n"                                                                               
       if self.nfrag == 1:
          if self.is_traceless: 
             log+= (" traceless form@,   origin: %s\n"  % str(self.origin*0.5291772086)[1:-1]).rjust(100)  
          else: 
             log+= (" primitive form@,   origin: %s \n" % str(self.origin*0.5291772086)[1:-1]).rjust(100) 
       else:
          if self.is_traceless: 
             log+= " traceless form@\n".rjust(100)  
          else: 
             log+= " primitive form@\n".rjust(100) 
       log+= "\n\n\n"
       return log


    # Intrinsic methods which will be converted into private in the future: user should not use them any longer explicitly

    def MAKE_FULL(self):
        """creates numpy.arrays of ordinary forms of distributed multipoles
        in full format."""
        CHARGES = numpy.array(self.DMA[0])
        DIPOLES = numpy.array(self.DMA[1])
        QDPOLES = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
        OCTPLES = numpy.zeros((self.nfrag,3,3,3),dtype=numpy.float64)
        if self.has_hexadecapoles:
           HEXPLES = numpy.zeros((self.nfrag,3,3,3,3),dtype=numpy.float64)
        for site in range(self.nfrag):

            # add quadrupoles
            q = numpy.array(self.DMA[2])
            QDPOLES[site,0,0] = q[site,0]
            QDPOLES[site,1,1] = q[site,1]
            QDPOLES[site,2,2] = q[site,2]
            QDPOLES[site,0,1] = q[site,3]
            QDPOLES[site,0,2] = q[site,4]
            QDPOLES[site,1,2] = q[site,5]
            #
            QDPOLES[site,1,0] = QDPOLES[site,0,1]
            QDPOLES[site,2,0] = QDPOLES[site,0,2]
            QDPOLES[site,2,1] = QDPOLES[site,1,2]
            
            # add octupoles
            o = numpy.array(self.DMA[3])
            OCTPLES[site,0,0,0] = o[site,0]
            OCTPLES[site,1,1,1] = o[site,1]
            OCTPLES[site,2,2,2] = o[site,2]
            OCTPLES[site,0,0,1] = o[site,3]
            OCTPLES[site,0,0,2] = o[site,4]
            OCTPLES[site,0,1,1] = o[site,5] 
            OCTPLES[site,1,1,2] = o[site,6]
            OCTPLES[site,0,2,2] = o[site,7]
            OCTPLES[site,1,2,2] = o[site,8]
            OCTPLES[site,0,1,2] = o[site,9]
            # 
            OCTPLES[site,0,1,0] = OCTPLES[site,0,0,1] 
            OCTPLES[site,1,0,0] = OCTPLES[site,0,0,1] 
            #
            OCTPLES[site,0,2,0] = OCTPLES[site,0,0,2] 
            OCTPLES[site,2,0,0] = OCTPLES[site,0,0,2] 
            #
            OCTPLES[site,1,0,1] = OCTPLES[site,0,1,1] 
            OCTPLES[site,1,1,0] = OCTPLES[site,0,1,1] 
            #
            OCTPLES[site,1,2,1] = OCTPLES[site,1,1,2]
            OCTPLES[site,2,1,1] = OCTPLES[site,1,1,2]
            #
            OCTPLES[site,2,0,2] = OCTPLES[site,0,2,2]
            OCTPLES[site,2,2,0] = OCTPLES[site,0,2,2]
            #
            OCTPLES[site,2,1,2] = OCTPLES[site,1,2,2]
            OCTPLES[site,2,2,1] = OCTPLES[site,1,2,2]
            #
            OCTPLES[site,0,2,1] = OCTPLES[site,0,1,2]
            OCTPLES[site,2,0,1] = OCTPLES[site,0,1,2]
            #
            OCTPLES[site,1,0,2] = OCTPLES[site,0,1,2]
            #
            OCTPLES[site,1,2,0] = OCTPLES[site,0,1,2]
            OCTPLES[site,2,1,0] = OCTPLES[site,0,1,2]
            
            # add hexadecapoles
            if self.has_hexadecapoles:
               h = numpy.array(self.DMA[4])    # REDUCED FORMAT COMPONENTS

               a = h[site,0]             # XXXX
               HEXPLES[site,0,0,0,0] = a
               
               a = h[site,1]             # YYYY
               HEXPLES[site,1,1,1,1] = a
               
               a = h[site,2]             # ZZZZ
               HEXPLES[site,2,2,2,2] = a
               
               a = h[site,3]             # XXXY
               HEXPLES[site,0,0,0,1] = a
               HEXPLES[site,0,0,1,0] = a
               HEXPLES[site,0,1,0,0] = a
               HEXPLES[site,1,0,0,0] = a
               
               a = h[site,4]             # XXXZ
               HEXPLES[site,0,0,0,2] = a
               HEXPLES[site,0,0,2,0] = a
               HEXPLES[site,0,2,0,0] = a
               HEXPLES[site,2,0,0,0] = a
               
               a = h[site,5]             # YYYX
               HEXPLES[site,0,1,1,1] = a
               HEXPLES[site,1,0,1,1] = a
               HEXPLES[site,1,1,0,1] = a
               HEXPLES[site,1,1,1,0] = a
               
               a = h[site,6]             # YYYZ
               HEXPLES[site,1,1,1,2] = a
               HEXPLES[site,1,1,2,1] = a
               HEXPLES[site,1,2,1,1] = a
               HEXPLES[site,2,1,1,1] = a
               
               a = h[site,7]             # ZZZX
               HEXPLES[site,0,2,2,2] = a
               HEXPLES[site,2,0,2,2] = a
               HEXPLES[site,2,2,0,2] = a
               HEXPLES[site,2,2,2,0] = a
               
               a = h[site,8]             # ZZZY
               HEXPLES[site,1,2,2,2] = a
               HEXPLES[site,2,1,2,2] = a
               HEXPLES[site,2,2,1,2] = a
               HEXPLES[site,2,2,2,1] = a
               
               a = h[site,9]             # XXYY
               HEXPLES[site,0,0,1,1] = a
               HEXPLES[site,0,1,0,1] = a
               HEXPLES[site,0,1,1,0] = a
               HEXPLES[site,1,0,0,1] = a
               HEXPLES[site,1,0,1,0] = a
               HEXPLES[site,1,1,0,0] = a
               
               a = h[site,10]            # XXZZ
               HEXPLES[site,0,0,2,2] = a
               HEXPLES[site,0,2,0,2] = a
               HEXPLES[site,0,2,2,0] = a
               HEXPLES[site,2,0,0,2] = a
               HEXPLES[site,2,0,2,0] = a
               HEXPLES[site,2,2,0,0] = a
               
               a = h[site,11]            # YYZZ
               HEXPLES[site,1,1,2,2] = a
               HEXPLES[site,1,2,1,2] = a
               HEXPLES[site,1,2,2,1] = a
               HEXPLES[site,2,1,1,2] = a
               HEXPLES[site,2,1,2,1] = a
               HEXPLES[site,2,2,1,1] = a
               
               a = h[site,12]            # XXYZ
               HEXPLES[site,0,0,1,2] = a
               HEXPLES[site,0,0,2,1] = a
               HEXPLES[site,0,1,0,2] = a
               HEXPLES[site,0,1,2,0] = a
               HEXPLES[site,0,2,0,1] = a
               HEXPLES[site,0,2,1,0] = a
               HEXPLES[site,1,0,0,2] = a
               HEXPLES[site,1,0,2,0] = a
               HEXPLES[site,1,2,0,0] = a
               HEXPLES[site,2,0,0,1] = a
               HEXPLES[site,2,0,1,0] = a
               HEXPLES[site,2,1,0,0] = a
               
               a = h[site,13]            # YYXZ
               HEXPLES[site,0,1,1,2] = a
               HEXPLES[site,0,1,2,1] = a
               HEXPLES[site,0,2,1,1] = a
               HEXPLES[site,1,0,1,2] = a
               HEXPLES[site,1,0,2,1] = a
               HEXPLES[site,1,1,0,2] = a
               HEXPLES[site,1,1,2,0] = a
               HEXPLES[site,1,2,0,1] = a
               HEXPLES[site,1,2,1,0] = a
               HEXPLES[site,2,0,1,1] = a
               HEXPLES[site,2,1,0,1] = a
               HEXPLES[site,2,1,1,0] = a
               
               a = h[site,14]            # ZZXY
               HEXPLES[site,0,1,2,2] = a
               HEXPLES[site,0,2,1,2] = a
               HEXPLES[site,0,2,2,1] = a
               HEXPLES[site,1,0,2,2] = a
               HEXPLES[site,1,2,0,2] = a
               HEXPLES[site,1,2,2,0] = a
               HEXPLES[site,2,0,1,2] = a
               HEXPLES[site,2,0,2,1] = a
               HEXPLES[site,2,1,0,2] = a
               HEXPLES[site,2,1,2,0] = a
               HEXPLES[site,2,2,0,1] = a
               HEXPLES[site,2,2,1,0] = a
        
        if self.has_hexadecapoles:    
           self.DMA_FULL = [ self.origin , CHARGES , DIPOLES , QDPOLES , OCTPLES , HEXPLES ]
        else:
           self.DMA_FULL = [ self.origin , CHARGES , DIPOLES , QDPOLES , OCTPLES ]
        self.full = True

    def makeDMAfromFULL(self):
        """Saves the reduced-formatted DMA data from full-formatted data"""
        self.pos         = numpy.array(self.DMA_FULL[0])
        self.DMA[0]      = numpy.array(self.DMA_FULL[1])
        self.DMA[1]      = numpy.array(self.DMA_FULL[2])
        #
        self.DMA[2][:,0] = numpy.array(self.DMA_FULL[3][:,0,0])
        self.DMA[2][:,1] = numpy.array(self.DMA_FULL[3][:,1,1])
        self.DMA[2][:,2] = numpy.array(self.DMA_FULL[3][:,2,2])
        self.DMA[2][:,3] = numpy.array(self.DMA_FULL[3][:,0,1])
        self.DMA[2][:,4] = numpy.array(self.DMA_FULL[3][:,0,2])
        self.DMA[2][:,5] = numpy.array(self.DMA_FULL[3][:,1,2])
        #
        self.DMA[3][:,0] = numpy.array(self.DMA_FULL[4][:,0,0,0])
        self.DMA[3][:,1] = numpy.array(self.DMA_FULL[4][:,1,1,1])
        self.DMA[3][:,2] = numpy.array(self.DMA_FULL[4][:,2,2,2])
        self.DMA[3][:,3] = numpy.array(self.DMA_FULL[4][:,0,0,1])
        self.DMA[3][:,4] = numpy.array(self.DMA_FULL[4][:,0,0,2])
        self.DMA[3][:,5] = numpy.array(self.DMA_FULL[4][:,0,1,1])
        #
        self.DMA[3][:,6] = numpy.array(self.DMA_FULL[4][:,1,1,2])
        self.DMA[3][:,7] = numpy.array(self.DMA_FULL[4][:,0,2,2])
        self.DMA[3][:,8] = numpy.array(self.DMA_FULL[4][:,1,2,2])
        self.DMA[3][:,9] = numpy.array(self.DMA_FULL[4][:,0,1,2])
        #
        if self.has_hexadecapoles:
           self.DMA[4][:, 0]  = numpy.array(self.DMA_FULL[5][:,0,0,0,0])  
           self.DMA[4][:, 1]  = numpy.array(self.DMA_FULL[5][:,1,1,1,1])
           self.DMA[4][:, 2]  = numpy.array(self.DMA_FULL[5][:,2,2,2,2])
           self.DMA[4][:, 3]  = numpy.array(self.DMA_FULL[5][:,0,0,0,1])
           self.DMA[4][:, 4]  = numpy.array(self.DMA_FULL[5][:,0,0,0,2])
           self.DMA[4][:, 5]  = numpy.array(self.DMA_FULL[5][:,1,1,1,0])
           self.DMA[4][:, 6]  = numpy.array(self.DMA_FULL[5][:,1,1,1,2])
           self.DMA[4][:, 7]  = numpy.array(self.DMA_FULL[5][:,2,2,2,0])
           self.DMA[4][:, 8]  = numpy.array(self.DMA_FULL[5][:,2,2,2,1])
           self.DMA[4][:, 9]  = numpy.array(self.DMA_FULL[5][:,0,0,1,1])
           self.DMA[4][:,10]  = numpy.array(self.DMA_FULL[5][:,0,0,2,2])
           self.DMA[4][:,11]  = numpy.array(self.DMA_FULL[5][:,1,1,2,2])
           self.DMA[4][:,12]  = numpy.array(self.DMA_FULL[5][:,0,0,1,2])
           self.DMA[4][:,13]  = numpy.array(self.DMA_FULL[5][:,1,1,0,2])
           self.DMA[4][:,14]  = numpy.array(self.DMA_FULL[5][:,2,2,0,1])
 
        del self.DMA_FULL
        self.full = False
 
    def MakeTraceless_wrong_code(self):
        """turns ordinary full-formatted DMA_FULL into traceless 
full-formatted DMA_FULL. The traceless form is taken to be
following Buckingham convention (0.5 * Jackson convention) /add citation!!!"""

        if self.is_traceless:
           raise  Exception("\n Error: DMA is already traceless!!! quitting ...\n")

        if self.full:
           P3 = [(0,1,2,3),(0,2,1,3),(0,3,1,2)]                                            # permutations for oct
           P4 = [(0,1,2,3,4),(0,1,3,2,4),(0,1,4,2,3),(0,2,3,1,4),(0,2,4,1,3),(0,3,4,1,2)]  # permutations for hex
           
           Q_P = self.DMA_FULL[3]
           O_P = self.DMA_FULL[4]
 
           # generate vector of identities
           d = numpy.identity(3,dtype=numpy.float64)
           I3= numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           for i in range(self.nfrag): I3[i] = d.copy()

           # compute traceless properties 
           Q_T = (3./2.) * Q_P - 0.500 * utilities2.array_outer_product_1_n( Q_P.trace(axis1=1,axis2=2), I3 )
           O_T = (5./2.) * O_P
           for perm in P3: O_T -=  0.500 * (utilities2.array_outer_product_1_2(O_P.trace(axis1=perm[2],axis2=perm[3]), I3)).transpose(perm)

           self.DMA_FULL[3] = Q_T
           self.DMA_FULL[4] = O_T

           if self.has_hexadecapoles:
              H_P = self.DMA_FULL[5]
              H_T = (35./8.) * H_P
              H_trace = H_P.trace(axis1=1,axis2=2)
              I33  = utilities2.array_outer_product_2_2(I3, I3)
              I33 += I33.transpose((0,1,3,2,4)) + I33.transpose((0,1,4,2,3))
              B = utilities2.array_outer_product_2_2(H_trace, I3)
              for perm in P4: H_T -= (5./8.) * B.transpose(perm)
              H_T += utilities2.array_outer_product_1_n( H_trace.trace(axis1=1,axis2=2) , I33) / 8.00

              self.DMA_FULL[5] = H_T

           self.is_traceless = True
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")
        return

    def MakeTraceless(self):
        """turns ordinary full-formatted DMA_FULL into traceless 
        full-formatted DMA_FULL. The traceless form is taken to be
        following Buckingham /add citation!!!"""
        
        if self.is_traceless:
           raise  Exception("\nerror: DMA is already traceless!!! quitting ...\n")
           
        if self.full:
           # Quadrupole moment
           for Q in self.DMA_FULL[3]:
               t = Q.trace()
               Q *= (3./2.)
               for i in [0,1,2]: Q[i,i]-=(1./2.)*t
           # Octapole moment
           for O in self.DMA_FULL[4]:
               W  = O.copy()
               O *= (5./2.)
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       for k in [0,1,2]:
                           Wt = W.trace()
                           if   i==j and i==k:
                                O[i,j,k]-= (1./2.) * (Wt[i] + Wt[j] + Wt[k])
                           elif i==j and i!=k:
                                O[i,j,k]-= (1./2.) * Wt[k]
                           elif j!=k and i==k:
                                O[i,j,k]-= (1./2.) * Wt[j]
                           elif i!=j and j==k:
                                O[i,j,k]-= (1./2.) * Wt[i]

           # Hexadecapole moment                                                                                       
           if self.has_hexadecapoles:
              for H in self.DMA_FULL[5]:
                  W = H.copy()
                  H *= (35./8.)
                  Wt = W.trace()
                  Wtt= W.trace().trace()
                  c = (5./8.)
                  p = (1./8.)
                  d = numpy.identity(3,numpy.float64)
                  for i in [0,1,2]:
                      for j in [0,1,2]:
                          for k in [0,1,2]:
                              for l in [0,1,2]:
                                  H[i,j,k,l] -= c * (d[i,j] * Wt[k,l] + d[i,k] * Wt[j,l] + d[i,l] * Wt[j,k] + 
                                                     d[j,k] * Wt[i,l] + d[j,l] * Wt[i,k] + d[k,l] * Wt[i,j] ) \
                                             -  p * Wtt * (d[i,j] * d[k,l] + d[i,k] * d[j,l] + d[i,l] * d[j,k] )
                            
           self.is_traceless = True
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def ChangeOrigin(self,new_origin_set=0,zero=False):
        """shifts origin(s) of ordinary full-formatted DMA_FULL
           into new one(s)"""

        if self.full:
           if zero:
              new_origin_set = numpy.zeros((self.nfrag,3),dtype=numpy.float64)
           old_origin_set = numpy.array(self.origin)#numpy.array(self.DMA_FULL[0])
           M0 = numpy.array(self.DMA_FULL[1])
           M1 = numpy.array(self.DMA_FULL[2]) 
           M2 = numpy.array(self.DMA_FULL[3]) 
           M3 = numpy.array(self.DMA_FULL[4])

           #A2 = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           #B2 = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           A2 = utilities2.array_outer_product(old_origin_set,old_origin_set)
           B2 = utilities2.array_outer_product(new_origin_set,new_origin_set)

           D = new_origin_set - old_origin_set
           
           M1_o = M1 + utilities2.array_outer_product_1_n(M0, old_origin_set)
           #AM = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           #MA = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           AM = utilities2.array_outer_product(old_origin_set,M1_o)
           MA = utilities2.array_outer_product(M1_o,old_origin_set)
           M2_o = M2 - utilities2.array_outer_product_1_n(M0, A2) + AM + MA
           
           # First moments
           new_M1 = M1 - utilities2.array_outer_product_1_n(M0, new_origin_set- old_origin_set)
           # Second moments
           #MD = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           MD = utilities2.array_outer_product(M1_o,D)
           #DM = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           DM =  utilities2.array_outer_product(D,M1_o)
           new_M2 = M2 + utilities2.array_outer_product_1_n(M0,(B2-A2)) - MD - DM

           # Third moments
           A3 = utilities2.array_outer_product_1_2(old_origin_set,A2)#
           B3 = utilities2.array_outer_product_1_2(new_origin_set,B2)#
           M1A= utilities2.array_outer_product(M1_o,old_origin_set)#
           M1B= utilities2.array_outer_product(M1_o,new_origin_set)#
           #BAM2 = numpy.zeros((self.nfrag,3,3,3),dtype=numpy.float64)
           BAM2  = utilities2.array_outer_product_1_2(D,M2_o)#
           B2A2M1= utilities2.array_outer_product_2_1(B2-A2,M1_o)#
           QB3A3 = utilities2.array_outer_product_1_n(M0,(B3-A3))#
           M2BA  = utilities2.array_outer_product_2_1(M2_o,D)#
           M1B2A2= utilities2.array_outer_product_1_2(M1_o,B2-A2)#
           AM1A  = utilities2.array_outer_product_1_2(old_origin_set,M1A)#
           BM1B  = utilities2.array_outer_product_1_2(new_origin_set,M1B)#
           #M1AM1 = utilities2.array_outer_product_2_1(M1A,old_origin_set)
           #M1BM1 = utilities2.array_outer_product_2_1(M1B,new_origin_set)
           BM2TR = numpy.transpose(utilities2.array_outer_product_1_2(new_origin_set,M2_o),
                                   axes=(0,2,1,3))#
           AM2TR = numpy.transpose(utilities2.array_outer_product_1_2(old_origin_set,M2_o),
                                   axes=(0,2,1,3))#

           new_M3 = M3   - BAM2 + B2A2M1 - QB3A3 - M2BA + M1B2A2 + \
                    BM1B - AM1A - BM2TR  + AM2TR
           
           # Fourth moments
           if self.has_hexadecapoles:
              M4 = numpy.array(self.DMA_FULL[5])   # hexadecapoles centered at old origins 
              #raise NotImplementedError, 'Recentering hexadecapoles is not implemented yet! Quitting...'
              # first find the dipoles, quadrupoles and octupoles at the origin [0,0,0]
              other = self.copy()
              other.has_hexadecapoles = False
              other.MAKE_FULL()
              other.ChangeOrigin(zero=True)
              other.MAKE_FULL()
              M10 = other.DMA_FULL[2]; M20 = other.DMA_FULL[3]; M30 = other.DMA_FULL[4]
              del other
              # next define the permutation sequences for generalized outer products [ M_k x D^(4-k) ] where k is rank of multipole at the origin [0,0,0]
              T1 = [(0,1,2,3,4),(0,2,1,3,4),(0,2,3,1,4),(0,2,3,4,1)]                          # dipole permutations
              T2 = [(0,1,2,3,4),(0,1,3,2,4),(0,3,1,4,2),(0,3,4,1,2),(0,1,3,4,2),(0,3,1,2,4)]  # quadrupole permutations
              T3 = [(0,1,2,3,4),(0,1,2,4,3),(0,1,4,2,3),(0,4,1,2,3)]                          # octupole permutations
              # then construct the basis generalized outer tensor products
              D2 = B2 - A2
              D3 = B3 - A3
              D4 = utilities2.array_outer_product_2_2(B2,B2) - utilities2.array_outer_product_2_2(A2,A2)
              Y1 = utilities2.array_outer_product_1_3(M10,D3)
              Y2 = utilities2.array_outer_product_2_2(M20,D2)
              Y3 = utilities2.array_outer_product_3_1(M30,D )
              # finally, compute the new hexadecapoles from lower moments
              new_M4 = M4 + utilities2.array_outer_product_1_n(M0,D4)
              for perm in T1: new_M4 -= Y1.transpose(perm)
              for perm in T2: new_M4 += Y2.transpose(perm)
              for perm in T3: new_M4 -= Y3.transpose(perm)


           if self.nfrag == 1:
              self.origin = new_origin_set[0]
         
           #self.DMA_FULL.append(new_origin_set)
           self.DMA_FULL[2] = new_M1
           self.DMA_FULL[3] = new_M2
           self.DMA_FULL[4] = new_M3
           if self.has_hexadecapoles:
              self.DMA_FULL[5] = new_M4
           
           self.makeDMAfromFULL()
           self.origin = numpy.array(new_origin_set)
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n")

    def Rotate(self,rotmat):
        """rotates the ordinary full-formatted DMA_FULL
           in molecular space based on rotation matrix
           given"""

        #self.MAKE_FULL()
        if self.full:
           ### transform the origins and positions
           for i in range(len(self.origin)):
               self.origin[i] = numpy.dot(self.origin[i],rotmat)
           for i in range(len(self.pos)):
               self.pos[i] = numpy.dot(self.pos[i],rotmat)
           ### transform the dipoles
           self.DMA_FULL[2] = numpy.dot(self.DMA[1],rotmat)

           ### transform the quadrupoles
           quadrupole = numpy.zeros((self.nfrag,3,3),dtype=numpy.float64)
           for i in range(self.nfrag):
               quadrupole[i] = numpy.dot(rotmat.T,numpy.dot( self.DMA_FULL[3][i] ,rotmat))

           self.DMA_FULL[3] = quadrupole
        
           ### transform the octupoles
           octupole = numpy.zeros((self.nfrag,3,3,3),dtype=numpy.float64)
           for i in range(self.nfrag):
               octupole[i] = numpy.tensordot(rotmat.T,numpy.tensordot(
                                       rotmat.T,numpy.tensordot(
                                       rotmat.T,self.DMA_FULL[4][i] , (1,2)
                                     ),(1,2)
                                     ),(1,2)
                                     )

           self.DMA_FULL[4] = octupole

           ### transform the hexadecapoles
           if self.has_hexadecapoles:
              hexadecapole = numpy.zeros((self.nfrag,3,3,3,3),dtype=numpy.float64)                       
              for i in range(self.nfrag):
                  hexadecapole[i] = numpy.tensordot(rotmat.T,numpy.tensordot(
                                              rotmat.T,numpy.tensordot(
                                              rotmat.T,numpy.tensordot(
                                              rotmat.T,self.DMA_FULL[5][i] , (1,3)
                                            ),(1,3)
                                            ),(1,3)
                                            ),(1,3)
                                            )

              self.DMA_FULL[5] = hexadecapole

           self.makeDMAfromFULL()
           
        else: raise Exception("\nerror: no FULL DMA object created! quitting...\n") 

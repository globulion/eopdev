#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------- #
#             EFFECTIVE FRAGMENT PARAMETER FORMAT MODULE                #
# --------------------------------------------------------------------- #

import sys, copy, os, re, math, numpy #, libbbg, solvshift.efprot, \
#       PyQuante.Ints
from ._libbbg_units import UNITS
from ..math.matrix import rotate_ao_matrix, Superimposer
import psi4
sys.stdout.flush()

__all__     = ['Fragment',]
__version__ = '2.0.1'

class Fragment(UNITS):
    """
 =============================================================================
                 Effective Fragment Parameter Format System                   
                                                                              
            Bartosz Błasiak                last revision: 21 Oct 2014
                                           adaptation to DMS and Psi4:
                                                          12 Mar 2020
 =============================================================================
                                                                              
 Represents Efective Fragment Potential (EFP) fragment. Designed to unify func
 tionality of the QM molecular object that contains information about all inte
 rmolecular interactions based on gas-phase molecular properties, mainly first
 order Coulomb electrostatics, second-order polarization forces, first-order c
 orrection due to overlap and Pauli repulsion, charge-transfer and dispersion.
 Supports the following functionalities for condensed phase physics modeling: 
 
   a) intermolecular interaction energies                                     
   b) nonlinear electrooptical properties                                     
   c) infra-red solvatochromism                                               
 
 It contains various memory sections necessary for evaluation of molecular agg
 regate property using coarse-grained models. The major groups of molecular pa
 rameters are:                                                                
 
   I  ) Basic molecular specification                                         
        a) atomic coordinates, numbers and masses                             
 
   II ) Electrostatic population parameters                                   
        a) distributed multipole moments                                      
        b) distributed dipole-dipole polarizabilities                         
 
   III) Nonlinear electrooptical properties                                   
        a) molecular charge and dipole moment                                 
        b) molecular polarizability                                           
        c) molecular first- and second-order hyperpolarizabilities            
 
   IV ) Frequency analysis data                                               
        a) harmonic frequencies, reduced masses and mass-weighted eigenvectors
        b) cubic anharmonic force constants                                   
 
   V  ) Wave function and QM energetics (CMO and LMO space)                   
        a) wave-function coefficients                                         
        b) LMO centroids                                                      
        c) Fock matrix                                                        
        d) derivatives of the above wrt electric field and normal coordinates 
 
 -----------------------------------------------------------------------------
 All of the above mentioned quantities are stored in ATOMIC UNITS!!!          
 -----------------------------------------------------------------------------
 
 Beneath, I list the further contents of this documentation box:              
    * How to use Frag class                                                   
    * Listing of memorials stored in the Frag instance                        
    * Information about storage dimensions of Frag memorials                  
 
 =============================================================================
  [1] USAGE (if '=': the argument is default)                     
 =============================================================================
 
  1) create the object (initialize)                                           
                                                                             
     frag = Frag()       initialize empty object or ...
     frag = Frag(name)   where name is the name of molecule or fragment writte
                         n in SLV standard fragment library $SLV_PATH/frg
                         The following are possible:
    -----------------------------------------------------------------------------------------------
       name             molecule                        EINT   FREQ   NLO      comment                           
    -----------------------------------------------------------------------------------------------                   
     o water            H2O                             YES    NO     NO       with DMA-5                        
     o water2           H2O                             YES    NO     NO       with CAMM-3                       
     o dmso             DMSO                            YES    NO     NO                                         
     o meoh             Methanol                        YES    NO     NO                                         
     o etoh             Ethanol                         YES    NO     NO
     o mesh             Methanothiol                    YES    NO     NO
     o dms              Dimethylsylfide                 YES    NO     NO
     o menh2            Methylamine                     YES    NO     NO
     o ccl4             Tetrachloromethane              YES    NO     NO                                         
     o chcl3            Chloroform                      YES    NO     NO                                         
     o dcm              Dichloromethane                 YES    NO     NO  
     o li+              Lithium Cation                  YES    NO     NO
     o na+              Sodium Cation                   YES    NO     NO
     o cl-              Chloride Anion                  YES    NO     NO
     o me-so3-          Methyl Sulphonate (-1)          YES    NO     NO
     o so3--            SO3 (2-) anion                  YES    NO     NO
     o meoac            Methyl acetate                  YES    YES    NO       C=O stretch                       
     o nma              N-Methyl Acetamide              YES    YES    NO       amide I mode                      
     o nma-d7           NMA-D7                          YES    YES    NO       amide I' mode (all H deuterated)
     o mescn            Methyl Thiocyanate              YES    YES    NO       C=N stretch
     o me-s-13-15n      Methyl Thiocyanate              YES    YES    NO       13C=15N (nitrile)
     o mecn             Methyl Nitrile                  YES    YES    NO       C=N stretch
     o et2coo           Diethyl Carbonate               YES    YES    NO       C=O stretch
     o imidazole        Imidazole                       YES    NO     NO
     o imidazolium+     Imidazolium Cation              YES    NO     NO
     o 4-me-imidazole   4-Methyl Imidazole              YES    NO     NO
     o 4-me-phenol      4-Methyl Phenol                 YES    NO     NO
     o chonh2           Formamide                       YES    NO     NO
     o chonhme          N-Methyl Formamide              YES    NO     NO
     o comenh2          Acetamide                       YES    NO     NO
     o mecooh           Acetic Acid                     YES    NO     NO
     o mecoo-           Acetate (-1) Anion              YES    NO     NO
     o menh3+           Methyl Ammonium (+1) Cation     YES    NO     NO
     o me-guanidinium+  Methyl Guanidinium (+1) Cation  YES    NO     NO
     o methane          Methane                         YES    NO     NO
     o ethane           Ethane                          YES    NO     NO
     o n-propane        n-Propane                       YES    NO     NO
     o cyclohexane      Cyclohexane                     YES    NO     NO
     o benzene          Benzene                         YES    NO     NO
     o phenole          Phenole                         YES    NO     NO
     o phelonate-       Phenolate Anion                 YES    NO     NO
     o dmf              Dimethylformamide               YES    NO     NO
     o thf              Tetrahydrofurane                YES    NO     NO
     o cho-ch-nh2-ch3   2-Aminopropanal                 YES    NO     NO
     o me-s-co-ch-ch2   S-Methyl-2-Propenethioate       YES    NO     NO
     o me-s-cho         Thiomethylaldehyde              YES    NO     NO
    -----------------------------------------------------------------------------------------------
     frag = Frag(file)   where file is the Frag format file (*.frg) with parame
                         ters created by a user. See slv_der-* for fragment par
                         ameter generators
     frag(file)          after initialization of Frag instance read the *.frg file
                                                                              
  2) set the data to frag from other source than Frag format file (*.frg)     
                                                                              
     frag.set(mol=None, anh=None, frag=None, dma=None, chelpg=None, esp=None) 
                         mol  - Molecule    class object                      
                         anh  - FREQ        class object                      
                         frag - FragFactory class object                      
                         dma  - DMA         class object                      
                         and other:                                           
                         * elpg - array of ChelpG charges                     
                         * esp    - array of ESP    charges                   
                         * other keywors that match list of memorials (see [2])
                                                                              
  3) write the parameters to a file                                           
                                                                              
     frag.write(file='slv.frg')                                               
                                                                              
  4) actions on Frag object                                                   
                                                                              
     frag.translate(transl)     # translate all tensors by cartesian vector   
     frag.rotate(rot)           # rotate all tensors by unitary 3x3 matrix    
     frag.sup(xyz)              # superimpose (translate + rotate) to the pos
                                # indicated by xyz Cartesian coordinates
     frag.reorder(ord)          # reorder the turn of atoms by specifying ord
                                # vector, e.g.: ord=[3,1,2] is equivalent to
                                #        A-B-C ----> B-C-A
                                #        1 2 3
     f_copy = frag.copy()       # create a deep copy of frag objects
     print frag                 # print the status of frag object

  5) return objects from frag object

     par = frag.get()           # get parameter dictionary. The keys are the
                                # 'Shortcut' names in section [1] of this docu
                                # mentation, e.g.:
     vecl1 = par['vecl1']       # gets LMO wave function first derivatives wrt
                                # all normal coordinates
     name= frag.get_name()      # return the name of fragment 
     pos = frag.get_pos()       # get position Cartesian coordinates of atoms
     bfs = frag.get_bfs()       # get BasisSet object                        
     qad, oct = frag.get_traceless() # get traceless quadrupoles and octupoles
                                                                              
 =============================================================================
  [2] LISTING OF MEMORIAL NAMES                                               
      The *.frg file has two parts:                                           
      * preambule (molecule section)                                          
      * parameters sections                                                   
 =============================================================================
  GR : Section name                      Shortcut       Dimension             
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : molecule                    · · · mol                                  
     :                             ·      name                                
     :                             ·      shortname                           
     :                             ·      basis                               
     :                             ·      charge
     :                             P      multiplicity                       
     :                             R      nbasis                              
     :                             E      method                              
     :                             A      atoms                               
     :                             M      natoms                              
  I  :                             B      ndma                                
     :                             U      npol                                
     :                             L      nmos                                
     :                             E      ncmos
     :                             ·      naocc                               
     :                             ·      nbocc                               
     :                             ·      navir                               
     :                             ·      nbvir                               
     :                             ·      mode
     :                             · · ·  nmodes                              
     : Atomic symbols                    atoms           natoms               
     : Atomic numbers                    atno            natoms               
     : Atomic masses                     atms            natoms               
     : Atomic coordinates                pos             natoms, 3            
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : ESP charges                       esp             natoms               
     : ChelpG charges                    chlpg           natoms               
     : DMTP centers                      rdma            ndma  , 3            
     : DMTP charges                      dmac            ndma                 
     : DMTP dipoles                      dmad            ndma  , 3            
 II  : DMTP quadrupoles                  dmaq            ndma  , 6            
     : DMTP octupoles                    dmao            ndma  , 10           
     : DMTP charges fder                 dmac1           nmodes, ndma         
     : DMTP dipoles fder                 dmad1           nmodes, ndma , 3      
     : DMTP quadrupoles fder             dmaq1           nmodes, ndma , 6      
     : DMTP octupoles fder               dmao1           nmodes, ndma , 10     
     : DMTP charges sder                 dmac2           nmodes_sder, ndma                 
     : DMTP dipoles sder                 dmad2           nmodes_sder, ndma  , 3              
     : DMTP quadrupoles sder             dmad2           nmodes_sder, ndma  , 6              
     : DMTP octupoles sder               dmao2           nmodes_sder, ndma  , 10             
     : Polarizable centers               rpol            npol  , 3            
     : Distributed polarizabilities      dpol            npol  , 9            
     : Distributed polarizabilities      dpol1           nmodes, npol, 9      
     : - first derivatives                                                    
     : Dist. pol. wrt imaginary freq     dpoli           npol  , 12  , 9
     : Dist. pol. wrt imaginary freq 
     : - first derivatives               dpoli1          nmodes, npol, 12, 9
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : Harmonic frequencies              freq            nmodes               
     : Reduced masses                    redmass         nmodes               
 IV  : Mass-weighted eigenvectors        lvec            nmodes, natoms, 3    
     : Cubic anharmonic constants        gijk            nmodes, nmodes, nmodes
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : LMO centroids                     lmoc            nmos  , 3            
     : LMO centroids                     lmoc1           nmodes, nmos  , 3    
     : - first derivatives                                                    
     : Fock matrix                       fock            nmos  , nmos         
     : Fock matrix                       fock1           nmodes, nmos  , nmos 
     : - first derivatives                                                    
     : Canonical Fock matrix             fckc            ncmos , ncmos        
     : Canonical Fock matrix             fckc1           nmodes, ncmos , ncmos
  V  : - first derivatives                                                    
     : AO->LMO matrix                    vecl            nmos  , nbasis       
     : AO->LMO matrix                    vecl1           nmodes, nmos  , nbasis
     : - first derivatives                                                    
     : AO->CMO matrix                    vecc            ncmos , nbasis       
     : AO->CMO matrix                    vecc1           nmodes, ncmos , nbasis
     : - first derivatives                                                    
     : One-particle density matrix       opdm            nbasis, nbasis
     : AO->CMO Occ matrix - Alpha        caocc           nbasis, naocc
     : AO->CMO Vir matrix - Alpha        cavir           nbasis, navir
 =============================================================================
  [3] TENSOR DIMENSIONS AND STORAGE CONVENTIONS                               
 =============================================================================
   All of the nonsymmetric tensors exept for DMA and DPOL analyses are in full 
   format. DMA and DPOL are stored in reduced formats. Other symmetric tensor 
   elements are stored in triangular matrix format (Fock matrices + derivs, cu
   bic anharmonic constants).

   Reduced format for electrostatic descriptors:
                                                                              
             CHG(i) .
                    1
             DIP(i) X   Y   Z
                    1   2   3
             QAD(i) XX  YY  ZZ  XY  XZ  YZ
                    1   2   3   4   5   6
             OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
                    1   2   3   4   5   6   7   8   9   10
             POL(i) XX  XY  XZ  YX  YY  YZ  ZX  ZY  ZZ 
                    1   2   3   4   5   6   7   8   9
                                                                              
   * Turn of indices in the case of DMA distributions is identical as GAMESS  
   * Turn of indices in distributed polarizabilities is different as in GAMESS
 =============================================================================
 Notes:                                                                       
   1) -ALL- the data in the instance of Frag are stored in ATOMIC UNITS!!!    
                                                                              
       [ Angstrom, kcal/mole, cm-1, C, amu, ... ] ----> [a.u.]                
                                                                              
 References:                                                                  
  * add Gordon's papers, M.Cho's papers, R.W.Góra's papers and your papers    
  * add MDAnalysis reference                                                  
  * add LibBBG, Coulomb.py, PyQuante, Gaussian09 and GAMESS references        
 =============================================================================
  Globulion@                                                                  
"""
    ### built-in parameters
    params = ['water'  ,'water2' ,'nma'     ,'nma-d7'   ,'mescn'   ,'meoac'   ,
              'meoh'   ,'dmso'   ,'dcm'     ,'chcl3'    ,'na+'     ,'so3--'   ,
              'me-so3-','li+'    ,'et2coo'  ,'mecn'     ,'4-me-imidazole'     , 
              '4-me-phenol'      ,'chonh2'  ,'chonhme'  ,'comenh2' ,'ethane'  ,
              'mecoo-' ,'menh3+' ,'methane' ,'n-propane','me-guanidinium+'    ,
              'ccl4'   ,'benzene','dmf'     ,'thf'      ,'cyclohexane'        ,
              'etoh'   ,'mesh'   ,'menh2'   ,'phenol'   ,'cho-ch-nh2-ch3'     , 
              'dms'    ,'imidazole'         ,'mecooh'   ,'me-s-co-ch-ch2'     ,
              'cl-'    ,'imidazolium+'      ,'me-s-cho' ,'phenolate-'         ,
              'me-s-13c-15n'     ,] 
               
    def __init__(self, file=None):
        self.__file = file
        self.__dir  = None
        self._create()
        #self._nlines = lambda n: n/5+bool(n%5)
        if file is not None: self.__call__(file)
    
    def _nlines(self,n):
        """auxiliary method"""
        return int(n/5+bool(n%5))
 
    def __call__(self, file):
        """parse the parameters from a file"""
        if file in self.params:
           self.__file = os.environ['SLV_DATA'] + '/frg/%s/%s.frg' % (file, file)
           self.__dir  = os.environ['SLV_DATA'] + '/frg/%s/'       %  file 
           if not os.path.isfile(self.__file):
              message = " No such file: <%s>. There is no *.frg file for that molecule,\n" % self.__file
              message+= " the name is misspelled or check if the $SLV_DATA variable is set properly\n"
              message+= " (echo $SLV_DATA)" 
              raise IOError(message)
        else: self.__file = file
        filef = open(self.__file,'r')
        text = filef.read()
        filef.close()
        templ = re.compile(r'\[',re.DOTALL)
        sections = re.split(templ,text)[1:]
        del text
        for section in sections:
            self._read(section)
        return
    
    def __repr__(self):
        """print the information about the status"""
        log = ' \n'
        log+= ' ================================= \n'
        log+= ' SLV SOLVATOCHROMIC EFP PARAMETERS \n'
        log+= ' ================================= \n\n'
        log+= ' MOLECULE SPECIFICATION \n'
        log+= ' name     %s    \n' %  self.__name
        log+= ' method   %s/%s \n' % (self.__method, self.__basis)
        log+= ' charge   %s    \n' %  self.__charge
        log+= ' multipli %s    \n' %  self.__multiplicity
        log+= ' natoms   %s    \n' %  self.__natoms
        log+= ' nbasis   %s    \n' %  self.__nbasis
        log+= ' naocc    %s    \n' %  self.__naocc
        log+= ' nbocc    %s    \n' %  self.__nbocc
        log+= ' nmodes   %s    \n' %  self.__nmodes
        log+= ' \n'
        return str(log)
   
    # --- public
    # SET METHODS
    
    def set(self, mol=None, method=None, basis=None, anh=None, frag=None, dma=None, dms=None, **kwargs):
        """set the parameters providing appropriate objects"""
        # molecular structure
        if mol is not None:
           self.__name   = mol.name()
           self.__charge = mol.molecular_charge()
           self.__multiplicity = mol.multiplicity()
           self.__natoms = mol.natom()
           self.__nmodes = 3 * self.__natoms - 6 # assumed non-linear molecule
           self.__pos    = mol.geometry().to_array(dense=True)
           self.__atoms  = numpy.array([mol.symbol(x) for x in range(self.__natoms)])
           self.__atno   = numpy.array([mol.Z     (x) for x in range(self.__natoms)])
           self.__atms   = numpy.array([mol.mass  (x) for x in range(self.__natoms)]) * self.AmuToElectronMass
           self.__nmos   = int((self.__atno.sum() - self.__charge)/2)  # assumed closed-shell
           self.__naocc  = self.__nmos
           self.__nbocc  = self.__nmos
           self.__npol   = self.__nmos
        # level of theoretical description
        if method is not None:
           self.__method = method
        # AO primary basis set
        if basis is not None:
           self.__basis  = basis
           self.__bfs    = self.get_bfs()
           self.__nbasis = self.__bfs.nbf()
           self.__ncmos  = self.__nbasis
           self.__navir  = self.__ncmos - self.__naocc
           self.__nbvir  = self.__ncmos - self.__nbocc
        # electrostatic data
        if dma is not None:
           assert not dma.if_traceless(), "DMA object is ALREADY traceless!!!"
           self.__pos    = dma.get_pos()
           self.__rdma   = dma.get_origin()
           self.__ndma   = len(self.__pos)
           self.__dmac   = dma[0]
           self.__dmad   = dma[1]
           self.__dmaq   = dma[2]
           self.__dmao   = dma[3]
        # anharmonic file object (FREQ)
        if anh is not None:
           assert anh.if_w(), 'Anharmonic object is in wrong units! Supply anh.w() object'
           self.__redmass = anh.redmass
           self.__freq = anh.freq
           self.__lvec = self._tr_lvec(anh.L, self.__nmodes, self.__natoms)
           self.__gijk = anh.K3
        # DMS object
        if dms is not None:
           self.__opdm = dms.M()
           self.__caocc= dms.Ca_occ()
           self.__cavir= dms.Ca_vir()
           if (1,0) in dms.available_orders(): self.__dms_10 = dms.B(1,0)
           if (2,0) in dms.available_orders(): self.__dms_20 = dms.B(2,0)
           if (0,1) in dms.available_orders(): self.__dms_01 = dms.B(0,1)
           if (0,2) in dms.available_orders(): self.__dms_02 = dms.B(0,2)
           if (0,3) in dms.available_orders(): self.__dms_03 = dms.B(0,3)
           if (0,4) in dms.available_orders(): self.__dms_04 = dms.B(0,4)
           if (0,5) in dms.available_orders(): self.__dms_05 = dms.B(0,5)
           if (0,6) in dms.available_orders(): self.__dms_06 = dms.B(0,6)
           if (1,1) in dms.available_orders(): self.__dms_11 = dms.B(1,1)
           if (2,1) in dms.available_orders(): self.__dms_21 = dms.B(2,1)
           if (1,2) in dms.available_orders(): self.__dms_12 = dms.B(1,2)
           if (2,2) in dms.available_orders(): self.__dms_22 = dms.B(2,2)

        # EFP fragment parameters
        if frag is not None:
           self._parse_dict( frag.get() )
        # Other special data
        if kwargs: self.set_par(**kwargs)
        return
    
    def set_par(self,**kwargs):
        """set the particular parameters"""
        for key, arg in kwargs.items():
            if key == 'rdma'   : self.__rdma  = arg      
            if key == 'rpol'   : self.__rpol  = arg
            if key == 'npol'   : self.__npol  = arg
            if key == 'dmac'   : self.__dmac  = arg
            if key == 'dmad'   : self.__dmad  = arg
            if key == 'dmaq'   : self.__dmaq  = arg
            if key == 'dmao'   : self.__dmao  = arg
            if key == 'dpol'   : self.__dpol  = arg
            if key == 'dpol1'  : self.__dpol1 = arg
            if key == 'dpoli'  : self.__dpoli = arg
            if key == 'dpoli1' : self.__dpoli1= arg
            if key == 'fock'   : self.__fock  = arg
            if key == 'fock1'  : self.__fock1 = arg
            if key == 'vecl'   : self.__vecl  = arg
            if key == 'vecl1'  : self.__vecl1 = arg
            if key == 'lmoc'   : self.__lmoc  = arg
            if key == 'lmoc1'  : self.__lmoc1 = arg
            if key == 'chlpg'  : self.__chlpg = arg
            if key == 'esp'    : self.__esp   = arg
            if key == 'dmac1'  : self.__dmac1 = arg
            if key == 'dmad1'  : self.__dmad1 = arg
            if key == 'dmaq1'  : self.__dmaq1 = arg
            if key == 'dmao1'  : self.__dmao1 = arg
            if key == 'dmac2'  : self.__dmac2 = arg
            if key == 'dmad2'  : self.__dmad2 = arg
            if key == 'dmaq2'  : self.__dmaq2 = arg
            if key == 'dmao2'  : self.__dmao2 = arg
            if key == 'lvec'   : self.__lvec  = arg
            if key == 'freq'   : self.__freq  = arg
            if key == 'redmass': self.__redmass = arg
            if key == 'redmss' : self.__redmass = arg
            if key == 'gijk'   : self.__gijk  = gijk
            if key == 'opdm'   : self.__opdm  = arg
            if key == 'caocc'  : self.__caocc = arg
            if key == 'cavir'  : self.__cavir = arg
            if key == 'dms_10' : self.__dms_10= arg
            if key == 'dms_20' : self.__dms_20= arg
            if key == 'dms_01' : self.__dms_01= arg
            if key == 'dms_02' : self.__dms_02= arg
            if key == 'dms_03' : self.__dms_03= arg
            if key == 'dms_04' : self.__dms_04= arg
            if key == 'dms_05' : self.__dms_05= arg
            if key == 'dms_06' : self.__dms_06= arg
            if key == 'dms_11' : self.__dms_11= arg
            if key == 'dms_12' : self.__dms_12= arg
            if key == 'dms_21' : self.__dms_21= arg
            if key == 'dms_22' : self.__dms_22= arg
            if key == 'mode'   : 
               self.__mode  = arg
               self.__nmodes_sder = len(arg)
        return

    # GET METHODS
            
    def get(self):
        """returns dictionary with parameters"""
        return self._make_dict()
   
    def get_dir(self):
        """Return the base directory of fragment storage"""
        return self.__dir
 
    def get_name(self):
        """return the name of the fragment"""
        return self.__name

    def get_nbf(self):
        """return the number of AO basis functions of the fragment"""
        return self.__nbasis
    
    def get_natoms(self):
        """return the number of atoms"""
        return self.__natoms

    def get_pos(self):
        """return position of a fragment"""
        return self.__pos

    def get_atno(self):
        """return atomic numbers"""
        return self.__atno

    def get_atms(self):
        """return atomic masses"""
        return self.__atms

    def _make_molecule(self):
        units = "Angs" # force Angstrom units
        c = 1.0; log_units = "units bohr\n"
        if units.lower().startswith("angs"):                                         
           c = psi4.constants.bohr2angstroms
           log_units = "units angstrom\n"
        log = "\n"
        log+= "%d %d\n" % (self.__charge, self.__multiplicity)
        for i in range(self.__natoms):
            log += "%s"     %  self.__atoms[i]
            log += "%16.6f" % (self.__pos[i,0] * c)
            log += "%16.6f" % (self.__pos[i,1] * c)
            log += "%16.6f" % (self.__pos[i,2] * c)
            log += "\n"
        log += log_units
        log += "symmetry c1\n"
        log += "no_reorient\n"
        log += "no_com\n"
        mol = psi4.geometry(log)
        return mol
   
    def get_mol(self):
        """return Psi4.core.Molecule object"""
        #mol = libbbg.utilities.MakeMol(self.__atno, self.__pos, name=self.__name)
        mol = self._make_molecule()
        return mol

    def get_bfs(self):
        """return basis set object"""
        #mol = libbbg.utilities.MakeMol(self.__atno, self.__pos, name=self.__name)
        #bfs = PyQuante.Ints.getbasis(mol, self.__basis)
        mol = self._make_molecule()
        bfs = psi4.core.BasisSet.build(mol, "BASIS", self.__basis, puream=False, quiet=True)
        return bfs

    def get_property(self, key, **kwargs):
        """
Evaluate various fragment-derived properties.
Available properties:

 Interaction energies
 o Dispersion            property='edisp6' , frag=<frg>

 Dispersion coefficients
 o C6    (LMO-distr)     property='c6d'    , frag=<frg>
 o C6    (total)         property='c6'     , frag=<frg>

 Solvatochromic properties
 o SolCAMM               property='solcamm', ea=True/False
 o SolC6 (LMO-distr)     property='solc6d' , frag=<frg>
 o SolC6 (total)         property='solc6'  , frag=<frg>
"""
        if key.lower() == 'c6d'    : prop = self._get_property_c6d    (**kwargs)
        if key.lower() == 'c6'     : prop = self._get_property_c6     (**kwargs)
        if key.lower() == 'edisp6' : prop = self._get_property_edisp6 (**kwargs)
        if key.lower() == 'solcamm': prop = self._get_property_solcamm(**kwargs)
        if key.lower() == 'solc6d' : prop = self._get_property_solc6d (**kwargs)
        if key.lower() == 'solc6'  : prop = self._get_property_solc6  (**kwargs)
        if key.lower() == 'gms_field': prop = self._get_property_gms_field(**kwargs) # TODO
        return prop
  
    def get_rotranrms(self):
        """Return rotation matrix, translation vector
and RMS from the last superimposition"""
        return self.__rot, self.__transl, self.__rms

    def get_traceless(self):
        """return traceless quadrupoles and octupoles"""
        dmaq, dmao = solvshift.efprot.tracls(self.__dmaq.copy(), self.__dmao.copy())
        return dmaq, dmao
    
#    def get_traceless_1(self, ravel=False):
    def get_traceless_1(self):
        """return traceless 1st derivatives wrt modes of quadrupoles and octupoles"""
        dmaq1 = self.__dmaq1.ravel()
        dmao1 = self.__dmao1.ravel()
        dmaq1, dmao1 = solvshift.efprot.tracl1(dmaq1.copy(), dmao1.copy(), self.__nmodes, self.__ndma)
        #if not ravel:
        #   self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes, self.__ndma, 6)
        #   self.__dmao1 = self.__dmao1.reshape(self.__nmodes, self.__ndma, 10)
        #print self.__dmao1.shape
        return dmaq1, dmao1
    
    def get_traceless_2(self, mode):
        """return traceless  2nd derivatives wrt modes of quadrupoles and octupoles"""
        mode_sder_index = list(self.__mode).index(mode-1)
        dmaq2, dmao2 = solvshift.efprot.tracls(self.__dmaq2[mode_sder_index].copy(), self.__dmao2[mode_sder_index].copy())
        return dmaq2, dmao2
   
    # UTILITY METHODS
 
    def write(self,file='slv.frg',par=None):
        """writes the parameters into a file"""
        f = open(file,'w')
        if par is not None:
           pass
        # basic molecular data
        if self.__name       is not None: self._write_preambule(f) 
        if self.__atoms      is not None: self._write_atoms(f)
        if self.__pos        is not None: self._write_pos(f)
        if self.__origin     is not None: self._write_origin(f)
        if self.__atno       is not None: self._write_atno(f)
        if self.__atms       is not None: self._write_atms(f)
        # frequency analysis
        if self.__redmass    is not None: self._write_redmass(f) 
        if self.__freq       is not None: self._write_freq(f)
        if self.__lvec       is not None: self._write_lvec(f)
        if self.__gijk       is not None: self._write_gijk(f)
        # population analysis
        if self.__esp        is not None: self._write_esp(f)   
        if self.__chlpg      is not None: self._write_chlpg(f)
        if self.__rdma       is not None: self._write_rdma(f)
        if self.__dmac       is not None: self._write_dmac(f)
        if self.__dmad       is not None: self._write_dmad(f)
        if self.__dmaq       is not None: self._write_dmaq(f)
        if self.__dmao       is not None: self._write_dmao(f)
        if self.__dmac1      is not None: self._write_dmac1(f)
        if self.__dmad1      is not None: self._write_dmad1(f)
        if self.__dmaq1      is not None: self._write_dmaq1(f)
        if self.__dmao1      is not None: self._write_dmao1(f)
        if self.__dmac2      is not None: self._write_dmac2(f)
        if self.__dmad2      is not None: self._write_dmad2(f)
        if self.__dmaq2      is not None: self._write_dmaq2(f)
        if self.__dmao2      is not None: self._write_dmao2(f)
        if self.__rpol       is not None: self._write_rpol(f)
        if self.__dpol       is not None: self._write_dpol(f)
        if self.__dpol1      is not None: self._write_dpol1(f)
        if self.__dpoli      is not None: self._write_dpoli(f)
        if self.__dpoli1     is not None: self._write_dpoli1(f)
        # EFP parameters
        if self.__lmoc       is not None: self._write_lmoc(f)  
        if self.__lmoc1      is not None: self._write_lmoc1(f)
        if self.__fock       is not None: self._write_fock(f)
        if self.__fock1      is not None: self._write_fock1(f)
        if self.__vecl       is not None: self._write_vecl(f)
        if self.__vecl1      is not None: self._write_vecl1(f)
        if self.__fckc       is not None: self._write_fckc(f)
        if self.__fckc1      is not None: self._write_fckc1(f)
        if self.__vecc       is not None: self._write_vecc(f)
        if self.__vecc1      is not None: self._write_vecc1(f)
        # DMS parameters
        if self.__opdm       is not None: self._write_opdm(f)
        if self.__caocc      is not None: self._write_caocc(f)
        if self.__cavir      is not None: self._write_cavir(f)
        if self.__dms_10     is not None: self._write_dms_10(f)
        if self.__dms_20     is not None: self._write_dms_20(f)
        if self.__dms_01     is not None: self._write_dms_01(f)
        if self.__dms_02     is not None: self._write_dms_02(f)
        if self.__dms_03     is not None: self._write_dms_03(f)
        if self.__dms_04     is not None: self._write_dms_04(f)
        if self.__dms_05     is not None: self._write_dms_05(f)
        if self.__dms_06     is not None: self._write_dms_06(f)
        if self.__dms_11     is not None: self._write_dms_11(f)
        if self.__dms_12     is not None: self._write_dms_12(f)
        if self.__dms_21     is not None: self._write_dms_21(f)
        if self.__dms_22     is not None: self._write_dms_22(f)
        #
        f.close()
        return
    
    def translate(self,transl):
        """translate tensors by <transl> cartesian displacement vector"""
        if self.__pos        is not None: self.__pos   += transl 
        if self.__lmoc       is not None: self.__lmoc  += transl
        if self.__rpol       is not None: self.__rpol  += transl
        if self.__rdma       is not None: self.__rdma  += transl
        return
    
    def rotate(self,rot):
        """rotate the tensors by <rot> unitary matrix"""
        # transform the atomic position static and dynamic information
        if self.__pos   is not None:
           self.__pos    = numpy.dot(self.__pos , rot)
        if self.__lmoc  is not None:
           self.__lmoc   = numpy.dot(self.__lmoc, rot)
        if self.__lmoc1 is not None:
           self.__lmoc1  = numpy.dot(self.__lmoc1,rot)
        if self.__lvec  is not None:
           self.__lvec   = numpy.dot(self.__lvec, rot)
        if self.__rpol  is not None:
           self.__rpol   = numpy.dot(self.__rpol, rot)
        if self.__rdma  is not None:
           self.__rdma   = numpy.dot(self.__rdma, rot)
        # transform dipoles, quadrupoles and octupoles!
        if self.__dmac  is not None:
           self.__dmad, self.__dmaq, self.__dmao = \
           solvshift.efprot.rotdma(self.__dmad, self.__dmaq, self.__dmao,rot)
        if self.__dmac1 is not None:
           self.__dmad1 = self.__dmad1.ravel()
           self.__dmaq1 = self.__dmaq1.ravel()
           self.__dmao1 = self.__dmao1.ravel()
           self.__dmad1, self.__dmaq1, self.__dmao1 = \
           solvshift.efprot.rotdm1(self.__dmad1, self.__dmaq1, self.__dmao1, rot,
                        self.__nmodes, self.__ndma)
           self.__dmad1 = self.__dmad1.reshape(self.__nmodes, self.__ndma, 3)
           self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes, self.__ndma, 6)
           self.__dmao1 = self.__dmao1.reshape(self.__nmodes, self.__ndma, 10)
        if self.__dmac2 is not None:
           for i in range(self.__nmodes_sder):
               self.__dmad2[i], self.__dmaq2[i], self.__dmao2[i] = \
               solvshift.efprot.rotdma(self.__dmad2[i], self.__dmaq2[i], self.__dmao2[i],rot)
        # transform distributed polarizabilities!
        if self.__dpol   is not None:
           for i in range(self.__npol):
               self.__dpol[i] = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpol[i], rot))
        if self.__dpol1  is not None:
           for i in range(self.__nmodes):
               for j in range(self.__npol):
                   self.__dpol1[i,j]  = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpol1[i,j], rot))
        if self.__dpoli  is not None:
           for i in range(self.__npol):
               for j in range(12):
                   self.__dpoli[i,j] = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpoli[i,j], rot))
        if self.__dpoli1 is not None:
           for i in range(self.__nmodes):
               for j in range(self.__npol):
                   for k in range(12):
                       self.__dpoli1[i,j,k]  = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpoli1[i,j,k], rot))
        # transform the wave function!
        if self.__vecl  is not None:
           bfs = self.get_bfs()
           typs= bfs.get_bfst().sum(axis=1)
           self.__vecl = solvshift.efprot.vecrot(self.__vecl, rot, typs)
           if self.__vecc  is not None:
              self.__vecc  = solvshift.efprot.vecrot(self.__vecc, rot, typs)
           if self.__vecl1 is not None:
              self.__vecl1 = solvshift.efprot.vc1rot(self.__vecl1, rot, typs)
           if self.__vecc1 is not None:
              self.__vecc1 = solvshift.efprot.vc1rot(self.__vecc1, rot, typs)
        # transform the DMS tensors
        if self.__opdm is not None:
           self.__bfs = self.get_bfs()
           m, R = rotate_ao_matrix(numpy.identity(self.__nbasis), rot, self.__bfs, return_rot=True, orbitals=True)
           self.__R = R; del m
           self.__opdm = self._rotate_dms_nn(self.__opdm, R)
           self.__caocc= self._rotate_aomo(self.__caocc, R)
           self.__cavir= self._rotate_aomo(self.__cavir, R)
        if self.__dms_10 is not None:
           self.__dms_10 = self._rotate_dms_nnN3 (self.__dms_10, R, rot)
        if self.__dms_20 is not None:                          
           self.__dms_20 = self._rotate_dms_nnN33(self.__dms_20, R, rot)
        if self.__dms_01 is not None:                          
           self.__dms_01 = self._rotate_dms_nn   (self.__dms_01, R     )
        if self.__dms_02 is not None:                          
           self.__dms_02 = self._rotate_dms_nn   (self.__dms_02, R     )
        if self.__dms_03 is not None:                          
           self.__dms_03 = self._rotate_dms_nn   (self.__dms_03, R     )
        if self.__dms_04 is not None:                          
           self.__dms_04 = self._rotate_dms_nn   (self.__dms_04, R     )
        if self.__dms_05 is not None:                          
           self.__dms_05 = self._rotate_dms_nn   (self.__dms_05, R     )
        if self.__dms_06 is not None:                          
           self.__dms_06 = self._rotate_dms_nn   (self.__dms_06, R     )

        if self.__dms_11 is not None:                          
           self.__dms_11 = self._rotate_dms_nnN3 (self.__dms_11, R, rot)
        if self.__dms_12 is not None:                          
           self.__dms_12 = self._rotate_dms_nnN3 (self.__dms_12, R, rot)
        if self.__dms_21 is not None:                          
           self.__dms_21 = self._rotate_dms_nnN33(self.__dms_21, R, rot)
        if self.__dms_22 is not None:                          
           self.__dms_22 = self._rotate_dms_nnN33(self.__dms_22, R, rot)

        return

    def _rotate_aomo(self, C, R):
        C_rot = numpy.einsum("bc,ba->ac", C, R)
        return C_rot
    def _rotate_dms_nn(self, B, R):
        B_rot = numpy.einsum("ab,ac,bd->cd", B, R, R)
        return B_rot
    def _rotate_dms_nnN3(self, B, R, r):
        B_rot = numpy.einsum("abiu,ac,bd,ux->cdix", B, R, R, r)
        return B_rot
    def _rotate_dms_nnN33(self, B, R, r):
        B_rot = numpy.einsum("abiuw,ac,bd,ux,wy->cdixy", B, R, R, r, r)
        return B_rot

    def superimpose(self, *args, **kwargs): return self.sup(*args, **kwargs) 

    def sup(self, xyz, suplist=None, dxyz=None, rotran=None):
        """superimpose structures <xyz> and <self.__pos>. Return rms in A.U."""
        # superimposition
        if rotran is None:
           if len(xyz)!=1:                                                             
              s = Superimposer()
              if dxyz is None:
                 if suplist is None: s.set(xyz,self.__pos)
                 else:               s.set(xyz[suplist],self.__pos[suplist])
                 s.run()                       
                 rms         = s.get_rms()
                 rot, transl = s.get_rotran()
              else:
                 #libbbg.utilities.PRINTL(self.__pos*self.BohrToAngstrom,'','')
                 if suplist is None: s.set( xyz-dxyz,self.__pos )
                 else:               s.set((xyz-dxyz)[suplist],(self.__pos)[suplist])
                 s.run()
                 rms         = s.get_rms()
                 rot, transl = s.get_rotran()
                 #pos = numpy.dot(self.__pos , rot) + transl
                 #s = libbbg.utilities.SVDSuperimposer()
                 #if suplist is None: s.set(pos ,self.__pos + dxyz )
                 #else:               s.set(pos[suplist],(self.__pos + dxyz)[suplist])
                 #s.run()
                 #rms         = s.get_rms()
                 #rot, transl = s.get_rotran()
                 #print "HUH", rot, transl, rms
           else: # the case of 1-atom "molecule" like Na+ cation
              rot    = numpy.identity(3)
              transl = xyz - self.__pos
              rms    = 0.0
        # user-defined rotation matrix and translation vector
        else:
           rot, transl = rotran
           rms = 0.0

        # save the actual transformation tensors
        self.__rot = rot
        self.__transl = transl
        self.__rms = rms
        # perform transformations
        self.rotate(rot)
        self.translate(transl)
        ##self.__pos       = s.get_transformed()
        #self.__pos = numpy.dot(self.__pos , rot) + transl
        #if self.__lmoc  is not None: self.__lmoc   = numpy.dot(self.__lmoc , rot) + transl
        #if self.__rpol  is not None: self.__rpol   = numpy.dot(self.__rpol , rot) + transl
        #if self.__rdma  is not None: self.__rdma   = numpy.dot(self.__rdma , rot) + transl
        #if self.__lmoc1 is not None: self.__lmoc1  = numpy.dot(self.__lmoc1, rot)
        #if self.__lvec  is not None: self.__lvec   = numpy.dot(self.__lvec , rot)
        ## transform dipoles, quadrupoles and octupoles!
        #if self.__dmac  is not None:
        #   self.__dmad, self.__dmaq, self.__dmao = \
        #   solvshift.efprot.rotdma(self.__dmad, self.__dmaq, self.__dmao, rot)
        #if self.__dmac1 is not None:
        #   self.__dmad1 = self.__dmad1.ravel()
        #   self.__dmaq1 = self.__dmaq1.ravel()
        #   self.__dmao1 = self.__dmao1.ravel()
        #   self.__dmad1, self.__dmaq1, self.__dmao1 = \
        #   solvshift.efprot.rotdm1(self.__dmad1, self.__dmaq1, self.__dmao1, rot,
        #                self.__nmodes, self.__ndma)
        #   self.__dmad1 = self.__dmad1.reshape(self.__nmodes, self.__ndma, 3)
        #   self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes, self.__ndma, 6)
        #   self.__dmao1 = self.__dmao1.reshape(self.__nmodes, self.__ndma, 10)
        #if self.__dmac2 is not None:
        #   for i in range(self.__nmodes_sder):
        #       self.__dmad2[i], self.__dmaq2[i], self.__dmao2[i] = \
        #       solvshift.efprot.rotdma(self.__dmad2[i], self.__dmaq2[i], self.__dmao2[i], rot)
        ## transform distributed polarizabilities!
        #if self.__dpol   is not None:
        #   for i in range(self.__npol):
        #       self.__dpol[i] = numpy.dot(numpy.transpose(rot),numpy.dot(self.__dpol[i],rot))
        #if self.__dpol1  is not None:
        #   for i in range(self.__nmodes):
        #       for j in range(self.__npol):
        #           self.__dpol1[i,j]  = numpy.dot(numpy.transpose(rot),numpy.dot(self.__dpol1[i,j],rot))
        #if self.__dpoli  is not None:
        #   for i in range(self.__npol):
        #       for j in range(12):
        #           self.__dpoli[i,j] = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpoli[i,j], rot))
        #if self.__dpoli1 is not None:
        #   for i in range(self.__nmodes):
        #       for j in range(self.__npol):
        #           for k in range(12):
        #               self.__dpoli1[i,j,k]  = numpy.dot(numpy.transpose(rot), numpy.dot(self.__dpoli1[i,j,k], rot))
        ## - wave function
        #if self.__vecl  is not None:
        #   bfs = self.get_bfs()
        #   typs= bfs.get_bfst().sum(axis=1)
        #   self.__vecl = solvshift.efprot.vecrot(self.__vecl, rot, typs)
        #   if self.__vecc  is not None:
        #      self.__vecc  = solvshift.efprot.vecrot(self.__vecc, rot, typs)
        #   if self.__vecl1 is not None:
        #      self.__vecl1 = solvshift.efprot.vc1rot(self.__vecl1, rot, typs)
        #   if self.__vecc1 is not None:
        #      self.__vecc1 = solvshift.efprot.vc1rot(self.__vecc1, rot, typs)
        return rms
    
    def copy(self):
        """return a deepcopy of me!"""
        return copy.deepcopy(self)
    
    def reorder(self,ord,dma=True): #TODO
        """Reorder the turn of atoms in the parameter objects according to <ord>.
<dma> option is designed for the case when dma is not atomic and mid-bond-based (then False).
In other words - dma=True if DMA, CAMM and CBAMM are used. False if LMTP and others ..."""
        LIST1    = self.get_bfs().get_bfsl() + 1
        sim      = numpy.zeros((self.__natoms, 2),int)
        sim[:,0] = numpy.arange(1, self.__natoms+1, 1)
        sim[:,1] = ord
        # reorder atomic positions
        self.__pos = libbbg.utilities.reorder(self.__pos,  sim)
        self.__atno= libbbg.utilities.reorder(self.__atno, sim)
        self.__atms= libbbg.utilities.reorder(self.__atms, sim)
        # reorder multipole moments
        if dma:
           if self.__rdma is not None:
              self.__rdma[:self.__natoms] = libbbg.utilities.reorder(self.__rdma[:self.__natoms], sim)
           if self.__dmac is not None:
              self.__dmac[:self.__natoms] = libbbg.utilities.reorder(self.__dmac[:self.__natoms], sim)
              self.__dmad[:self.__natoms] = libbbg.utilities.reorder(self.__dmad[:self.__natoms], sim)
              self.__dmaq[:self.__natoms] = libbbg.utilities.reorder(self.__dmaq[:self.__natoms], sim)
              self.__dmao[:self.__natoms] = libbbg.utilities.reorder(self.__dmao[:self.__natoms], sim)
           #
           if self.__dmac1 is not None:
              self.__dmac1[:self.__natoms] = libbbg.utilities.reorder(self.__dmac1[:self.__natoms], sim, axis=1)
              self.__dmad1[:self.__natoms] = libbbg.utilities.reorder(self.__dmad1[:self.__natoms], sim, axis=1)
              self.__dmaq1[:self.__natoms] = libbbg.utilities.reorder(self.__dmaq1[:self.__natoms], sim, axis=1)
              self.__dmao1[:self.__natoms] = libbbg.utilities.reorder(self.__dmao1[:self.__natoms], sim, axis=1)
              #
              self.__dmac2[:self.__natoms] = libbbg.utilities.reorder(self.__dmac2[:self.__natoms], sim, axis=1)
              self.__dmad2[:self.__natoms] = libbbg.utilities.reorder(self.__dmad2[:self.__natoms], sim, axis=1)
              self.__dmaq2[:self.__natoms] = libbbg.utilities.reorder(self.__dmaq2[:self.__natoms], sim, axis=1)
              self.__dmao2[:self.__natoms] = libbbg.utilities.reorder(self.__dmao2[:self.__natoms], sim, axis=1)
        # reorder eigenvectors
        if self.__lvec is not None: 
           self.__lvec = libbbg.utilities.reorder(self.__lvec, sim, axis=1)
        # reorder the wavefunction
        if self.__vecl is not None:
           self._reorder_wfn(sim, LIST1)
        return
    
    def xyz(self, units='Angstrom'):
        """generate the xyz file contents (without 2 first lines)"""
        log = ''
        r = self.__pos.copy()
        Z = self.__atno
        if units.lower().startswith('angs'): r *= self.BohrToAngstrom
        for i,x in enumerate(r):
            log += '%3s' % self.__atoms[i]
            log += 3*'%16.6f' % tuple(x)
            log += '\n'
        return log 
   
    # --- protected
    
    def _create(self):
        """creates the namespace of parameter variables"""
        self.__name ,self.__natoms, self.__nbasis = None, None, None
        self.__nmos ,self.__nmodes, self.__basis  = None, None, None
        self.__naocc,self.__nbocc                 = None, None
        self.__navir,self.__nbvir                 = None, None
        self.__atoms,self.__shortname             = None, None
        self.__pos  ,self.__origin, self.__nsites = None, None, None
        self.__atno ,self.__atms  , self.__mode   = None, None, None
        self.__charge, self.__multiplicity        = None, None
        #
        self.__redmass, self.__freq, self.__lvec  = None, None, None
        self.__gijk                               = None
        #
        self.__fock ,self.__lmoc ,self.__vecl     = None, None, None 
        self.__fock1,self.__lmoc1,self.__vecl1    = None, None, None
        self.__ncmos,self.__vecc ,self.__vecc1    = None, None, None
        self.__fckc ,self.__fckc1                 = None, None
        #                                                            
        self.__esp  ,self.__chlpg,self.__npol     = None, None, None
        self.__dpol ,self.__dpol1,self.__rpol     = None, None, None
        self.__dpoli,self.__dpoli1                = None, None
        #                                                            
        self.__rdma, self.__dmac, self.__dmad     = None, None, None
        self.__dmaq, self.__dmao, self.__ndma     = None, None, None
        self.__dmac1,self.__dmad1,self.__dmaq1    = None, None, None
        self.__dmao1                              = None
        self.__dmac2,self.__dmad2,self.__dmaq2    = None, None, None
        self.__dmao2                              = None
        #
        self.__opdm,self.__caocc,self.__cavir     = None, None, None
        self.__dms_10,self.__dms_20,self.__dms_02 = None, None, None
        self.__dms_12,self.__dms_22               = None, None
        self.__dms_01,self.__dms_11,self.__dms_21 = None, None, None
        self.__dms_03,self.__dms_04               = None, None
        self.__dms_05,self.__dms_06               = None, None
        #
        mol_names = ('name' ,'basis' ,'method','natoms'   ,'nbasis',
                     'nmos' ,'nmodes','atoms' ,'shortname','nsites',
                     'ncmos','npol'  ,'ndma'  ,'mode'     ,'charge',
                     'naocc','nbocc' ,'navir' ,'nbvir'    ,'multiplicity',)
        sec_names = {  'mol': '[ molecule ]'                                        ,
                     'atoms': '[ Atomic symbols ]'                                  ,
                      'atno': '[ Atomic numbers ]'                                  ,
                      'atms': '[ Atomic masses ]'                                   ,
                       'pos': '[ Atomic coordinates ]'                              ,
                      'freq': '[ Harmonic frequencies ]'                            ,
                   'redmass': '[ Reduced masses ]'                                  ,
                      'lvec': '[ Mass-weighted eigenvectors ]'                      ,
                      'gijk': '[ Cubic anharmonic constants ]'                      ,
                    'origin': '[ DMTP origins ]'                                    ,
                       'esp': '[ ESP charges ]'                                     ,
                     'chlpg': '[ ChelpG charges ]'                                  ,
                      'rdma': '[ DMTP centers ]'                                    ,
                      'dmac': '[ DMTP charges ]'                                    ,
                      'dmad': '[ DMTP dipoles ]'                                    ,
                      'dmaq': '[ DMTP quadrupoles ]'                                ,
                      'dmao': '[ DMTP octupoles ]'                                  ,
                     'dmac1': '[ DMTP charges - first derivatives ]'                ,
                     'dmad1': '[ DMTP dipoles - first derivatives ]'                ,
                     'dmaq1': '[ DMTP quadrupoles - first derivatives ]'            ,
                     'dmao1': '[ DMTP octupoles - first derivatives ]'              ,
                     'dmac2': '[ DMTP charges - second derivatives wrt mode]'       ,
                     'dmad2': '[ DMTP dipoles - second derivatives wrt mode]'       ,
                     'dmaq2': '[ DMTP quadrupoles - second derivatives wrt mode]'   ,
                     'dmao2': '[ DMTP octupoles - second derivatives wrt mode]'     ,
                      'rpol': '[ Polarizable centers ]'                             ,
                      'dpol': '[ Distributed polarizabilities ]'                    ,
                     'dpol1': '[ Distributed polarizabilities - first derivatives ]',
                     'dpoli': '[ Distr. pol. wrt imaginary freq ]'                  ,
                    'dpoli1': '[ Distr. pol. wrt imaginary freq - 1st derivatives ]',
                     'lmoc' : '[ LMO centroids ]'                                   ,
                     'lmoc1': '[ LMO centroids - first derivatives ]'               ,
                     'fock' : '[ Fock matrix ]'                                     ,
                     'fock1': '[ Fock matrix - first derivatives ]'                 ,
                     'fckc' : '[ Canonical Fock matrix ]'                           ,
                     'fckc1': '[ Canonical Fock matrix - first derivatives ]'       ,
                     'vecl' : '[ AO->LMO matrix ]'                                  ,
                     'vecl1': '[ AO->LMO matrix - first derivatives ]'              ,
                     'vecc' : '[ AO->CMO matrix ]'                                  ,
                     'vecc1': '[ AO->CMO matrix - first derivatives ]'              ,
                     'opdm' : '[ One-particle density matrix - Alpha ]'             ,
                    'caocc' : '[ SCF AO-MO matrix - Occ, Alpha ]'                   ,
                    'cavir' : '[ SCF AO-MO matrix - Vir, Alpha ]'                   ,
                    'dms_10': '[ DMS 10 ]'                                          ,
                    'dms_20': '[ DMS 20 ]'                                          ,
                    'dms_01': '[ DMS 01 ]'                                          ,
                    'dms_02': '[ DMS 02 ]'                                          ,
                    'dms_03': '[ DMS 03 ]'                                          ,
                    'dms_04': '[ DMS 04 ]'                                          ,
                    'dms_05': '[ DMS 05 ]'                                          ,
                    'dms_06': '[ DMS 06 ]'                                          ,
                    'dms_11': '[ DMS 11 ]'                                          ,
                    'dms_12': '[ DMS 12 ]'                                          ,
                    'dms_21': '[ DMS 21 ]'                                          ,
                    'dms_22': '[ DMS 22 ]'                                          ,}
        self.__mol_names = mol_names
        self.__sec_names = sec_names
       #self.__atomic_symbols = {1: 'H', 2: 'He', 3: 'Li', 6:'C', 7:'N', 8:'O', 9:'F', 11:'Na', 12:'Mg', 16:'S', 17:'Cl'}
        return

    def _parse_dict(self,par):
        """save the memorials from par dictionary"""
        for key, val in par.items():
            if   key == 'fock' : self.__fock  = val
            elif key == 'fock1': self.__fock1 = val
            elif key == 'lmoc' : self.__lmoc  = val
            elif key == 'lmoc1': self.__lmoc1 = val
            elif key == 'vecl' : self.__vecl  = val
            elif key == 'vecl1': self.__vecl1 = val
            elif key == 'vecc' : self.__vecc  = val
            elif key == 'vecc1': self.__vecc1 = val
            elif key == 'fckc' : self.__fckc  = val
            elif key == 'fckc1': self.__fckc1 = val
            elif key == 'opdm' : self.__opdm  = val
            elif key == 'caocc': self.__caocc = val
            elif key == 'cavir': self.__cavir = val
            elif key =='dms_10': self.__dms_10= val
            elif key =='dms_20': self.__dms_20= val
            elif key =='dms_01': self.__dms_01= val
            elif key =='dms_02': self.__dms_02= val
            elif key =='dms_03': self.__dms_03= val
            elif key =='dms_04': self.__dms_04= val
            elif key =='dms_05': self.__dms_05= val
            elif key =='dms_06': self.__dms_06= val
            elif key =='dms_11': self.__dms_11= val
            elif key =='dms_12': self.__dms_12= val
            elif key =='dms_21': self.__dms_21= val
            elif key =='dms_22': self.__dms_22= val
        return

    def _make_dict(self):
        """returns dictionary with parameters"""
        par = {}
        # basic molecular data
        if self.__pos    is not None: par['pos'   ] = self.__pos
        if self.__origin is not None: par['origin'] = self.__origin
        if self.__atno   is not None: par['atno'  ] = self.__atno
        if self.__atms   is not None: par['atms'  ] = self.__atms
        if self.__ndma   is not None: par['ndma'  ] = self.__ndma
        if self.__npol   is not None: par['npol'  ] = self.__npol
        if self.__nmodes is not None: par['nmodes'] = self.__nmodes
        if self.__natoms is not None: par['natoms'] = self.__natoms
        if self.__name   is not None: par['name'  ] = self.__name
        if self.__basis  is not None: par['basis' ] = self.__basis
        if self.__nbasis is not None: par['nbasis'] = self.__nbasis
        if self.__nmos   is not None: par['nmos'  ] = self.__nmos
        if self.__naocc  is not None: par['naocc' ] = self.__naocc
        if self.__navir  is not None: par['navir' ] = self.__navir
        if self.__nbocc  is not None: par['nbocc' ] = self.__nbocc
        if self.__nbvir  is not None: par['nbvir' ] = self.__nbvir
        if self.__mode   is not None: par['mode'  ] = self.__mode
        if self.__charge is not None: par['charge'] = self.__charge
        if self.__multiplicity is not None: par['multiplicity'] = self.__multiplicity
        # frequency analysis
        if self.__redmass  is not None: par['redmass'] = self.__redmass
        if self.__freq     is not None: par['freq'   ] = self.__freq
        if self.__lvec     is not None: par['lvec'   ] = self.__lvec
        if self.__gijk     is not None: par['gijk'   ] = self.__gijk
        # EFP parameters
        if self.__lmoc  is not None: par['lmoc' ] = self.__lmoc
        if self.__lmoc1 is not None: par['lmoc1'] = self.__lmoc1
        if self.__fock  is not None: par['fock' ] = self.__fock
        if self.__fock1 is not None: par['fock1'] = self.__fock1
        if self.__vecl  is not None: par['vecl' ] = self.__vecl
        if self.__vecl1 is not None: par['vecl1'] = self.__vecl1
        if self.__vecc  is not None: par['vecc' ] = self.__vecc
        if self.__vecc1 is not None: par['vecc1'] = self.__vecc1
        if self.__fckc  is not None: par['fckc' ] = self.__fckc
        if self.__fckc1 is not None: par['fckc1'] = self.__fckc1
        # Population analysis
        if self.__esp   is not None: par['esp'  ] = self.__esp
        if self.__chlpg is not None: par['chlpg'] = self.__chlpg
        if self.__rdma  is not None: par['rdma' ] = self.__rdma
        if self.__dmac  is not None: par['dmac' ] = self.__dmac
        if self.__dmad  is not None: par['dmad' ] = self.__dmad
        if self.__dmaq  is not None: par['dmaq' ] = self.__dmaq
        if self.__dmao  is not None: par['dmao' ] = self.__dmao
        if self.__dmac1 is not None: par['dmac1'] = self.__dmac1
        if self.__dmad1 is not None: par['dmad1'] = self.__dmad1
        if self.__dmaq1 is not None: par['dmaq1'] = self.__dmaq1
        if self.__dmao1 is not None: par['dmao1'] = self.__dmao1
        if self.__dmac2 is not None: par['dmac2'] = self.__dmac2
        if self.__dmad2 is not None: par['dmad2'] = self.__dmad2
        if self.__dmaq2 is not None: par['dmaq2'] = self.__dmaq2
        if self.__dmao2 is not None: par['dmao2'] = self.__dmao2
        if self.__rpol  is not None: par['rpol' ] = self.__rpol
        if self.__dpol  is not None: par['dpol' ] = self.__dpol
        if self.__dpol1 is not None: par['dpol1'] = self.__dpol1
        if self.__dpoli  is not None: par['dpoli' ] = self.__dpoli
        if self.__dpoli1 is not None: par['dpoli1'] = self.__dpoli1
        # DMS tensors
        if self.__opdm   is not None: par['opdm']   = self.__opdm
        if self.__caocc  is not None: par['caocc']  = self.__caocc
        if self.__cavir  is not None: par['cavir']  = self.__cavir
        if self.__dms_10 is not None: par['dms_10'] = self.__dms_10
        if self.__dms_20 is not None: par['dms_20'] = self.__dms_20
        if self.__dms_01 is not None: par['dms_01'] = self.__dms_01
        if self.__dms_02 is not None: par['dms_02'] = self.__dms_02
        if self.__dms_03 is not None: par['dms_03'] = self.__dms_03
        if self.__dms_04 is not None: par['dms_04'] = self.__dms_04
        if self.__dms_05 is not None: par['dms_05'] = self.__dms_05
        if self.__dms_06 is not None: par['dms_06'] = self.__dms_06
        if self.__dms_11 is not None: par['dms_11'] = self.__dms_11
        if self.__dms_12 is not None: par['dms_12'] = self.__dms_12
        if self.__dms_21 is not None: par['dms_21'] = self.__dms_21
        if self.__dms_22 is not None: par['dms_22'] = self.__dms_22
        #
        return par
        
    def _reorder_wfn(self, sim, LIST1):
        """reorders the turn of basis functions according to <ord> list of atomic interchange"""
        P1, P2 = [], []
        #
        T1 = numpy.transpose(self.__vecl,axes=(1,0))
        #
        b=[]; n=1; curr= LIST1[0]
        for i in range(len(LIST1)-1):
            if curr==LIST1[i+1]: n+=1
            else:
                b.append(n)
                n=1
            curr=LIST1[i+1]
        b.append(n)
        #
        c = [ sum(b[:x]) for x in range(len(b)) ]
        c.append(len(LIST1))
        #
        for i in range(self.__natoms):
            P1.append( T1[c[i]:c[i+1]] )
        P1 = numpy.array(P1)
        #
        N1 = numpy.zeros(self.__natoms, dtype=object)
        #
        for i,j in sim:
            N1[i-1] = P1[j-1]
        #
        self.__vecl = numpy.concatenate(tuple(N1))
        self.__vecl = numpy.transpose(self.__vecl, axes=(1,0))
        ### DERIVATIVES
        if self.__vecl1 is not None:
           T2 = numpy.transpose(self.__vecl1, axes=(2,1,0))
           #
           for i in range(self.__natoms):
               P2.append( T2[c[i]:c[i+1]] )
           P2 = numpy.array(P2)
           #
           N2 = numpy.zeros(self.__natoms, dtype=object)
           for i,j in sim:
               N2[i-1] = P2[j-1]
           #
           self.__vecl1 = numpy.concatenate(tuple(N2))
           self.__vecl1 = numpy.transpose(self.__vecl1, axes=(2,1,0))
        return
    
    def _tr_lvec(self, lvec, nmodes, natoms):
        """Transpose axes and reshape the L-matrix from FREQ instance"""
        return lvec.transpose().reshape(nmodes,natoms,3)

    def _gauss_legendre_12(self, f, v=0.3000):
        """Gauss-Legendre quadrature in 12 points (numerical) from 0 to infinity:

\int_0^{\infty} f dq = \sum_{i=1}^{12} w_i * 2v/(1-t_i)**2 * f_i

It uses the substitution of variables:

              q = v(1+t)/(1-t)   dq = 2v dt / (1-t)**2

where v=0.3. Argument f is an array with the first dimension equal to 12
"""
        g_12 = """\
        11      0.0471753363865118     -0.9815606342467192
        9       0.1069393259953184     -0.9041172563704749
        7       0.1600783285433462     -0.7699026741943047
        5       0.2031674267230659     -0.5873179542866175
        3       0.2334925365383548     -0.3678314989981802
        1       0.2491470458134028     -0.1252334085114689
        2       0.2491470458134028      0.1252334085114689
        4       0.2334925365383548      0.3678314989981802
        6       0.2031674267230659      0.5873179542866175
        8       0.1600783285433462      0.7699026741943047
        10      0.1069393259953184      0.9041172563704749
        12      0.0471753363865118      0.9815606342467192
        """
        g_12 = numpy.array(g_12.split(), numpy.float64).reshape(12,3)
        w_12 = g_12[:,1]
        t_12 = g_12[:,2]
        integral = numpy.zeros(list(f.shape)[1:], numpy.float64)
        for i in range(12): integral += w_12[i] * f[i] / (1.0 - t_12[i])**2
        return 2.0 * v * integral

    # --------------------------------------------------------- #
    #            P R O P E R T Y   E V A L U A T O R S          #
    # --------------------------------------------------------- #

    def _get_property_solcamm(self, mode=None, ea=True, dma=True):
        """
Compute SolCAMM parameters. If dma=False it returns a tuple:

rdma, pos, solcamm_c, solcamm_d, solcamm_q, solcamm_o

If dma=True it returns Coulomb.py DMA object
"""
        assert mode is not None, " You MUST specify <mode> to compute solvatochromic parameters (in normal numbers, helico)"
        if ea: assert (mode in self.__mode+1), " The second derivatives wrt mode %d are not evaluated for %s!" % (mode, self.__name)

        par    = self.get()
        pos    = par['pos']
        rdma   = par['rdma']
        redmass= par['redmass']
        freq   = par['freq']
        gijj   = par['gijk'][:,mode-1,mode-1]
        dmac1  = par['dmac1']; dmad1  = par['dmad1']
        dmaq1  = par['dmaq1']; dmao1  = par['dmao1']
        if ea:
           dmac2  = par['dmac2']; dmad2  = par['dmad2']
           dmaq2  = par['dmaq2']; dmao2  = par['dmao2']

        # vibrational weights
        g_wi= gijj/(redmass*freq**2)
        g_wj= 2.00 * redmass[mode-1] * freq[mode-1]

        # mechanical anharmonicity approximation
        solcamm_c = -numpy.tensordot( g_wi, dmac1, (0,0))
        solcamm_d = -numpy.tensordot( g_wi, dmad1, (0,0))
        solcamm_q = -numpy.tensordot( g_wi, dmaq1, (0,0))
        solcamm_o = -numpy.tensordot( g_wi, dmao1, (0,0))

        # electronic anharmonicity
        if ea:
           m = list(self.__mode).index(mode-1)
           solcamm_c += dmac2[m]; solcamm_d += dmad2[m]
           solcamm_q += dmaq2[m]; solcamm_o += dmao2[m]
        solcamm_c /= g_wj; solcamm_d /= g_wj
        solcamm_q /= g_wj; solcamm_o /= g_wj

        if dma:
           name = '%s - SolCAMM parameters for mode %d' % (self.__name, mode)
           if ea: name+= ' (MEA+EA approximation)'
           else : name+= ' (MEA approximation)'
           solcamm = libbbg.dma.DMA(name=name, pos=pos, origin=rdma, traceless=False, atoms=None,
                                    q=solcamm_c, m=solcamm_d, T=solcamm_q, O=solcamm_o, H=None)
           return solcamm
        else:  return rdma, pos, solcamm_c, solcamm_d, solcamm_q, solcamm_o

    def _get_property_c6d(self, frag):
        """Compute LMO distributed isotropic (scalar) C6 coefficients between self and frag"""
        par_self   = self.get()        ; par_frag   = frag.get()
        dpoli_self = par_self['dpoli'] ; dpoli_frag = par_frag['dpoli']
        nmos_self  = par_self['nmos' ] ; nmos_frag  = par_frag['nmos']
        c6d = numpy.zeros((nmos_self, nmos_frag), numpy.float64)
        for i in range(nmos_self):
            a_i = dpoli_self[i]
            a_oi= a_i.trace(axis1=1,axis2=2)/3.000
            for j in range(nmos_frag):
                a_j = dpoli_frag[j]
                a_oj= a_j.trace(axis1=1,axis2=2)/3.000
                c   = self._gauss_legendre_12(a_oi * a_oj) * 3.00 / math.pi
                c6d[i,j] = c
        return c6d
 
    def _get_property_c6(self, frag):
        """Compute LMO total isotropic (scalar) C6 coefficient between self and frag"""
        c6 = self._get_property_c6d(frag).sum()
        return c6

    def _get_property_edisp6(self, frag):
        """Computes dispersion energy between self and frag"""
        lmoc_self = self.get()['lmoc']
        lmoc_frag = frag.get()['lmoc']
        c6d = self._get_property_c6d(frag)
        E = 0
        nmosa, nmosb = c6d.shape
        R = lmoc_self[:,:,numpy.newaxis] - lmoc_frag.T[numpy.newaxis,:,:]
        R = R*R; R = numpy.sum(R, axis=1)
        R = numpy.sqrt(R)
        R6= R**6
        E = - numpy.sum( c6d / R6, axis=None)
        #for i in range(nmosa):
        #    for j in range(nmodb):
        #        r6  = (lmoc_self[i] - lmoc_frag[j])**6
        #        E  -= c6d[i,j] / r6
        return E

    def _get_property_solc6d(self, frag):  # in the future add mode=... keyword
        """Compute solvatochromic LMO distributed isotropic (scalar) C6 coefficients between self and frag"""
        par_self   = self.get()        ; par_frag   = frag.get()
        dpoli_self = par_self['dpoli1']; dpoli_frag = par_frag['dpoli']
        nmos_self  = par_self['nmos' ] ; nmos_frag  = par_frag['nmos']
        nmodes     = par_self['nmodes']
        redmass= par_self['redmass']
        freq   = par_self['freq']
        gijk   = par_self['gijk']
        c6d = numpy.zeros((nmodes, nmos_self, nmos_frag), numpy.float64)
        for M in range(nmodes):
            vib_m = 1./(2.0 * redmass[M] * freq[M])
            for i in range(nmos_self):                                              
                 for j in range(nmos_frag):
                    c = 0.0
                    a_j = dpoli_frag[j]
                    a_oj= a_j.trace(axis1=1,axis2=2)/3.000
                    for k in range(nmodes):
                        a_i = dpoli_self[k,i]
                        a_oi= a_i.trace(axis1=1,axis2=2)/3.000
                        vib_w = gijk[k,M,M] / (redmass[k] * freq[k]**2)
                        c  -= vib_w * self._gauss_legendre_12(a_oi * a_oj) * 3.00 / math.pi
                    c6d[M,i,j] = c * vib_m
        return c6d

    def _get_property_solc6(self, frag):  # in the future add mode=... keyword
        """Compute solvatochromic total isotropic (scalar) C6 coefficient between self and frag"""
        solc6 = self._get_property_solc6d(frag).sum(axis=1).sum(axis=1)
        return solc6
   
    # --------------------------------------------------------- #
    #            R E A D I N G    P R O C E D U R E S           #
    # --------------------------------------------------------- #
    
    def _read(self, section):
        """read the appropriate field from section and save it in SLVPAR instance"""
        section = section.split('\n')
        querry  = section[0]
        for key, val in self.__sec_names.items():
            if val[2:] in querry: break

        ### basic molecular specification
        if key == 'mol':
            data = section[1:]
            for line in data:
                if '=' in line:
                    name, arg = line.split('=')
                    name = name.strip(); arg = arg.strip()
                    if name == 'name':
                        self.__name = arg
                    if name == 'shortname':
                        self.__shortname = arg
                    if name == 'atoms':
                        self.__atoms = arg.split(',')
                    if name == 'charge':
                        self.__charge= int(arg)
                    if name.startswith('multiplic'):
                        self.__multiplicity= int(arg)
                    if name == 'natoms':
                        self.__natoms = int(arg)
                    if name == 'nsites':
                        self.__nsites = int(arg)
                    if name == 'basis':
                        self.__method, self.__basis = arg.split('/')
                    if name == 'nbasis':
                        self.__nbasis = int(arg)
                    if name == 'nmos':
                        self.__nmos = int(arg)
                    if name == 'naocc':
                        self.__naocc = int(arg)
                    if name == 'navir':
                        self.__navir = int(arg)
                    if name == 'nbocc':
                        self.__nbocc = int(arg)
                    if name == 'nbvir':
                        self.__nbvir = int(arg)
                    if name == 'ncmos':
                        self.__ncmos = int(arg)
                    if name == 'nmodes':
                        self.__nmodes = int(arg)
                    if name == 'npol':
                        self.__npol = int(arg)
                    if name == 'ndma':
                        self.__ndma = int(arg)
                    if name == 'mode':
                        self.__mode = libbbg.utilities.text_to_list(arg, delimiter=',') - 1
                        self.__nmodes_sder = len(self.__mode)

        ### more advanced information
        else:
            data = [] ; N = int(querry.split('N=')[1])
            n = self._nlines(N)
            for i in range(n): data += section[i+1].split()
            dtype= str if 'Atomic symbols' in querry else numpy.float64 
            data = numpy.array(data, dtype=dtype)
            # ------------------------------------ MOL --------------------------------------------
            # Atomic symbols
            if   key == 'atoms':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic symbols ]!'
                 assert self.__natoms == N, merror
                 self.__atoms = data
            # Atomic coordinates
            elif key == 'pos':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic coordinates ]!'
                 assert self.__natoms == N/3, merror
                 data = data.reshape(self.__natoms, 3)
                 self.__pos = data
            # DMTP origins
            elif key == 'origin':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP origins ]!'
                 assert 1==1, merror
                 data = data.reshape(N/3, 3)
                 self.__origin = data
            # Atomic numbers
            elif key == 'atno':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic numbers ]!'
                 assert self.__natoms == N, merror
                 self.__atno = numpy.array(data, int)
            # Atomic masses
            elif key == 'atms':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic masses ]!'
                 assert self.__natoms == N, merror
                 self.__atms = data
            # ------------------------------------ DMTP -------------------------------------------
            # ESP charges
            elif key == 'esp':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ ESP charges ]!'
                 assert self.__natoms == N, merror
                 self.__esp = data
            # ChelpG charges
            elif key == 'chlpg':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ ChelpG charges ]!'
                 assert self.__natoms == N, merror
                 self.__chlpg = data
            # DMA/CAMM
            # distributed center coordinates
            elif key == 'rdma':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP centers ]!'
                 assert self.__ndma == N/3, merror
                 data = data.reshape(self.__ndma, 3)
                 self.__rdma = data
            # distributed charges
            elif key == 'dmac':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma == N, merror
                 self.__dmac = data
            # distributed dipoles
            elif key == 'dmad':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma == N/3, merror
                 data = data.reshape(self.__ndma, 3)
                 self.__dmad = data
            # distributed quadrupoles
            elif key == 'dmaq':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma == N/6, merror
                 data = data.reshape(self.__ndma, 6)
                 self.__dmaq = data
            # distributed octupoles
            elif key == 'dmao':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma == N/10, merror
                 data = data.reshape(self.__ndma, 10)
                 self.__dmao = data
            # - DMTP derivatives -
            # first derivatives
            # distributed charge 1st derivatives
            elif key == 'dmac1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma*self.__nmodes == N, merror
                 self.__dmac1 = data.reshape(self.__nmodes, self.__ndma)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis]
                 self.__dmac1 = temp * self.__dmac1
            # distributed dipole 1st derivatives
            elif key == 'dmad1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma*self.__nmodes == N/3, merror
                 data = data.reshape(self.__nmodes, self.__ndma, 3)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                 self.__dmad1 = temp * data
            # distributed quadrupoles 1st derivatives
            elif key == 'dmaq1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma*self.__nmodes == N/6, merror
                 data = data.reshape(self.__nmodes, self.__ndma, 6)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                 self.__dmaq1 = temp * data
            # distributed octupoles 1st derivatives
            elif key == 'dmao1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma*self.__nmodes == N/10, merror
                 data = data.reshape(self.__nmodes, self.__ndma, 10)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                 self.__dmao1 = temp * data
            # second derivatives
            # distributed charges 2nd derivatives
            elif key == 'dmac2':
                 merror = 'ndma and mode in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma == N/self.__nmodes_sder, merror
                 data = data.reshape(self.__nmodes_sder, self.__ndma)
                 temp = self.__redmass[self.__mode][:, numpy.newaxis]
                 self.__dmac2 = temp * data
            # distributed dipoles 2nd derivatives
            elif key == 'dmad2':
                 merror = 'ndma and mode in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma == N/(3*self.__nmodes_sder), merror
                 data = data.reshape(self.__nmodes_sder, self.__ndma, 3)
                 temp = self.__redmass[self.__mode][:, numpy.newaxis, numpy.newaxis]
                 self.__dmad2 = temp * data
            # distributed quadrupoles 2nd derivatives
            elif key == 'dmaq2':
                 merror = 'ndma and mode in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma == N/(6*self.__nmodes_sder), merror
                 data = data.reshape(self.__nmodes_sder, self.__ndma, 6)
                 temp = self.__redmass[self.__mode][:, numpy.newaxis, numpy.newaxis]
                 self.__dmaq2 = temp * data
            # distributed octupoles 2nd derivatives
            elif key == 'dmao2':
                 merror = 'ndma and mode in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma == N/(10*self.__nmodes_sder), merror
                 data = data.reshape(self.__nmodes_sder, self.__ndma, 10)
                 temp = self.__redmass[self.__mode][:, numpy.newaxis, numpy.newaxis]
                 self.__dmao2 = temp * data
            # CABMM
            elif key == 'cabmm':
                 pass
            # LMTP
            elif key == 'lmtp':
                 pass
            # ------------------------------------ DPOL -------------------------------------------
            elif key == 'rpol':
                 merror = None
                 self.__rpol = data.reshape(self.__npol, 3)
            elif key == 'dpol':
                 merror = 'npol in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distributed polarizabilities ]!'
                 assert self.__npol == (N/9), merror
                 self.__dpol = data.reshape(self.__npol, 3, 3)
            elif key == 'dpol1':
                 merror = 'npol and nmodes in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distributed polarizabilities ]!'
                 assert self.__npol == (N/(9*self.__nmodes)), merror
                 self.__dpol1 = data.reshape(self.__nmodes, self.__npol, 3, 3)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis, numpy.newaxis]
                 self.__dpol1 = temp * self.__dpol1
            elif key == 'dpoli':
                 merror = 'npol in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distr. pol. wrt imaginary freq ]'
                 assert self.__npol == (N/(9*12)), merror
                 self.__dpoli = data.reshape(self.__npol, 12, 3, 3)
            elif key == 'dpoli1':
                 merror = 'npol and nmodes in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distr. pol. wrt imaginary freq - 1st derivatives ]'
                 assert self.__npol == (N/(9*12*self.__nmodes)), merror
                 self.__dpoli1 = data.reshape(self.__nmodes, self.__npol, 12, 3, 3)
                 temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis, numpy.newaxis, numpy.newaxis]
                 self.__dpoli1 = temp * self.__dpoli1
            # ------------------------------------ FREQ -------------------------------------------
            # Harmonic frequencies
            elif key == 'freq':
                 merror = 'nmodes in section [ molecule ] '
                 merror+= 'is not consistent with section [ Harmonic frequencies ]!'
                 assert self.__nmodes == N, merror
                 self.__freq = data
            # Reduced masses
            elif  key == 'redmass':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Reduced masses ]!'
                  assert self.__nmodes == N, merror
                  self.__redmass = data
            # Mass-weighted eigenvectors
            elif  key == 'lvec':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Mass-weighted eigenvectors ]!'
                  assert self.__nmodes == N/(self.__natoms*3), merror
                  self.__lvec = data.reshape(self.__nmodes, self.__natoms, 3)
            # Cubic anharmonic constants
            elif  key == 'gijk':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Cubic anharmonic constants ]!'
                  assert 1==1, merror
                  gijk = numpy.zeros((self.__nmodes, self.__nmodes, self.__nmodes), dtype=numpy.float64)
                  K = 0
                  for i in range(self.__nmodes):
                      for j in range(i+1):
                          for k in range(j+1):
                              d = data[K]
                              gijk[i,j,k] = d
                              gijk[i,k,j] = d
                              gijk[j,i,k] = d
                              gijk[j,k,i] = d
                              gijk[k,i,j] = d
                              gijk[k,j,i] = d
                              K += 1                
                  self.__gijk = gijk
            # ------------------------------------ WFN --------------------------------------------
            # LMO centroids
            elif  key == 'lmoc':
                  merror = 'nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ LMO centroids ]!'
                  assert self.__nmos == N/3, merror
                  data = data.reshape(self.__nmos, 3)
                  self.__lmoc = data
            # LMO centroids first derivatives
            elif  key == 'lmoc1':
                  merror = 'nmos and nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ LMO centroids - first derivatives ]!'
                  assert self.__nmos == N/(3*self.__nmodes), merror
                  data = data.reshape(self.__nmodes, self.__nmos, 3)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                  data = temp * data
                  self.__lmoc1 = data
            # Fock matrix
            elif  key == 'fock':
                  merror = 'nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Fock matrix ]!'
                  assert self.__nmos == (math.sqrt(1+8*N)-1)/2, merror
                  fock = numpy.zeros((self.__nmos, self.__nmos), dtype=numpy.float64)
                  K = 0
                  for i in range(self.__nmos):
                      for j in range(i+1):
                          d = data[K]
                          fock[i,j] = d
                          fock[j,i] = d
                          K += 1
                  self.__fock = fock
            # Canonical Fock matrix
            elif  key == 'fckc':
                  merror = 'ncmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Canonical Fock matrix ]!'
                  assert self.__ncmos == (math.sqrt(1+8*N)-1)/2, merror
                  fckc = numpy.zeros((self.__ncmos, self.__ncmos), dtype=numpy.float64)
                  K = 0
                  for i in range(self.__ncmos):
                      for j in range(i+1):
                          d = data[K]
                          fckc[i,j] = d
                          fckc[j,i] = d
                          K += 1
                  self.__fckc = fckc
            # Fock matrix first derivatives
            elif  key == 'fock1':
                  merror = 'nbasis and nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Fock matrix - first derivatives ]!'
                  assert self.__nmos == (math.sqrt(1+8*(N/self.__nmodes))-1)/2, merror
                  fock1 = numpy.zeros((self.__nmodes, self.__nmos, self.__nmos), dtype=numpy.float64)
                  data = data.reshape(self.__nmodes, N/self.__nmodes)
                  for i in range(self.__nmodes):
                      K = 0
                      for j in range(self.__nmos):
                          for k in range(j+1):
                              d = data[i,K]
                              fock1[i,j,k] = d
                              fock1[i,k,j] = d
                              K += 1
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                  fock1= temp * fock1
                  self.__fock1 = fock1
            # Canonical Fock matrix first derivatives
            elif  key == 'fckc1':
                  merror = 'nbasis and ncmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Canonical Fock matrix - first derivatives ]!'
                  assert self.__ncmos == (math.sqrt(1+8*(N/self.__nmodes))-1)/2, merror
                  fckc1 = numpy.zeros((self.__nmodes, self.__ncmos, self.__ncmos), dtype=numpy.float64)
                  data  = data.reshape(self.__nmodes, N/self.__nmodes)
                  for i in range(self.__nmodes):
                      K = 0
                      for j in range(self.__ncmos):
                          for k in range(j+1):
                              d = data[i,K]
                              fckc1[i,j,k] = d
                              fckc1[i,k,j] = d
                              K += 1
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                  fckc1= temp * fckc1
                  self.__fckc1 = fckc1
            # AO to LMO transformation matrix
            elif  key == 'vecl':
                  merror = 'nmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->LMO matrix ]!'
                  assert self.__nmos == N/self.__nbasis, merror
                  data = data.reshape(self.__nmos, self.__nbasis)
                  self.__vecl = data
            # AO to CMO transformation matrix
            elif  key == 'vecc':
                  merror = 'ncmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->CMO matrix ]!'
                  assert self.__ncmos == N/self.__nbasis, merror
                  data = data.reshape(self.__ncmos, self.__nbasis)
                  self.__vecc = data
            # AO to LMO transformation matrix first derivatives
            elif  key == 'vecl1':
                  merror = 'nmodes, nmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->LMO matrix - first derivatives ]!'
                  assert self.__nmos == N/(self.__nbasis*self.__nmodes)
                  data = data.reshape(self.__nmodes, self.__nmos, self.__nbasis)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                  data = temp * data
                  self.__vecl1 = data
            # AO to CMO transformation matrix first derivatives
            elif  key == 'vecc1':
                  merror = 'nmodes, ncmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->CMO matrix - first derivatives ]!'
                  assert self.__ncmos == N/(self.__nbasis*self.__nmodes)
                  data = data.reshape(self.__nmodes, self.__ncmos, self.__nbasis)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = numpy.sqrt(self.__redmass)[:, numpy.newaxis, numpy.newaxis]
                  data = temp * data
                  self.__vecc1 = data
            # ------------------------------------ DMS --------------------------------------------
            # OPDM-Alpha
            elif  key == 'opdm':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__opdm = data
            # Ca(OCC)-Alpha
            elif  key == 'caocc':
                  merror = ''
                  assert N == self.__nbasis*self.__naocc, merror
                  data = data.reshape(self.__nbasis, self.__naocc)
                  self.__caocc = data
            # Ca(OCC)-Alpha
            elif  key == 'cavir':
                  merror = ''
                  assert N == self.__nbasis*self.__navir, merror
                  data = data.reshape(self.__nbasis, self.__navir)
                  self.__cavir = data
            # DMS 10
            elif  key == 'dms_10':
                  merror = ''
                  assert N == 3*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3)
                  self.__dms_10 = data
            # DMS 20
            elif  key == 'dms_20':
                  merror = ''
                  assert N == 9*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3, 3)
                  self.__dms_20 = data
            # DMS 01
            elif  key == 'dms_01':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_01 = data
            # DMS 02
            elif  key == 'dms_02':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_02 = data
            # DMS 03
            elif  key == 'dms_03':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_03 = data
            # DMS 04
            elif  key == 'dms_04':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_04 = data
            # DMS 05
            elif  key == 'dms_05':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_05 = data
            # DMS 06
            elif  key == 'dms_06':
                  merror = ''
                  assert N == self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis)
                  self.__dms_06 = data
            # DMS 11
            elif  key == 'dms_11':
                  merror = ''
                  assert N == 3*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3)
                  self.__dms_11 = data
            # DMS 12
            elif  key == 'dms_12':
                  merror = ''
                  assert N == 3*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3)
                  self.__dms_12 = data
            # DMS 21
            elif  key == 'dms_21':
                  merror = ''
                  assert N == 9*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3, 3)
                  self.__dms_21 = data
            # DMS 22
            elif  key == 'dms_22':
                  merror = ''
                  assert N == 9*self.__natoms*self.__nbasis**2, merror
                  data = data.reshape(self.__nbasis, self.__nbasis, self.__natoms, 3, 3)
                  self.__dms_22 = data
           
            # other not-programmed section
            else:
                raise Exception("Non-standard section name '%s' detected! Check the format of your file!"%key)
        return

    # --------------------------------------------------------- #
    #            W R I T I N G    P R O C E D U R E S           #
    # --------------------------------------------------------- #
    
    def _write_preambule(self,file):
        """write the preambule of the parameter file"""
        log    = ' %s\n' % self.__sec_names['mol'].ljust(40)
        log   += '   name       = %s\n'    % self.__name
        log   += '   basis      = %s/%s\n' %(self.__method, self.__basis)
        log   += '   charge     = %s\n'    % self.__charge
        log   += '   multiplic  = %s\n'    % self.__multiplicity
        log   += '   natoms     = %s\n'    % self.__natoms
        log   += '   nbasis     = %s\n'    % self.__nbasis
        log   += '   nmodes     = %s\n'    % self.__nmodes
        if self.__ndma is not None:
           log+= '   ndma       = %s\n'    % self.__ndma
        if self.__npol is not None:
           log+= '   npol       = %s\n'    % self.__npol
        log   += '   nmos       = %s\n'    % self.__nmos
        log   += '   naocc      = %s\n'    % self.__naocc
        log   += '   nbocc      = %s\n'    % self.__nbocc
        log   += '   navir      = %s\n'    % self.__navir
        log   += '   nbvir      = %s\n'    % self.__nbvir
        log   += '   ncmos      = %s\n'    % self.__ncmos
        if self.__mode is not None:
           log+= '   mode       = %s\n'    % ','.join(['%3d' % x for x in self.__mode])
        log   += ' \n'
        file.write(log)
        return
    
    def _write_pos(self,file):
        """write atomic coordinates"""
        natoms = self.__pos.shape[0]
        N = natoms * 3
        log = ' %s %s= %d\n' % (self.__sec_names['pos'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            for j in range(3):
                log+= "%20.10E" % self.__pos[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_atoms(self,file):
        """write atomic symbols"""
        natoms = self.__atoms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['atoms'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            log+= "%20s" % self.__atoms[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_atno(self,file):
        """write atomic numbers"""
        natoms = self.__atno.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['atno'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            log+= "%20d" % self.__atno[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_atms(self,file):
        """write atomic masses"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['atms'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            log+= "%20.10E" % self.__atms[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_origin(self,file):
        """write DMTP origin coordinates"""
        nsites = self.__origin.shape[0]
        N = nsites * 3
        log = ' %s %s= %d\n' % (self.__sec_names['origin'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nsites):
            for j in range(3):
                log+= "%20.10E" % self.__origin[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_esp(self, file):
        """write ESP charges"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['esp'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            log+= "%20.10E" % self.__esp[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_chlpg(self, file):
        """write ESP charges"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['chlpg'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(natoms):
            log+= "%20.10E" % self.__chlpg[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_rdma(self, file):
        """write DMTP center coordinates"""
        ndma = self.__rdma.shape[0]
        N = ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['rdma'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ndma):
            for j in range(3):
                log+= "%20.10E" % self.__rdma[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dmac(self, file):
        """write DMTP charges"""
        ndma = self.__dmac.shape[0]
        N = ndma
        log = ' %s %s= %d\n' % (self.__sec_names['dmac'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ndma):
            log+= "%20.10E" % self.__dmac[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmac1(self, file):
        """write DMTP charges 1st der wrt nmodes"""
        nmodes, ndma = self.__dmac1.shape
        N = ndma * nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['dmac1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ndma):
                log+= "%20.10E" % self.__dmac1[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmac2(self, file):
        """write DMTP charges 2nd derivatives wrt <self.__mode> 1D array"""
        ndma = self.__dmac2.shape[1]
        N = self.__nmodes_sder * ndma
        log = ' %s %s= %d\n' % (self.__sec_names['dmac2'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for m in range(self.__nmodes_sder):
            for i in range(ndma):                
                log+= "%20.10E" % self.__dmac2[m,i]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad(self, file):
        """write DMTP dipoles"""
        ndma = self.__dmad.shape[0]
        N = ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ndma):
            for j in range(3):
                log+= "%20.10E" % self.__dmad[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad1(self, file):
        """write DMTP dipoles 1st der wrt nmodes"""
        nmodes, ndma, c = self.__dmad1.shape
        N = nmodes * ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ndma):
                for k in range(3):
                    log+= "%20.10E" % self.__dmad1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad2(self, file):
        """write DMTP dipoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmad2.shape[1]
        N = self.__nmodes_sder * ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad2'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for m in range(self.__nmodes_sder):
            for i in range(ndma):                      
                for j in range(3):
                    log+= "%20.10E" % self.__dmad2[m,i,j]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq(self, file):
        """write DMTP quadrupoles"""
        ndma = self.__dmaq.shape[0]
        N = ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ndma):
            for j in range(6):
                log+= "%20.10E" % self.__dmaq[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq1(self, file):
        """write DMTP quadrupoles 1st der wrt nmodes"""
        nmodes, ndma, c = self.__dmaq1.shape
        N = nmodes * ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ndma):
                for k in range(6):
                    log+= "%20.10E" % self.__dmaq1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq2(self, file):
        """write DMTP quadrupoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmaq2.shape[1]
        N = self.__nmodes_sder * ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq2'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for m in range(self.__nmodes_sder):
            for i in range(ndma):                      
                for j in range(6):
                    log+= "%20.10E" % self.__dmaq2[m,i,j]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao(self, file):
        """write DMTP octupoles"""
        ndma = self.__dmao.shape[0]
        N = ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ndma):
            for j in range(10):
                log+= "%20.10E" % self.__dmao[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao1(self, file):
        """write DMTP octupoles 1st der wrt nmodes"""
        nmodes, ndma, c = self.__dmao1.shape
        N = nmodes * ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ndma):
                for k in range(10):
                    log+= "%20.10E" % self.__dmao1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao2(self, file):
        """write DMTP octupoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmao2.shape[1]
        N = self.__nmodes_sder * ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao2'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for m in range(self.__nmodes_sder):
            for i in range(ndma):                      
                for j in range(10):
                    log+= "%20.10E" % self.__dmao2[m,i,j]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_rpol(self, file):
        """write polarizable center coordinates"""
        npol = self.__rpol.shape[0]
        N = npol * 3
        log = ' %s %s= %d\n' % (self.__sec_names['rpol'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(npol):
            for j in range(3):
                log+= "%20.10E" % self.__rpol[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_dpol(self, file):
        """write distributed polarizabilities"""
        nmos = self.__dpol.shape[0]
        N = nmos*9
        log = ' %s %s= %d\n' % (self.__sec_names['dpol'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmos):
            for j in [0,1,2]:
                for k in [0,1,2]:
                    log+= "%20.10E" % self.__dpol[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_dpol1(self, file):
        """write distributed polarizabilities first derivatives wrt modes"""
        nmodes = self.__dpol1.shape[0]
        nmos   = self.__dpol1.shape[1]
        N = nmodes*nmos*9
        log = ' %s %s= %d\n' % (self.__sec_names['dpol1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(nmos):
                for k in [0,1,2]:
                    for l in [0,1,2]:
                        log+= "%20.10E" % self.__dpol1[i,j,k,l]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dpoli(self, file):
        """write distributed polarizabilities wrt imaginary frequencies"""
        nmos = self.__dpoli.shape[0]
        N = nmos*9*12
        log = ' %s %s= %d\n' % (self.__sec_names['dpoli'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmos):
            for j in range(12):
                for k in [0,1,2]:
                    for l in [0,1,2]:
                        log+= "%20.10E" % self.__dpoli[i,j,k,l]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dpoli1(self, file):
        """write distributed polarizabilities wrt imag freq - first derivatives wrt modes"""
        nmodes = self.__dpoli1.shape[0]
        nmos   = self.__dpoli1.shape[1]
        N = nmodes*nmos*9*12
        log = ' %s %s= %d\n' % (self.__sec_names['dpoli1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(nmos):
                for k in range(12):
                    for l in [0,1,2]:
                        for m in [0,1,2]:
                            log+= "%20.10E" % self.__dpoli1[i,j,k,l,m]
                            if not n%5: log+= '\n'
                            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
   
    def _write_freq(self, file):
        """write harmonic frequencies"""
        nmodes = self.__freq.shape[0]
        N = nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['freq'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            log+= "%20.10E" % self.__freq[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_redmass(self, file):
        """write reduced masses"""
        nmodes = self.__redmass.shape[0]
        N = nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['redmass'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            log+= "%20.10E" % self.__redmass[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_lvec(self, file):
        """write mass-weighted eigenvectors"""
        nmodes, natoms, x = self.__lvec.shape 
        N = nmodes * natoms * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lvec'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(natoms):
                for k in [0,1,2]:
                    log+= "%20.10E" % self.__lvec[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_gijk(self, file):
        """write cubic anharmonic constants"""
        nmodes = self.__gijk.shape[0]
        N = nmodes*(nmodes+1)*(nmodes+2)/6
        log = ' %s %s= %d\n' % (self.__sec_names['gijk'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(i+1):
                for k in range(j+1):
                    log+= "%20.10E" % self.__gijk[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_lmoc(self, file):
        """write LMO centroid coordinates"""
        nmos = self.__lmoc.shape[0]
        N = nmos * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lmoc'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmos):
            for j in range(3):
                log+= "%20.10E" % self.__lmoc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_lmoc1(self, file):
        """write LMO centroid first derivatives wrt nmodes"""
        nmodes, nmos, n = self.__lmoc1.shape
        N = nmodes * nmos * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lmoc1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(nmos):
                for k in range(3):
                    log+= "%20.10E" % self.__lmoc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_fock(self, file):
        """write Fock matrix elements"""
        nmos = self.__fock.shape[0]
        N = (nmos**2 - nmos) / 2 + nmos
        log = ' %s %s= %d\n' % (self.__sec_names['fock'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmos):
            for j in range(i+1):
                log+= "%20.10E" % self.__fock[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_fckc(self, file):
        """write canonical Fock matrix elements"""
        ncmos = self.__fckc.shape[0]
        N = (ncmos**2 - ncmos) / 2 + ncmos
        log = ' %s %s= %d\n' % (self.__sec_names['fckc'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ncmos):
            for j in range(i+1):
                log+= "%20.10E" % self.__fckc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_fock1(self, file):
        """write Fock matrix element first derivatives wrt nmodes"""
        nmodes = self.__fock1.shape[0]
        nmos = self.__fock1.shape[1]
        N = nmodes * ( (nmos**2 - nmos) / 2 + nmos )
        log = ' %s %s= %d\n' % (self.__sec_names['fock1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(nmos):
                for k in range(j+1):
                    log+= "%20.10E" % self.__fock1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_fckc1(self, file):
        """write canonical Fock matrix element first derivatives wrt nmodes"""
        nmodes = self.__fckc1.shape[0]
        ncmos = self.__fckc1.shape[1]
        N = nmodes * ( (ncmos**2 - ncmos) / 2 + ncmos )
        log = ' %s %s= %d\n' % (self.__sec_names['fckc1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ncmos):
                for k in range(j+1):
                    log+= "%20.10E" % self.__fckc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_vecl(self, file):
        """write AO-LMO transformation matrix elements"""
        nmos, nbasis = self.__vecl.shape
        N = nmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecl'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmos):
            for j in range(nbasis):
                log+= "%20.10E" % self.__vecl[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_vecc(self, file):
        """write AO-CMO transformation matrix elements"""
        ncmos, nbasis = self.__vecc.shape
        N = ncmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecc'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(ncmos):
            for j in range(nbasis):
                log+= "%20.10E" % self.__vecc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_vecl1(self, file):
        """write AO-LMO transformation matrix element first derivatives wrt nmodes"""
        nmodes, nmos, nbasis = self.__vecl1.shape
        N = nmodes * nmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecl1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(nmos):
                for k in range(nbasis):
                    log+= "%20.10E" % self.__vecl1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_vecc1(self, file):
        """write AO-CMO transformation matrix element first derivatives wrt nmodes"""
        nmodes, ncmos, nbasis = self.__vecc1.shape
        N = nmodes * ncmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecc1'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for i in range(nmodes):
            for j in range(ncmos):
                for k in range(nbasis):
                    log+= "%20.10E" % self.__vecc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_opdm(self, file):
        """write OPDM"""
        N = self.__opdm.size
        log = ' %s %s= %d\n' % (self.__sec_names['opdm'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                log+= "%20.10E" % self.__opdm[a,b]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_caocc(self, file):
        """write CaOCC"""
        N = self.__caocc.size
        log = ' %s %s= %d\n' % (self.__sec_names['caocc'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__naocc):
                log+= "%20.10E" % self.__caocc[a,b]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_cavir(self, file):
        """write CaVIR"""
        N = self.__cavir.size
        log = ' %s %s= %d\n' % (self.__sec_names['cavir'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__navir):
                log+= "%20.10E" % self.__cavir[a,b]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_cbocc(self, file):
        """write CbOCC"""
        N = self.__cbocc.size
        log = ' %s %s= %d\n' % (self.__sec_names['cbocc'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbocc):
                log+= "%20.10E" % self.__cbocc[a,b]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_cbvir(self, file):
        """write CbVIR"""
        N = self.__cbvir.size
        log = ' %s %s= %d\n' % (self.__sec_names['cbvir'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbvir):
                log+= "%20.10E" % self.__cbvir[a,b]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_10(self, file):
        """write DMS of order 10"""
        N = self.__dms_10.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_10'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                        log+= "%20.10E" % self.__dms_10[a,b,i,x]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_20(self, file):
        """write DMS of order 20"""
        N = self.__dms_20.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_20'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                      for y in range(3):
                        log+= "%20.10E" % self.__dms_20[a,b,i,x,y]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_01(self, file):
        """write DMS of order 01"""
        N = self.__dms_01.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_01'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_01[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_02(self, file):
        """write DMS of order 02"""
        N = self.__dms_02.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_02'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_02[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_03(self, file):
        """write DMS of order 03"""
        N = self.__dms_03.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_03'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_03[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_04(self, file):
        """write DMS of order 04"""
        N = self.__dms_04.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_04'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_04[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_05(self, file):
        """write DMS of order 05"""
        N = self.__dms_05.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_05'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_05[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_06(self, file):
        """write DMS of order 06"""
        N = self.__dms_06.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_06'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                  log+= "%20.10E" % self.__dms_06[a,b]
                  if not n%5: log+= '\n'
                  n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_11(self, file):
        """write DMS of order 11"""
        N = self.__dms_11.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_11'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                        log+= "%20.10E" % self.__dms_11[a,b,i,x]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_22(self, file):
        """write DMS of order 22"""
        N = self.__dms_22.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_22'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                      for y in range(3):
                        log+= "%20.10E" % self.__dms_22[a,b,i,x,y]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_12(self, file):
        """write DMS of order 12"""
        N = self.__dms_12.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_12'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                        log+= "%20.10E" % self.__dms_12[a,b,i,x]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dms_21(self, file):
        """write DMS of order 21"""
        N = self.__dms_21.size
        log = ' %s %s= %d\n' % (self.__sec_names['dms_21'].ljust(40), 'N'.rjust(10), N)
        n = 1
        for a in range(self.__nbasis):
            for b in range(self.__nbasis):
                for i in range(self.__natoms):
                    for x in range(3):
                      for y in range(3):
                        log+= "%20.10E" % self.__dms_21[a,b,i,x,y]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

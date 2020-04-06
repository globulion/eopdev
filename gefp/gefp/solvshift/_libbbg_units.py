# --------------------------------------------------------------- #
#    U N I T S   A N D   S T A N D A R D   C O N V E R T E R S    #
# --------------------------------------------------------------- #

__all__=['UNITS','Atom']
__version__ = '3.0.1'

import numpy, math

class Atom:
      """Represents an atom. Stores, its symbol, atomic number, mass and vdW radius
      basing on a symbol given as an argument.
      
      Usage:
      
      Atom('symbol') 
      
      where symbol = Na, H, etc... 
      
      If no symbol Atom() a dummy atom is created"""
      
      def __init__(self,symbol='X'):
          self.symbol = symbol
          self.mass = self._mass(symbol)
          self.atno = self._atno(symbol)
          self.radius  = self._wdW_radius(symbol)
          
      def _mass(self,symbol='X'):
          """returns mass by symbol. In AMU.""" 
          
          masses = { 'H': 1.0078250,  'C':12.0000000,  'N':14.0030740,
                     'O':15.9949146,  'F':18.9984033,  'S':32.0666   ,
                    'Fe':55.847    , 'Na':22.9897697,  'X':0.0000    ,
                    'Cl':34.9688527,  'P':30.9737634, 'Li':7.0000    ,
                    'Mg':23.985045 }
          # isotopes
          masses.update({'D': 2.0141018, 'H(iso=2)': 2.0141018, 'Bq': 0.00000,
                         'C(iso=13)': 13.003355, 'N(iso=15)': 15.000109, })
          
          return masses[symbol]

      def _atno(self,symbol='X'):
          """returns atomic number by symbol. In Electron"""
          
          atomic_numbers = { 'H': 1, 'He': 2, 'C': 6,
                             'N': 7,  'O': 8, 'F': 9,
                             'S':16, 'Fe':26, 'Na': 11, 'X': 0,
                            'Cl':17,  'P':15, 'Li': 3,
                            'Mg':12, 'H(iso=2)': 1, 'C(iso=13)': 6, 'N(iso=15)': 7 }
          # isotopes
          atomic_numbers.update({'D': 1, 'H(iso=2)': 1, 'Bq': 0, 'C(iso=13)': 6, 'N(iso=15)': 7})
                             
          return atomic_numbers[symbol]

      def _wdW_radius(self,symbol='X'):
          """returns atomic vdW radii by symbol. In Bohr"""
          
          wdW_radius = { 'H': 0.5000, 'He': 0.5000, 'C': 2.0000,
                             'N': 2,  'O': 2.0000, 'F': 2.0000,
                             'S':2.0000, 'Fe':3.0000, 'Na': 2.0000, 'X': 0 ,
                            'Cl':1.9000, 'P' :1.1, 'Li': 1.7,
                            'Mg':2.9000 }
          # isotopes
          wdW_radius.update({'D': 0.5000, 'H(iso=2)': 0.5000, 'Bq': 0.0001,
                             'C(iso=13)': 2.0000, 'N(iso=15)': 2.0000})
                             
          return wdW_radius[symbol]

      def __repr__(self):
          """print me! I am a happy lovely atom!"""
          log=" Atom: %8s, atno: %6i, my mass is: %12.7f" %(self.symbol,
                                                            self.atno,
                                                            self.mass)
          return str(log)


class UNITS:
      """common unit contained and converter"""
      CmRecToHartree         = 4.556335281212E-06
      HartreeToCmRec         = 219474.63
      HzToAuAngFreq          = 2.41888433412929E-17
      BohrElectronToDebye    = 2.541746
      DebyeToBohrElectron    = 0.39343034276438327
      AngstromToBohr         = 1.889725989
      AngstromToNanometer    = 0.100000000
      NanometerToBohr        = 18.89725989 
      NanometerToAngstrom    = 10.00000000
      BohrToAngstrom         = 0.5291772086
      BohrToNanometer        = 0.05291772086      
      ElectronMassToAmu      = 1./1822.8889
      AmuToElectronMass      = 1822.8889      
      KmMToIrInt             = 1./42.2561
      IrIntToKmM             = 42.2561
      JouleToHartree         = 1./4.3597439E-18
      HartreeToJoule         = 4.3597439E-18
      AttoJouleToHartree     = 1./4.3597439
      HartreeToAttoJoule     = 4.3597439
      ToRedCubForceConst     = 9.85501E+06
      BohrToMeter            = 5.291772086E-11
      MeterToBohr            = 1.889725989E+10
      ElectronChargeToCoulomb= 1.602176487E-19
      CmRecToMeterRec        = 100.00000000000
      SpeedOfLight           = 2.99792458E+08           ### m/s
      AmuToKg                = 1.66053878E-27
      PlanckConstant         = 6.62606896E-34
      KgToAmu                = 1./1.66053878E-27
      ElectronMassToKg       = 9.1093800615056683E-31
      KgToElectronMass       = 1.0977695444125671E+30
      CmRecToHz              = 2.99792458E+10 * 2 * numpy.pi
      HzToCmRec              = 1./CmRecToHz
      CmToBohr               = 1.889725989E+08
      BohrToCm               = 5.291772086E-09
      CmRecToHartreePerHbar  = CmRecToHz * HzToAuAngFreq
      HartreePerHbarToCmRec  = 1./(CmRecToHz * HzToAuAngFreq)
      BohrElectronToEsuCm    = 2.541746E-18
      AuToStarkTunningRate   = 42.6810
      AuToQuadrupolarTunningRate    = 42.6810 * BohrToCm
      AuToOctupolarTunningRate      = 42.6810 * BohrToCm**2
      HartreeToKcalPerMole   = 627.509
      HartreeToElectronVolt  = 27.2116
      MuToG09Mu              = BohrElectronToDebye * math.sqrt(IrIntToKmM/ElectronMassToAmu) / BohrToAngstrom # converts [e a_0]/[a_0 * m_e^1/2] -> const* [D]/[Angs * AMU^1/2]
      
      ### masses of atoms in AMU
      #mass   = { 1: 1.0078250,  6:12.0000000,  7:14.0030740,
      #           8:15.9949146,  9:18.9984033, 16:32.0666   ,
      #          26:55.847    , 11:22.9897697}
      #mass   = { 1: 2.0141018,  6:12.0000000,  7:14.0030740, 
      #           8:15.9949146,  9:18.9984033, 16:32.0666   ,
      #          26:55.847  }
                
      ### atomic numbers of atoms in AMU
      atomic_numbers = { 'H': 1, 'He': 2, 'C': 6,
                         'N': 7,  'O': 8, 'F': 9,
                         'S':16, 'Fe':26, 'Na': 11,
                         'H(iso=2)':1, 'C(iso=13)':6, 'N(iso=15)':7 }

      ### non-global methods (object-specific unit conversion etc...)
      def Convert():
          """converts units of DMA object or other objects (not written yet)"""
          pass

      def __repr__(self):
          """print unit info"""
          log = "\n"
          log+= " ---- UNITS FROM LIBBBG LIBRARY ----\n\n"
          for mem in dir(UNITS):
              if not mem.startswith('_'):
                 log+= " %25s \n" % ( mem )

          return str(log)

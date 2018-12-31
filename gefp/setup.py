#!/usr/bin/python3
# -*- coding: utf-8 -*-

# ---------------------------------------------------- #
from setuptools import setup, find_packages, Extension
import sys
# ---------------------------------------------------- #

def check_python_version():
    python_version = sys.version_info
    print(" Python version is ", python_version[:3])
    if not python_version[:2] in [(3, 3), (3, 4), (3, 5), (3, 6), (3, 7)]: 
       print(" VersionError: GEFP-OEP requires Python 3.3 or higher!")
       sys.exit(-1)
    return

def Main(argv):
    #check_python_version()
    classifiers = ['Development Status :: 1 - Planning',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: Freeware',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Unix',
                   'Programming Language :: Python',
                   'Programming Language :: C++',
                   'Programming Language :: Fortran',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ]
    extensions  = []
    print(" Installing has started!")
    setup(name                         =  'GEFP-OEP'  , 
          version                      =  '1.0.0'     ,
          description                  =  'Generalized Effective Fragment Potentials (GEFP) derived from One-Electron Potential (OEP) Method',
          long_description             =   open("README.md").read(),
          author                       =  'Bartosz Błasiak'      ,
          author_email                 =  'blasiak.bartosz@gmail.com',
          url                          =  '',
          license                      =   open("LICENSE.md").read(),
          requires                     = ['numpy (>=1.10.0)', 'psi4', 'oepdev'],
          packages                     = find_packages(),
          classifiers                  =   classifiers,
          ext_modules                  =   extensions,
          install_requires             = [],
          test_suite                   =  '',
          tests_require                = [],
          )
          

if __name__ == "__main__": Main(sys.argv[1:])

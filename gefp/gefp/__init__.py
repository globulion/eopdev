#*-* coding: utf-8 *-*
"""
Generalized Effective Fragment Potential Package.

Bartosz Błasiak, Gundelfingen, Dec 2018
"""

from . import density
from . import basis
from . import math
from . import core

__author__ =  'Bartosz Błasiak'
__version__=  '1.0.1'

__all__    = ['core', 'basis', 'density', 'math']

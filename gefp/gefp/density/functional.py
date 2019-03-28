#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT and DFT Functional module.
 Bartosz Błasiak, Gundelfingen, Mar 2019
"""

import os
import sys
import math
import numpy
from abc import ABC, abstractmethod

__all__ = ["XCFunctional"]


class XCFunctional(ABC):
    """\
 The Exchange-Correlation DMFT functional.
"""
    default = 'hf'

    def __init__(self):
        super(XCFunctional, self).__init__()

    @classmethod
    def create(cls, name=default, **kwargs):
        """\
 Create a functional. Available functionals:
"""
        if   name.lower() == 'hf' :   xc_functional = HF_XCFunctional()
        elif name.lower() == 'mbb':   xc_functional = MBB_XCFunctional()
        else: raise ValueError("Chosen XC functional is not available! Mistyped?")
        return xc_functional

    @staticmethod
    @abstractmethod
    def fij(n): pass

    @staticmethod
    @abstractmethod
    def name(): pass

    @property
    def abbr(self): return default.upper()

    @classmethod
    def _correct_negative_occupancies(cls, n):
        "Remove negative values of occupancies."
        ns = n.copy()
        ns[ns<0.0] = 0.0
        return ns


class HF_XCFunctional(XCFunctional):
    """
 Hartree-Fock Exchange-Correlation Functional
"""
    def __init__(self):
        super(HF_XCFunctional, self).__init__()
    @staticmethod
    def fij(n): return numpy.outer(n, n)
    @staticmethod
    def name(): return "Hartree-Fock Functional XC for closed-shell systems"
    @property
    def abbr(self): return "HF"

class MBB_XCFunctional(XCFunctional):
    """
 Muller-Buijse-Baerends Exchange-Correlation Functional
"""
    def __init__(self):
        super(MBB_XCFunctional, self).__init__()
    @staticmethod
    def fij(n): 
        ns = numpy.sqrt(self.correct_negative_occupancies(n))
        return numpy.outer(ns, ns)
    @staticmethod
    def name(): return "Muller-Buijse-Baerends XC Functional for closed-shell systems"
    @property
    def abbr(self): return "MBB"




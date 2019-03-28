#!/usr/bin/python3
#*-* coding: utf-8 *-*
"""
 DMFT Solver module.
 Bartosz Błasiak, Gundelfingen, Mar 2019
"""

import os
import sys
import math
import numpy
from abc import ABC, abstractmethod
from .functional import XCFunctional

__all__ = ["DMFT"]


class DMFT(ABC):
    def __init__(self, wfn):
        self.wfn = wfn
        self.basis_name = 'sto-3g'.upper()
        super(DMFT, self).__init__()

    @classmethod
    def create(cls, wfn, algorithm='proj-d', **kwargs):
        if algorithm.lower() == 'proj-d':
           solver = DMFT_ProjD(wfn, **kwargs)
        else: raise ValueError("Chosen algorithm is not available! Mistyped?")
        return solver

    def run(self, xc_functional, verbose=True):
        if xc_functional.abbr.lower() != 'hf' and self.abbr.lower() == 'dmft-projd':
           print(""" Warning! The D-set is not Lipshitz with %s functional. 
 This will probably result in lack of convergence. Use P-set instead.""" % xc_functional.abbr.upper())

        if verbose:
           print(" Running %s:%s/%s" % (self.abbr, xc_functional.abbr, self.basis_name))

    @staticmethod
    @abstractmethod
    def name(): pass

    @abstractmethod
    def gradient(self): pass

    @property
    def abbr(self): pass


class DMFT_ProjD(DMFT):
    def __init__(self, wfn):
        super(DMFT_ProjD, self).__init__(wfn)


    @staticmethod
    def name(): return "DMFT with Projection on D Set"

    def gradient(self): pass

    @property
    def abbr(self): return "DMFT-ProjD"



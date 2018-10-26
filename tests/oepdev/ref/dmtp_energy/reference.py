#!/usr/bin/python
"""
 Interaction energy between two H2O molecules at HF/sto-3g.
 For test in OEPDev.
 
                  Created:      Gundelfingen, 26 Oct 2018
"""
import libbbg
d1 = libbbg.dma.DMA('w1.par')
d2 = libbbg.dma.DMA('w2.par')
e  = libbbg.utilities.get_elmtp(d1,d2,converter=1.0)
print """
 Interaction energy [A.U.]:
  * R-1  %10.6f
  * R-2  %10.6f
  * R-3  %10.6f
  * R-4  %10.6f
  * R-5  %10.6f
""" % e


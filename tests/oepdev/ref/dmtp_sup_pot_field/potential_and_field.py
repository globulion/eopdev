#!/usr/bin/python2.7
import dma, numpy

def Vr_dma(dma,Rb,is_full=False):
    """calculates electrostatic potential in point Rb 
from dma distribution.

COPIED FROM LIBBBG.UTILITIES
"""

    dmac = dma.copy()
    if not is_full:
       dmac.MAKE_FULL()
       dmac.MakeTraceless()
    Ra,qa,Da,Qa,Oa,Ha = dmac.DMA_FULL
    V=0
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
        V+=qa[i]/Rab
        V+=numpy.tensordot(Da[i],R,(0,0))/Rab**3 
        V+=numpy.tensordot(R,numpy.tensordot(Qa[i],R,(1,0)),(0,0))/Rab**5 
        V+=numpy.tensordot(R,numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7 
    return V

def ElectricField(dma,Rb,is_full=False):
    """calculates electrostatic field in point Rb
from dma distribution. Usage:
ElectricField(dma,R,full=False/True)
if not is_full - the dma object is turned into traceless object

COPIED FROM LIBBBG.UTILITIES
"""

    dmac = dma.copy()
    if not is_full:
       dmac.MAKE_FULL()
       dmac.MakeTraceless()
    Ra,qa,Da,Qa,Oa,Ha = dmac.DMA_FULL
    field=numpy.zeros(3,dtype=numpy.float64)
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
        field+= qa[i] * R / Rab**3
        
        field+= 3 * R * numpy.dot(R,Da[i]) / Rab**5
        field-= Da[i]/Rab**3
        
        t  =numpy.tensordot(R,numpy.tensordot(Qa[i],R,(0,0)),(0,0))
        field+= 5* t * R / Rab**7
        field-= 2* numpy.tensordot(Qa[i],R,(0,0)) / Rab**5
        
        c=numpy.tensordot(R,numpy.tensordot(Oa[i],R,(0,0)),(0,0))
        g=numpy.tensordot(R,c,(0,0))
        field+= 7 * g * R / Rab**9
        field-= 3 * c / Rab**7
        
    return field

dmtp = dma.DMA('w1.par')
R = numpy.array([2.0,-6.0, -2.0])

potential = Vr_dma(dmtp, R)
field = ElectricField(dmtp, R)

print(" Pot= ", potential)
print(" Efield=", field)

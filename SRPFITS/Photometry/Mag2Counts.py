""" Utility functions and classes for SRP

Context : SRP
Module  : Photometry
Version : 1.0.0
Author  : Stefano Covino
Date    : 15/02/2012
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : 
    
History : (15/02/2012) First version.

"""


import math
import numpy


def Mag2Counts (mag, emag, zp=0.0):
    cnt = 10**(-0.4*(numpy.array(mag)-zp))
    ecnt = (numpy.array(emag)*cnt)/(2.5/math.log(10.0))
    return cnt, ecnt
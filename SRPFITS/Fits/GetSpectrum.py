""" Utility functions and classes for SRP

Context : SRP
Module  : Spectroscopy
Version : 1.0.1
Author  : Stefano Covino
Date    : 18/05/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

History : (22/03/2013) First version.
        : (18/05/2017) Minor update.
"""

import numpy
from SRPFITS.Fits.GetData import GetData
from SRPFITS.Fits.GetHeaderValue import GetHeaderValue



def GetSpectrum (filename, extension=0):
    data = GetData(filename, extension)[0]
    npix = GetHeaderValue(filename,'NAXIS1',extension)[0]
    refpix = GetHeaderValue(filename,'CRPIX1',extension)[0]
    reflmb = GetHeaderValue(filename,'CRVAL1',extension)[0]
    refdl = GetHeaderValue(filename,'CDELT1',extension)[0]
    if data != None and refpix != None and reflmb != None and refdl != None:
        pxl = numpy.linspace(1,npix,npix)
        lmbd = (pxl-refpix)*refdl+reflmb
        return lmbd, data
    else:
        return None, None

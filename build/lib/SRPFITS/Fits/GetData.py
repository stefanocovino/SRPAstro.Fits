""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.1.0
Author  : Stefano Covino
Date    : 31/07/2015
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (21/05/2010) First version.
        : (31/07/2015) python3 porting.
"""

from astropy.io import fits
from . import FitsConstants

def GetData (fitsfile, extension=0):
    try:
        hdr = fits.open(fitsfile)
    except IOError:
        return None,FitsConstants.FitsFileNotFound
    try:
        dataval = hdr[extension].data
    except IndexError:
        return None,FitsConstants.FitsDataSetNotFound
    hdr.close()
    return dataval,FitsConstants.FitsDataSetFound
    

    

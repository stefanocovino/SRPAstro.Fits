""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.3.0
Author  : Stefano Covino
Date    : 31/07/2015
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (21/05/2010) First version.
        : (14/02/2012) Possibility to read from other extensions.
        : (31/07/2015) python3 porting.
"""

from astropy.io import fits
from . import FitsConstants

def GetHeaderValue (fitsfile, header, ext=0):
    try:
        hdr = fits.open(fitsfile)
    except IOError:
        return None,FitsConstants.FitsFileNotFound
    hdr[ext].verify('silentfix+ignore')
    try:
        headval = hdr[ext].header[header]
    except KeyError:
        return None,FitsConstants.FitsHeaderNotFound
    hdr.close()
    return headval,FitsConstants.FitsHeaderFound
    

    

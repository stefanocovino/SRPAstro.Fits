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

History : (23/05/2010) First version.
        : (14/02/2012) Possibility to read extensions.
        : (16/12/2014) Manage wrong FITS headers.
        : (31/07/2015) python3 porting.
"""

from astropy.io import fits
from . import FitsConstants

def GetHeader (fitsfile, ext=0):
    try:
        hdr = fits.open(fitsfile)
    except IOError:
        return None,FitsConstants.FitsFileNotFound
    hdr[ext].verify('silentfix+ignore')
    heder = hdr[ext].header
    hdr.close()
    return heder,FitsConstants.FitsHeaderFound
    

    

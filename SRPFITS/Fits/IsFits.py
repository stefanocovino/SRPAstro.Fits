""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.1.1
Author  : Stefano Covino
Date    : 18/05/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (27/05/2010) First version.
        : (31/07/2015) python3 porting.
        : (18/05/2017) Close open file anyway.
"""

from astropy.io import fits
import tempfile
from . import FitsConstants

def IsFits(fitsfile):
    try:
        hdr = fits.open(fitsfile)
        fitst = True
    except (IOError,OSError):
        hdr = tempfile.TemporaryFile()
        fitst = False
    finally:
        hdr.close()
        return fitst

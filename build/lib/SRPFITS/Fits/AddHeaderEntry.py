""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.4.0
Author  : Stefano Covino
Date    : 08/07/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : AddHeaderEntry (fitsfile, keylist, entrylist, commentlist, outfilename=None)
            "fitfile" is the FITS file name.
            "keylist" is a list of keywords to be added to the FITS file headers.
            "entrylist" is a list of the same length than "keylist" with keywords values.
            "commentlist" is a list of the same length of "keylist" with comments (units, etc.)
                to the new headers.
            "outfilename" optional filename for output. Else input file file is overwritten.
            
            Function returns two values: (res, code). If res is False code reports the problem 
                (codes are in SRPFITS.FitsConstants). Else res is True.

Remarks :

History : (12/10/2010) First version.
        : (25/04/2011) Output and input file names can be different.
        : (27/04/2011) Always list in output. Keywords longer than 8 characters.
        : (14/02/2012) Possibilty to write in other extensions.
        : (31/03/2013) Value considered float unless it is not possible to convert it.
        : (31/07/2015) python3 porting.
        : (18/05/2017) astropy.io.fits.
        : (08/07/2021) Better management of the verify options in creating output FITS files.
"""

import warnings

from astropy.io import fits
from . import FitsConstants

def AddHeaderEntry (fitsfile, keylist, entrylist, commentlist, outfilename=None, ext=0):
    try:
        hdr = fits.open(fitsfile)
    except IOError:
        return False,FitsConstants.FitsFileNotFound
    heder = hdr[ext].header    
    for i,l,m in zip(keylist,entrylist,commentlist):
        try:
            lf = float(l)
        except ValueError:
            lf = l
        try:
            heder.set(i,lf,m)
        except ValueError:
            heder.set(i.split()[0][:8],l,m)
    #
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    warnings.filterwarnings('ignore', category=ResourceWarning, append=True)
    if outfilename == None:
        hdr.writeto(fitsfile,overwrite=True,output_verify='ignore')
    else:
        hdr.writeto(outfilename,overwrite=True,output_verify='ignore')
    warnings.resetwarnings() 
    warnings.filterwarnings('always', category=UserWarning, append=True)
    warnings.filterwarnings('always', category=ResourceWarning, append=True)
    return True,FitsConstants.FitsOk

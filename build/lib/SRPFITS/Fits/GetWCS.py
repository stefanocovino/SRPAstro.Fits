""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 2.1.0
Author  : Stefano Covino
Date    : 30/11/2020
E-mail  : stefano.covino@inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (15/07/2010) First version.
        : (20/12/2010) Possibility of badly formatted ECS header included.
        : (23/08/2013) NUMPY_MODE in wcs.
        : (16/12/2014) Manage wrong FITS headers.
        : (07/03/2020) astropy.wcs
		: (20/11/2020) astLib WCS again available.
"""



from astropy.wcs import WCS
from SRPFITS.Fits import FitsConstants
from SRPFITS.Fits.IsFits import IsFits
from .GetHeader import GetHeader
import astLib.astWCS as aw



def GetWCS (fitsfile,mode='apy'):
    if mode == 'apy':
        if IsFits(fitsfile):
            wcs = WCS(fitsfile)
            return wcs,FitsConstants.FitsWCSFound
        else:
            return None,FitsConstants.FitsFileNotFound
    else:
        heder = GetHeader(fitsfile)
        if heder[1] == FitsConstants.FitsHeaderFound:
            try:
                wcs = aw.WCS(heder[0],mode='pyfits')
            except IOError:
                return None,FitsConstants.FitsFileNotFound
            except ValueError:
                return None,FitsConstants.FitsWCSNotKnown
            return wcs,FitsConstants.FitsWCSFound
        else:
            return None,FitsConstants.FitsFileNotFound

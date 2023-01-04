""" Utility functions and classes for SRP

Context : SRP
Module  : Frames.py
Version : 2.0.0
Author  : Stefano Covino
Date    : 14/05/2020
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (23/05/2010) First version.
		: (14/05/2020) astropy version.
"""

import astLib.astWCS as astWCS
from astropy.wcs import WCS as aWCS
import math


def Pixel2WCS (fitsheader, poslist, mode='apy'):
    outlist = []
    if mode == 'astlib':
        WCS = astWCS.WCS(fitsheader, mode="pyfits")
        for i in poslist:
            outlist.append(WCS.pix2wcs(i[0],i[1]))
        return outlist
    elif mode == 'apy':
        w = aWCS(fitsheader)
        for i in poslist:
            outlist.append(w.wcs_pix2world(i[0],i[1],0,ra_dec_order=True))
        return outlist

    
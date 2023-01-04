""" Utility functions and classes for SRP

Context : SRP
Module  : Frames.py
Version : 1.1.0
Author  : Stefano Covino
Date    : 21/11/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (22/09/2010) First version.
        : (21/11/2017) Update wioth astropy.
"""

import astLib.astWCS as astWCS
from astropy.wcs import WCS as aWCS
import math


def WCS2Pixel (fitsheader, poslist, mode='apy'):
    outlist = []
    #
    if mode == 'astlib':
        WCS = astWCS.WCS(fitsheader, mode="pyfits")
        for i in poslist:
            outlist.append(WCS.wcs2pix(i[0],i[1]))
        return outlist
    elif mode == 'apy':
        w = aWCS(fitsheader)
        for i in poslist:
            outlist.append(w.wcs_world2pix(i[0],i[1],0,ra_dec_order=True))
        return outlist

""" Utility functions and classes for SRP

Module  : Frames.py
Version : 1.1.0
Author  : Stefano Covino
Date    : 07/03/2020
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (31/08/2012) First version.
        : (15/05/2017) Minor update.
        : (07/03/2020) astropy wcs tools.
"""

from SRPFITS.Fits.GetWCS import GetWCS


def getCenterRADEC (fitsfile):
    wcs = GetWCS(fitsfile)[0]
    if wcs != None:
        xpix, ypix = wcs.pixel_shape
        return wcs.all_pix2world(xpix/2.,ypix/2.,1,ra_dec_order=True)
    else:
        return None

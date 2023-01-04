""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.0.0
Author  : Stefano Covino
Date    : 01/04/2015
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : input pars are astropy.io.fits.hdu.table.BinTableHDU
        :   t = astropy.io.fits.open('test.fits')
        :   par is t[1]

History : (01/04/2015) First version.

"""

from astropy.io import fits
#import FitsConstants

def FitsTabsAppend(orgtbl,apptbl):
    nrows1 = orgtbl.data.shape[0]
    nrows2 = apptbl.data.shape[0]
    nrows = nrows1 + nrows2
    hdu = fits.BinTableHDU.from_columns(orgtbl.columns, nrows=nrows)
    for colname in orgtbl.columns.names:
        try:
            hdu.data[colname][nrows1:] = apptbl.data[colname]
        except KeyError:
            return None
    return hdu
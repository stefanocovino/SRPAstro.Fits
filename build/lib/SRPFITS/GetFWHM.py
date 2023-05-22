""" Utility functions and classes for SRP

Context : SRP
Module  : SRPGW
Version : 1.1.0
Author  : Stefano Covino
Date    : 22/05/2023
E-mail  : stefano.covino@inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : GetApPhot

Remarks :

History : (18/12/2016) First version.
        : (21/12/2016) Minor bug.
        : (03/07/2017) size must be integer in computations.
        : (22/05/2023) Deprecated numpy feature.
"""

import numpy as np
from astropy.io.fits import getdata



def GetFWHM(x,y,fname,size=20):
    fwhm = []
    data = getdata(fname)
    szint = int(size)
    for ii,jj in zip(x,y):
        a = np.rint(ii)-1
        b = np.rint(jj)-1
        i = a.astype(int)
        l = b.astype(int)
        image = data[l-szint:l+szint,i-szint:i+szint]
        dat=image.flatten()
        try:
            maxi = image.max()
            floor = np.median(image)
            height = maxi - floor
            #
            fwhm.append(np.sqrt(np.sum(image>floor+height/2.).flatten())[0])
        except ValueError:
            fwhm.append(-99.)
    return fwhm

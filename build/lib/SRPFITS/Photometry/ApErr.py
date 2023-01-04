""" Utility functions and classes for SRP

Context : SRP
Module  : Photometry
Version : 1.0.0
Author  : Stefano Covino
Date    : 19/04/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : 
            gain: gain in "electron / pixel data units"
            NA: number of pixel in source aperture
            NB: number of pixels in background annulus
            ni: depth - of - coverage at pixel i; N i > 1 if image is a co-add
            SI: signal per pixel in image data units
            BM: estimated background per pixel in annulus (either mean or median) 
                    if mean than k = 1
                    else if median  k = np.pi/2
                    else if no background k = 0
            SB2: variance in sky background annulus in [image units]^2 / pixel

Remarks : adapted from a text by F. Masci, v. 1.0.
    
History : (19/04/2017) First version.

"""


import numpy as np


def ApErr (SI, NA, NB, BM, SB2, ni=1, gain=1, k=1):
    ft = (1./gain)*NA*(SI-BM)/ni
    st = (NA+k*NA**2/NB)*SB2
    ftw = np.where(ft < 0, 0., ft)
    return ftw+st

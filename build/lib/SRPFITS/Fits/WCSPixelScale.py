""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.0.0
Author  : Stefano Covino
Date    : 14/01/2020
E-mail  : stefano.covino@inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : returns pixel scales in degrees

History : (14/01/2021) First version.
"""


import numpy as np

def WCSPixelScale(WCS):
    scl = WCS.proj_plane_pixel_scales()
    return scl

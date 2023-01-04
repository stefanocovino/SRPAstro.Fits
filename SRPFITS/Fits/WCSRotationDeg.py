""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.0.0
Author  : Stefano Covino
Date    : 13/01/2020
E-mail  : stefano.covino@inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : returns angle in degrees

History : (13/01/2021) First version.
"""


import numpy as np

from SRPFITS.Fits.WCSPixelScale import WCSPixelScale


def WCSRotationDeg(WCS):
    scl = WCSPixelScale(WCS)
    xs = scl[0].value
    #
    ang = np.degrees(np.arcsin(WCS.pixel_scale_matrix[1][0]/xs))
    return ang
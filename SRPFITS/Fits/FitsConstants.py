""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.2.1
Author  : Stefano Covino
Date    : 27/04/2011
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (24/05/2010) First version.
        : (15/07/2010) New constants.
        : (30/08/2010) Astrometry headers.
        : (04/09/2010) New headers.
        : (27/09/2010) Eclipse constants.
        : (29/09/2010) No more eclipse constants.
        : (17/11/2010) Comments for WCS header.
        : (20/12/2010) Not known WCS entry.
        : (27/04/2011) Constant for "no problem".

"""

# Header
FitsFileNotFound    =   -1
FitsHeaderNotFound  =   -2
FitsHeaderFound     =   1
FitsOk              =   0

# Data
FitsDataSetFound    =   2
FitsDataSetNotFound =   -3

# WCS
FitsWCSFound    = 3
FitsWCSNotKnown = -3


# Finding-chart output
FCExt = '_fc.png'

# Source output
RegData = '_src.dat'
SkyData = '_src.skycat'

# Size Headers
NAxis       =   'NAXIS'

# Astrometry Headers
CRPIX1  =   'CRPIX1'
CRPIX2  =   'CRPIX2'
CRVAL1  =   'CRVAL1'
CRVAL2  =   'CRVAL2'
CTYPE1  =   'CTYPE1'
CTYPE2  =   'CTYPE2'
RAType  =   'RA---TAN'           
DECType =   'DEC--TAN'    
CDELT1  =   'CDELT1'
CDELT2  =   'CDELT2'
PC11    =   'PC1_1'
PC12    =   'PC1_2'
PC21    =   'PC2_1'
PC22    =   'PC2_2'
CD11    =   'CD1_1'
CD12    =   'CD1_2'
CD21    =   'CD2_1'
CD22    =   'CD2_2'
RADECSys    =   'RADECSYS'
RADECSysVal =   'FK5'
AngUnit1    =   'CUNIT1'
AngUnit2    =   'CUNIT2'
AngUnit1Val =   'Degrees '
AngUnit2Val =   'Degrees '
InAxesUnit  =   'PIXEL'
RA      =   'RA'
DEC     =   'DEC'
#
CRPIX1Cmt   =   'reference X pixel'
CRPIX2Cmt   =   'reference Y pixel'
CRVALCmt    =   'sky coords at reference pixel'
CTYPECmt    =   'sky projection'
CDELT1Cmt   =   'increment for X pixel'
CDELT2Cmt   =   'increment for Y pixel'
PCCmt       =   'rotation matrix values'
AngUnCmt    =   'sky coords units'
#
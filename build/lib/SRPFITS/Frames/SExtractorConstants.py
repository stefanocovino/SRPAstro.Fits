""" Utility functions and classes for SRP

Context : SRP
Module  : Frames
Version : 1.1.0
Author  : Stefano Covino
Date    : 14/05/2020
E-mail  : stefano.covino@inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (28/09/2010) First version.
        : (01/10/2010) Sextractor command names.
        : (05/02/2012) Better path management.
        : (18/08/2013) ROS2 parameters.
        : (30/05/2017) Path for FITS package.
        : (14/05/2020) New constant.
"""

import os.path

import SRP
from SRPFITS.SRPFITSPath import SRPFITSPath


# SExtractor parameter sets for photometry

BasePath        = os.path.join(SRPFITSPath(),'Data','SExtractor')

SexFName        = ("SRP.param", "SRP.conv", "SRP.nnw", "SRP.sex")

GenParSet       = ("SRPParamIn", "SRPConvIn", "SRPNnwIn", "SRPSexIn")

SRPSEXPARDICT   = {"LS-LASC":('SRPParamIn','SRPConvLSLASCIn','SRPNnwIn','SRPSexLSLASCIn'),
                    "REM-ROSS":('SRPParamIn','SRPConvIn','SRPNnwIn','SRPSexREMROSSIn'),
                    "REM-ROS2":('SRPParamIn','SRPConvIn','SRPNnwIn','SRPSexREMROS2In'),
                    "REM-REMIR":('SRPParamIn','SRPConvIn','SRPNnwIn','SRPSexREMREMIRIn'),
                    "TNG-LRS":('SRPParamIn','SRPConvIn','SRPNnwIn','SRPSexTNGLRSIn')}


# SExtractor commands
SRPsex	       = "sex"
SRPsex_cyg     = "sex.exe"
SRPsex_new     = "sextractor"

#

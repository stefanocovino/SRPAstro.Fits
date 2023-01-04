""" Utility functions and classes for SRP

Context : SRP
Module  : Spectroscopy
Version : 1.0.0
Author  : Stefano Covino
Date    : 02/04/2011
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

History : (02/04/2013) First version.

"""

import numpy



def GetSpectrumPosition (data):
    if not 1 <= data.ndim <= 2:
        return None
    elif data.ndim == 2:
        datasum = data.sum(axis=1)
    else:
        datasum = data
    num = datasum*list(range(len(datasum)))
    return num.sum()/datasum.sum()
        
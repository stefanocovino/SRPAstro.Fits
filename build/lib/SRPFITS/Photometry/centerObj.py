""" Utility functions and classes for SRP

Context : SRP
Module  : Photometry
Version : 1.0.1
Author  : Stefano Covino
Date    : 19/04/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : 
    
History : (10/11/2015) First version.
        : (19/04/2017) SRPFITS added.

"""

import numpy as np
from SRPFITS.Photometry.centerMoment import centerMoment
from photutils.centroids import centroid_com
from photutils.centroids import centroid_1dg
from photutils.centroids import centroid_2dg



def centerObj (table, x, y, r, mode='my'):
    if mode == 'my':
        return centerMoment(table, x, y, r)
    else:
        if np.floor(x-r) < 1:
            xmin = 1
        else:
            xmin = int(np.floor(x-r))
        if np.ceil(x+r) > table.shape[1]:
            xmax = table.shape[1]
        else:
            xmax = int(np.ceil(x+r))
        if np.floor(y-r) < 1:
            ymin = 1
        else:
            ymin = int(np.floor(y-r))
        if np.ceil(y+r) > table.shape[0]:
            ymax = table.shape[0]
        else:
            ymax = int(np.ceil(y+r))
        #
        field = table[(ymin-1):ymax,(xmin-1):xmax]
        if mode == 'apymom':
            xr,yr = centroid_com(field)
            return xr+xmin, yr+ymin
        if mode == 'apy1dg':
            xr,yr = centroid_1dg(field)
            return xr+xmin, yr+ymin
        if mode == 'apy2dg':
            xr,yr = centroid_2dg(field)
            return xr+xmin, yr+ymin
        else:
            return x,y

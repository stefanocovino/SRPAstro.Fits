""" Utility functions and classes for SRP

Context : SRP
Module  : Photometry
Version : 1.0.3
Author  : Stefano Covino
Date    : 07/03/2020
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks : 
    
History : (18/11/2015) First version.
        : (19/04/2017) SRPFITS added.
        : (12/05/2017) Import corrected.
        : (07/03/2020) Non-tuple sequence.
"""

import numpy
from SRPFITS.Photometry.getBackground import getBackground
from SRPFITS.Photometry.MinMax import MinMax



def centerMoment (table, x, y, r):
    bck = getBackground(table,x,y,0,r)
    minmax = MinMax(table,x,y,r)
    if numpy.floor(x-r) < 1:
        xmin = 1
    else:
        xmin = int(numpy.floor(x-r))
    if numpy.ceil(x+r) > table.shape[1]:
        xmax = table.shape[1]
    else:
        xmax = int(numpy.ceil(x+r))
    if numpy.floor(y-r) < 1:
        ymin = 1
    else:
        ymin = int(numpy.floor(y-r))
    if numpy.ceil(y+r) > table.shape[0]:
        ymax = table.shape[0]
    else:
        ymax = int(numpy.ceil(y+r))
    #
    field = table[(ymin-1):ymax,(xmin-1):xmax]
    bckmin1 = bck[0]+2*bck[1]
    bckmin2 = minmax[1]
    if bckmin1 > minmax[1]:
        bckmin1 = minmax[1]
    #
    nrows, ncols = field.shape
    yf, xf = numpy.mgrid[:nrows, :ncols]
    #
    l = [(field >= bckmin1) & (field <= bckmin2)]
    #
    xc = ((field[tuple(l)]-bck[0])**2*xf[tuple(l)]).sum()
    yc = ((field[tuple(l)]-bck[0])**2*yf[tuple(l)]).sum()
    tot = ((field[tuple(l)]-bck[0])**2).sum()
    #
    if tot != 0.0:
            return xc/tot+xmin, yc/tot+ymin
    else:
        return x, y


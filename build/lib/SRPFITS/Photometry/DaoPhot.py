""" Utility functions and classes for SRP

Context : SRP
Module  : SRPFITS
Version : 1.1.0
Author  : Stefano Covino
Date    : 40/08/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   :

Remarks :

History : (22/08/2017) First version.
        : (30/08/2017) More sky algorithms.

"""


import numpy
from PythonPhot import aper


def DaoPhot(x,y,data,rds=(5,10,15),backgr=True,gain=1.,ron=0.,skyalg='mmm'):
    if backgr:
        mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = aper.aper(data,x,y,phpadu=gain,apr=rds[0],zeropoint=30,skyrad=[rds[1],rds[2]],readnoise=ron,skyalgorithm=skyalg)
    else:
        mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = aper.aper(data,x,y,phpadu=gain,apr=rds[0],zeropoint=30,skyrad=[rds[1],rds[2]],setskyval=0.0,readnoise=ron)
    #
    return flux.flatten(),fluxerr.flatten()
    



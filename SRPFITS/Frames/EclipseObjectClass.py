""" Utility functions and classes for SRP

Context : SRP
Module  : Frames.py
Version : 1.1.1
Author  : Stefano Covino
Date    : 16/05/2017
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (27/09/2010) First version.
        : (29/09/2010) More ordered importing.
        : (31/07/2015) python3 porting.
        : (16/05/2017) SRPFits porting.
"""

import os

from SRPFITS.Frames import EclipseConstants
from SRP.SRPSystem.Pipe import Pipe
from SRP.SRPSystem.Which import Which



# Peak data
class EclipseObjects:
    class Object:
        def __init__ (self, dlista):
            self.Id = dlista[0]
            self.X = float(dlista[1])
            self.Y = float(dlista[2])
            self.pix = float(dlista[3])
            self.mean = float(dlista[4])
            self.dev = float(dlista[5])
            self.med = float(dlista[6])
            self.min = float(dlista[7])
            self.max = float(dlista[8])
            self.fx = float(dlista[9])
            self.fy = float(dlista[10])
            self.FWHM = float(dlista[11])
            self.flux = float(dlista[12])
            self.RA = self.X
            self.DEC = self.Y

        def __str__ (self):
            ostr = "%5s\t%7.2f\t%7.2f\t%5d\t" % (self.Id, self.X, self.Y, self.pix)
            ostr = ostr + "%8.2f\t%8.2f\t" % (self.min, self.max)
            ostr = ostr + "%7.2f\t%7.2f\t%7.2f\t%10.3f\t" % (self.fx, self.fy, self.FWHM, self.flux)
            ostr = ostr + "%15.5f\t%15.5f" % (self.RA, self.DEC)
            return ostr+os.linesep

        def __lt__ (self, other):
            return self.flux > other.flux


    def __init__ (self, fitsfile, level=2.0):
        self.FitsFile = fitsfile
        self.Level = level
        self.ListEntries = []
        

    def FindEclipseObjects (self):
        # Eclipse required
        if Which(EclipseConstants.SRPpeak) == None:
            return None,EclipseConstants.EclipseNotFound
        # star list
        ListEntries = []
        argl = " -f '5 10 15' -F -P '5 10 15' -m clip -k %s " % self.Level
        cmd = EclipseConstants.SRPpeak+argl+self.FitsFile
        res = Pipe(cmd)
        if res == None:
            return None,EclipseConstants.EclipseNotExecutable
        # parse output
        for l in res.decode().split(os.linesep):
            if len(l) > 2:
                try:
                    ListEntries.append(self.Object(l.split()))
                except ValueError:
                    pass
        # 
        self.ListEntries = ListEntries
        #
        return len(self.ListEntries)

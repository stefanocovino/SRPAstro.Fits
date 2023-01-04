""" Utility functions and classes for SRP

Context : SRP
Module  : Frames
Version : 1.2.0
Author  : Stefano Covino
Date    : 26/05/2015
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (27/10/2014) First version.
        : (31/07/2015) python3 porting.
        : (16/05/2017) SRPFits porting.
        : (26/05/2017) Room for FWHM added.
"""

import os

from SRPFITS.Frames import SexConstants
#from SRP.SRPSystem.Pipe import Pipe
#from SRP.SRPSystem.Which import Which

from astropy.io.fits import getdata
import numpy
import sep


# Table data
class SexObjects:
    class Object:
        def __init__ (self, dlista):
            self.Id = dlista[0]
            self.X = float(dlista[1])
            self.Y = float(dlista[2])
            self.npix = float(dlista[3])
            self.ellip = 1-(float(dlista[4])/float(dlista[5]))
            self.flux = float(dlista[6])
            self.peak = float(dlista[7])
            self.flag = float(dlista[8])
            self.FWHM = -99.
            self.RA = self.X
            self.DEC = self.Y

        def __str__ (self):
            ostr = "%5s\t%7.2f\t%7.2f\t%5d\t" % (self.Id, self.X, self.Y, self.npix)
            ostr = ostr + "%5.2f\t%8.2f\t" % (self.ellip, self.flux)
            ostr = ostr + "%8.2f\t%5d\t" % (self.peak, self.flag)
            ostr = ostr + "%5.2f\t" % self.FWHM
            ostr = ostr + "%15.5f\t%15.5f" % (self.RA, self.DEC)
            return ostr+os.linesep

        def __lt__ (self, other):
            return self.flux > other.flux


    def __init__ (self, fitsfile, level=3.0):
        self.FitsFile = fitsfile
        self.Level = level
        self.ListEntries = []
        

    def FindSexObjects (self):
        # star list
        ListEntries = []
        #
        try:
            dt = getdata(self.FitsFile)
        except IOError:
            return None,SexConstants.SexFrameNotFound
        # correc binary format
        data = dt.astype(numpy.float32)
        # background evaluation and subtraction
        bkg = sep.Background(data)
        bkg.subfrom(data)
        # source extraction
        objs = sep.extract(data, self.Level * bkg.globalrms)
        # parse output
        for l in range(len(objs)):
            ListEntries.append(self.Object((l+1,objs['x'][l],objs['y'][l],objs['npix'][l],objs['b'][l],objs['a'][l],objs['flux'][l],objs['peak'][l],objs['flag'][l])))
        #
        self.ListEntries = ListEntries
        #
        return len(self.ListEntries)

""" Utility functions and classes for SRP

Context : SRP
Module  : Frames
Version : 1.2.0
Author  : Stefano Covino
Date    : 02/03/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (16/05/2017) First version.
        : (26/05/2017) Room for FWHM added.
		: (24/06/2020) Maxiters rather than iters in sigma_clipped_stats.
        : (02/03/2021) Sorted output.
"""

import os

from astropy.io.fits import getdata
from astropy.stats import sigma_clipped_stats
import numpy as np
from photutils import DAOStarFinder


# Table data
class DAOObjects:
    class Object:
        def __init__ (self, dlista):
            self.Id = dlista[0]
            self.X = float(dlista[1])
            self.Y = float(dlista[2])
            self.Sharpness = float(dlista[3])
            self.Roundness1 = float(dlista[4])
            self.Roundness2 = float(dlista[5])
            self.Npix = float(dlista[6])
            self.Sky = float(dlista[7])
            self.Peak = float(dlista[8])
            self.Flux = float(dlista[9])
            self.FWHM = -99
            self.RA = self.X
            self.DEC = self.Y

        def __str__ (self):
            ostr = "%5s\t%7.2f\t%7.2f\t%5.2f\t" % (self.Id, self.X, self.Y, self.Sharpness)
            ostr = ostr + "%5.2f\t%5.2f\t" % (self.Roundness1, self.Roundness2)
            ostr = ostr + "%5d\t%7.2f\t" % (self.Npix, self.Sky)
            ostr = ostr + "%7.2f\t%7.2g\t" % (self.Peak, self.Flux)
            ostr = ostr + "%5.2f\t" % (self.FWHM)
            ostr = ostr + "%15.5f\t%15.5f" % (self.RA, self.DEC)
            return ostr+os.linesep

        def __lt__ (self, other):
            return self.Flux > other.Flux


    def __init__ (self, fitsfile, level=3.0, fwhm=5):
        self.FitsFile = fitsfile
        self.Level = level
        self.FWHM = fwhm
        self.ListEntries = []
        

    def FindDAOObjects (self):
        # star list
        ListEntries = []
        #
        try:
            dt = getdata(self.FitsFile)
        except IOError:
            return None,-1
        #
        mean, median, std = sigma_clipped_stats(dt, sigma=3, maxiters=5)
        #
        daofind = DAOStarFinder(fwhm=self.FWHM, threshold=self.Level*std)
        objs = daofind(dt - median)
        #
        # parse output
        #objs.sort('flux')
        #objs.reverse()
        #
        for l in range(len(objs)):
            ListEntries.append(self.Object((l+1,objs['xcentroid'][l],objs['ycentroid'][l],
                objs['sharpness'][l],objs['roundness1'][l],objs['roundness2'][l],
                objs['npix'][l],objs['sky'][l],objs['peak'][l],objs['flux'][l])))
        #
        self.ListEntries = ListEntries
        #
        return len(self.ListEntries)


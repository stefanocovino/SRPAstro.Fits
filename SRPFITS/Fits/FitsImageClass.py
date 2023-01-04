""" Utility functions and classes for SRP

Context : SRP
Module  : Fits.py
Version : 1.9.7
Author  : Stefano Covino
Date    : 02/03/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (23/05/2010) First version.
        : (15/07/2010) WCS added.
        : (24/08/2010) Minor correction and source list sorting.
        : (25/08/2010) Frame size.
        : (30/08/2010) Fits Headers.
        : (04/09/2010) New importing rules.
        : (27/09/2010) Eclipse sources added.
        : (28/09/2010) Named changed from ObjList to FitsImage and adapted to native source class. 
        : (29/09/2010) Naming rationalized.
        : (04/10/2010) Possibility to remove stars close to frame border.
        : (25/10/2010) Basic statistics for frames.
        : (19/08/2013) Possiblity to read different extensions and original BITPIX saved.
        : (16/10/2013) Frame center in pixel.
        : (21/01/2014) Possibility to get statistixsw for just a subregion.
        : (27/10/2014) Sextractor source finding.
        : (16/05/2017) Minor update.
        : (26/05/2017) DAOPHOT source finding.
		: (30/11/2020) Minor GetWCS syntax correction.
        : (03/12/2020) More WCS corrections.
		: (11/12/2020) Again minor corrections.
        : (13/01/2020) Astropy WCS library.
		: (28/01/2021) Matt Hilton corrected the bugs of astLib.
        : (02/02/2021) Better management of astropy.log.
        : (02/03/2021) Better sorting.
"""

import os

import numpy

from .GetData import GetData
from .GetHeader import GetHeader
from SRPFITS.Fits.GetWCS import GetWCS
from .GetHeaderValue import GetHeaderValue
from . import FitsConstants as FitsConstants
from SRPFITS.GetFWHM import GetFWHM

from SRPFITS.Frames.SourceObjectsClass import SourceObjects
from SRPFITS.Frames.SexObjectClass import SexObjects
from SRPFITS.Frames.DAOObjectClass import DAOObjects
from SRPFITS.Frames.Pixel2WCS import Pixel2WCS
from astropy import log
log.setLevel('WARNING')

class FitsImage:
    def __init__ (self, fitsfile, extension=0):
        self.Name = fitsfile
        self.Header = GetHeader(fitsfile,extension)[0]
        self.Data = GetData(fitsfile,extension)[0]
        self.WCS = GetWCS(fitsfile,'astlib')[0]
        self.BITPIX = GetHeaderValue(fitsfile,'BITPIX',extension)[0]
        self.NAXIS = GetHeaderValue(fitsfile,'NAXIS',extension)[0]
        self.List = []
        self.NativeSourcesFlag = False
        self.DAOSourcesFlag = False
        self.SexSourcesFlag = False


    def Sources(self, threshold=5.0, filtsing=3):
        slist = SourceObjects(self.Name)
        slist.FindObjects(self.Data, threshold, filtsing)
        srclist = []
        for i in slist.ListEntries:
            srclist.append((i.X,i.Y))
        coolist = Pixel2WCS(self.Header,srclist,'astlib')
        for i,l in zip(slist.ListEntries,coolist):
            i.RA = l[0]
            i.DEC = l[1]
        self.List = slist.ListEntries
        self.NativeSourcesFlag = True
        self.DAOSourcesFlag = False
        self.SexSourcesFlag = False
        return len(self.List)
        
        
    def SortSourceList(self):
        if self.NativeSourcesFlag:
            self.List.sort()
            #self.List.reverse()
        elif self.DAOSourcesFlag:
            self.List.sort()
            #self.List.reverse()
        elif self.SexSourcesFlag:
            self.List.sort()
            #self.List.reverse()


    def GetFWHM (self):
        X = []
        Y = []
        for i in self.List:
            X.append(i.X)
            Y.append(i.Y)
        fwhml = GetFWHM(X,Y,self.Name)
        for i,l in zip(self.List,fwhml):
            i.FWHM = l
        return len(self.List)

        
    def GetFrameSizePix(self):
        base = FitsConstants.NAxis
        sizelist = []
        try:
            naxis = self.Header[base]
        except KeyError:
            return None
        for i in range(naxis):
            baseax = '%s%s' % (base,i+1)
            try:
                value = self.Header[baseax]
            except KeyError:
                value = None
            sizelist.append(value)
        return sizelist



    def SexSources(self,level=3.0):
        elist = SexObjects(self.Name,level)
        elist.FindSexObjects()
        srclist = []
        for i in elist.ListEntries:
            srclist.append((i.X,i.Y))
        try:
            coolist = Pixel2WCS(self.Header,srclist,'astlib' )
        except AttributeError:
            coolist = []
            for ii in range(len(srclist)):
                coolist.append((0.,0.))
        for i,l in zip(elist.ListEntries,coolist):
            i.RA = l[0]
            i.DEC = l[1]
        self.List = elist.ListEntries
        self.SexSourcesFlag = True
        self.DAOSourcesFlag = False
        self.NativeSourcesFlag = False
        return len(self.List)



    def DAOSources(self,level=3.0):
        elist = DAOObjects(self.Name,level)
        elist.FindDAOObjects()
        srclist = []
        for i in elist.ListEntries:
            srclist.append((i.X,i.Y))
        try:
            coolist = Pixel2WCS(self.Header,srclist,'astlib')
        except AttributeError:
            coolist = []
            for ii in range(len(srclist)):
                coolist.append((0.,0.))
        for i,l in zip(elist.ListEntries,coolist):
            i.RA = l[0]
            i.DEC = l[1]
        self.List = elist.ListEntries
        self.SexSourcesFlag = False
        self.DAOSourcesFlag = True
        self.NativeSourcesFlag = False
        return len(self.List)



    def CleanBorderSources(self,cleaningpercentage=1):
        framesize = self.GetFrameSizePix()
        avoidzone = (framesize[0]*cleaningpercentage/100.,framesize[1]*cleaningpercentage/100.)
        CleanList = []
        oldnobj = len(self.List)
        for entry in self.List:
            if (1+avoidzone[0] <= entry.X <= framesize[0]-avoidzone[0]) and (1+avoidzone[1] <= entry.Y <= framesize[1]-avoidzone[1]):
                CleanList.append(entry)
        self.List = CleanList
        return len(self.List),oldnobj
        
        
    def __str__(self):
        msg = ''
        # Native
        if self.NativeSourcesFlag:
            for i in range(len(self.List)):
                msg = msg + str(self.List[i])
        # daophot
        elif self.DAOSourcesFlag:
            for i in range(len(self.List)):
                msg = msg + str(self.List[i])
        # Sex
        elif self.SexSourcesFlag:
            for i in range(len(self.List)):
                msg = msg + str(self.List[i])
        return msg
    

    def Skycat(self, outname='SRP.cat'):
        msg = ''
        msg = msg + "long_name: SRP catalog for file %s\n" % (self.Name)
        msg = msg + "short_name: %s\n" % (outname)
        msg = msg + "url: ./%s\n" % (outname)
        msg = msg + "symbol: {} {circle blue} 4\n"
        msg = msg + "id_col: 0\n"
        msg = msg + "x_col: 1\n"
        msg = msg + "y_col: 2\n"
        # Native
        if self.NativeSourcesFlag:
            msg = msg + "Id\tX\tY\tRA\tDEC\tNPix\tMag\tFWHM\n"
        # daophot
        elif self.DAOSourcesFlag:
            msg = msg + "Id\tX\tY\tSharpness\tRoundness1\tRoundness2\tNpix\tSky\tPeak\tFlux\tFWHM\tRA\tDEC\n"
        # Sex
        elif self.SexSourcesFlag:
            msg = msg + "Id\tX\tY\tNpix\tEllipticity\tFlux\tPeak\tFlag\tFWHM\tRA\tDEC\n"
        
        msg = msg + "---------\n"
        msg = msg + str(self)
        msg = msg + "EOD\n"
        return msg
        
        
    def GetStats(self,region=None):
        """
        region is a list of four pixel position (lower left is 1,1)
        They are leftx, bottomy, rightx, uppery.
        """
        if region == None or len(region) != 4:
            mean = numpy.mean(self.Data)
            std = numpy.std(self.Data)
            median = numpy.median(self.Data)
            max = numpy.max(self.Data)
        else:
            try:
                mean = numpy.mean(self.Data[region[1]+1:region[3]+1,region[0]+1:region[2]+1])
                std = numpy.std(self.Data[region[1]+1:region[3]+1,region[0]+1:region[2]+1])
                median = numpy.median(self.Data[region[1]+1:region[3]+1,region[0]+1:region[2]+1])
                max = numpy.max(self.Data[region[1]+1:region[3]+1,region[0]+1:region[2]+1])
            except ValueError:
                return numpy.nan,numpy.nan,numpy.nan,numpy.nan
        #
        return mean,std,median,max



    def GetFrameCenterPix (self):
        fs = self.GetFrameSizePix()
        return int(fs[0]/2.), int(fs[1]/2.)


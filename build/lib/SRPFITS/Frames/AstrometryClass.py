""" Utility functions and classes for SRP

Context : SRP
Module  : Frames.py
Version : 1.13.0
Author  : Stefano Covino
Date    : 17/03/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (27/09/2010) First version.
        : (28/09/2010) Minor name changes.
        : (29/09/2010) Minor correction.
        : (01/10/2010) Warning in saving FITS file suppressed.
        : (03/10/2010) Number of stars for astrometry. RA and DEC shift computed.
        : (04/10/2010) Sources close to frame border cleaned.
        : (13/10/2010) Larger search radius for catalogue object and no more maximum 
        :               number of objects.
        : (14/10/2010) Better estimate of frame size.
        : (24/10/2010) USNO-A2 catalogue also be used.
        : (25/10/2010) Flexible threshold for bad quality frames.
        : (27/10/2010) Simpler comment output.
        : (03/11/2010) Better comment formatting.
        : (17/11/2010) WCS header comments.
        : (20/11/2010) Better header formatting.
        : (31/12/2010) Minor bug correction.
        : (26/08/2011) Consistency check for unreliable CRVAL1,2 and CDELT1,2.
        : (04/11/2011) More flexibility in header coordinate format.
        : (22/01/2012) Check that at least one catalogue returns data.
        : (19/08/2013) Data saved with the original format.
        : (29/08/2013) Tunable tolerance for triangle match. Do not update RA,DEC if they already exist.
        : (03/09/2013) Better search of catalogue objects.
        : (16/10/2013) Better computation of frame center.
        : (16/07/2014) Deal with non standard FITS headers.
        : (30/10/2014) Sextractor source finding.
        : (09/12/2014) Better choice of threshold for sextractor source finding.
        : (12/01/2015) No more del in dictionaries.
        : (31/07/2015) python3 porting.
        : (14/03/2017) Less iterations in astrometry if rotation and/or pizel size is provided.
        : (16/05/2017) Minor update.
        : (18/05/2017) Minor update.
        : (31/10/2018) 'astlib' in computation.
        : (03/12/2020) Minor bugs in calling 'astlib'.
        : (13/01/2021) Astropy WCS adopted.
        : (28/01/2021) Matt Hilton corrected the bugs in astLib.
        : (09/02/2021) DAOPHOT added astroalign algorithm.
        : (17/03/2021) astroalign integrated a solution hunter.
"""


import copy, math, warnings

import astropy.io.fits as aif
import astLib.astCoords as aLaC
import numpy
import scipy.optimize as so
import astroalign as aa
import SRP
from SRP.SRPMath.AstroCoordInput import AstroCoordInput
from SRP.SRPCatalogue.TWOMASSClass import TWOMASS
from SRP.SRPCatalogue.USNOClass import USNO
from SRPFITS.Fits.FitsImageClass import FitsImage
from SRPFITS.Fits import FitsConstants
from SRP.SRPMath.AngularDistance import AngularDistance
from SRP.SRPMath.AstroAngleInput import AstroAngleInput
from SRP.SRPMath.AstroCoordInput import AstroCoordInput
from SRP.SRPMath.TriangleClass import Triangle
from SRPFITS.Math.TriangleMatch import TriangleMatch
from SRPSTATS.CoordDistanceMinimization import CoordDistanceMinimization
from SRPSTATS.SkyDistanceMinimization import SkyDistanceMinimization
from SRPSTATS.SkyDistanceMinimization import SkyDistSum
from SRPFITS.Frames.WCS2Pixel import WCS2Pixel
from SRPFITS.Frames.Pixel2WCS import Pixel2WCS
from SRPFITS.Fits.WCSRotationDeg import WCSRotationDeg


class Astrometry:
    def __init__ (self, fitsfile, center=None, pixcenter=None, point=None, pixsize=None, rotangle=None, framesize=None, maxres=3.0):
            # Fits file
            self.FitsFrame = FitsImage(fitsfile)
            # WCS
            # Pointing direction
            if point != None:
                pcoord = AstroCoordInput(point[0],point[1])
                self.RA = pcoord.RA
                self.DEC = pcoord.DEC
            else:
                try:
                    heacoord = AstroCoordInput(self.FitsFrame.Header[FitsConstants.RA],self.FitsFrame.Header[FitsConstants.DEC])
                except KeyError:
                    heacoord = AstroCoordInput(0.,0.)
                self.RA = heacoord.RA
                self.DEC = heacoord.DEC
            # Pixel reference
            if pixcenter != None:
                self.CRPIX1 = pixcenter[0]
                self.CRPIX2 = pixcenter[1]
            else:
                try:
                    self.CRPIX1 = self.FitsFrame.Header[FitsConstants.CRPIX1]
                    self.CRPIX2 = self.FitsFrame.Header[FitsConstants.CRPIX2]
                except KeyError:
                    self.CRPIX1 = 0.
                    self.CRPIX2 = 0.
            # Pixel size
            if pixsize != None:
                self.CDELT1 = pixsize[0]
                self.CDELT2 = pixsize[1]
                self.pixsizeprov = True
            else:
                self.pixsizeprov = False
                try:
                    self.CDELT1 = self.FitsFrame.Header[FitsConstants.CDELT1]
                    self.CDELT2 = self.FitsFrame.Header[FitsConstants.CDELT2]                    
                except KeyError:
                    self.CDELT1 = -1./3600.
                    self.CDELT2 = 1./3600.
            # check consistency of pixel size, it should not be larger than 0.1deg
            if math.fabs(self.CDELT1) > 0.1:
                self.CDELT1 = -1./3600.
            if math.fabs(self.CDELT2) > 0.1:
                self.CDELT2 = 1./3600.
            # Rotation angle
            if rotangle != None:
                angrad = math.radians(rotangle)
                self.RotationAngle = rotangle
                self.PC11 = math.cos(angrad)
                self.PC12 = -math.sin(angrad)
                self.PC21 = math.sin(angrad)
                self.PC22 = math.cos(angrad)
                self.pc = True
                self.rotangleprov = True
            else:
                self.rotangleprov = False
                try:
                    self.RotationAngle = self.FitsFrame.WCS.getRotationDeg()
                except AttributeError:
                    self.RotationAngle = 0.0
                try:
                    self.PC11 = self.FitsFrame.Header[FitsConstants.PC11]
                    self.PC12 = self.FitsFrame.Header[FitsConstants.PC12]
                    self.PC21 = self.FitsFrame.Header[FitsConstants.PC21]
                    self.PC22 = self.FitsFrame.Header[FitsConstants.PC22]
                    self.pc = True
                except KeyError:
                    self.PC11 = 1.
                    self.PC12 = 0.
                    self.PC21 = 0.
                    self.PC22 = 1.
                    self.pc = False
                try:
                    self.CD11 = self.FitsFrame.Header[FitsConstants.CD11]
                    self.CD12 = self.FitsFrame.Header[FitsConstants.CD12]
                    self.CD21 = self.FitsFrame.Header[FitsConstants.CD21]
                    self.CD22 = self.FitsFrame.Header[FitsConstants.CD22]
                    self.cd = True
                except KeyError:
                    self.cd = False
            if self.pc == False and self.cd == True:
                self.PC11 = self.CD11/self.CDELT1
                self.PC12 = self.CD12/self.CDELT2
                self.PC21 = self.CD21/self.CDELT1
                self.PC22 = self.CD22/self.CDELT2
            # Coordinate reference
            if center != None:
                coord = AstroCoordInput(center[0],center[1])
                self.CRVAL1 = coord.RA
                self.CRVAL2 = coord.DEC
            else:
                try:
                    heacoord = AstroCoordInput(self.FitsFrame.Header[FitsConstants.CRVAL1],self.FitsFrame.Header[FitsConstants.CRVAL2])
                    self.CRVAL1 = heacoord.RA
                    self.CRVAL2 = heacoord.DEC
                except KeyError:
                    self.CRVAL1 = self.RA
                    self.CRVAL2 = self.DEC
            # Check consistency between RA,DEC and CRVAL1,2
            if AngularDistance((self.RA,self.DEC),(self.CRVAL1,self.CRVAL2)) > 1.5:
                self.CRVAL1 = self.RA
                self.CRVAL2 = self.DEC
            # Image size
            sizelistpix = self.FitsFrame.GetFrameSizePix()
            if framesize != None:
                self.HalfSize = framesize/2.
            else:
                pxs = sum([math.fabs(i) for i in (self.CDELT1,self.CDELT2)])/2.0
                self.HalfSize = sum([i*pxs/2.0 for i in sizelistpix])/2.0
            #
            # Axes type
            try:
                self.CTYPE1 = self.FitsFrame.Header[FitsConstants.CTYPE1]
                self.CTYPE2 = self.FitsFrame.Header[FitsConstants.CTYPE2]
            except KeyError:
                self.CTYPE1 = FitsConstants.RAType
                self.CTYPE2 = FitsConstants.DECType
            # Strange CTYPES
            if self.CTYPE1 != FitsConstants.RAType or self.CTYPE1 != FitsConstants.DECType and self.CTYPE2 != FitsConstants.RAType or self.CTYPE2 != FitsConstants.DECType:
                self.CTYPE1 = FitsConstants.RAType
                self.CTYPE2 = FitsConstants.DECType   
            # Frame center
            self.FrameCenter = [(i/2.)+1.0 for i in sizelistpix]
            #
            self.List = []
            self.CatList = []
            self.Astrometrized = False
            self.Residuals = -1
            #
            self.MaxRes = maxres
            self.FitAstro = False
            self.TriAstro = False
            # Algorithms
            self.NativeSources = False
            self.DaoSources = False
            self.SexSources = False
            #
            self.DefineWCSHeader()
  
  
    def GetFrameCenterCoords(self):
            cc = self.FitsFrame.GetFrameCenterPix()
            cpos = Pixel2WCS(self.FitsFrame.Header,[cc],'astlib')
            return cpos[0]
  
  
    def GetSources(self, maxobjs=15, thre=23, mincnnt=3, cleanperc=1.):
            soglia = copy.copy(thre)
            while len(self.List) < maxobjs and soglia > 1: 
                self.FitsFrame.Sources(soglia,mincnnt)
                self.FitsFrame.CleanBorderSources(cleanperc)
                self.FitsFrame.SortSourceList()
                if len(self.FitsFrame.List) >= maxobjs:
                    self.List = self.FitsFrame.List[:maxobjs]
                else:
                    if soglia-5 > 3:
                        soglia = soglia - 5
                    elif soglia > 3:
                        soglia = 3
                    elif 2 < soglia <= 3:
                        soglia = soglia - 0.5
                    else:
                        soglia = soglia - 0.25
                    self.List = self.FitsFrame.List
            self.NativeSources = True
            self.DaoSources = False
            self.SexSources = False


    def GetDaoSources(self, maxobjs=15, thre=50, cleanperc=1.):
            soglia = copy.copy(thre)
            while len(self.List) < maxobjs and soglia > 0.1:
                self.FitsFrame.DAOSources(soglia)
                self.FitsFrame.CleanBorderSources(cleanperc)
                self.FitsFrame.SortSourceList()
                if len(self.FitsFrame.List) >= maxobjs:
                    self.List = self.FitsFrame.List[:maxobjs]
                else:
                    if soglia >= 10:
                        soglia = soglia - 5
                    elif 3 <= soglia < 10:
                        soglia = soglia - 1
                    elif 1 <= soglia < 3:
                        soglia = soglia - 0.5
                    else:
                        soglia = soglia - 0.1
                    self.List = self.FitsFrame.List
            self.DaoSources = True
            self.NativeSources = False
            self.SexSources = False



    def GetSexSources(self, maxobjs=15, thre=100, cleanperc=1.):
            soglia = copy.copy(thre)
            while len(self.List) < maxobjs and soglia > 0.1:
                self.FitsFrame.SexSources(soglia)
                self.FitsFrame.CleanBorderSources(cleanperc)
                self.FitsFrame.SortSourceList()
                if len(self.FitsFrame.List) >= maxobjs:
                    self.List = self.FitsFrame.List[:maxobjs]
                else:
                    if soglia >= 10:
                        soglia = soglia - 5
                    elif 3 <= soglia < 10:
                        soglia = soglia - 1
                    elif 1 <= soglia < 3:
                        soglia = soglia - 0.5
                    else:
                        soglia = soglia - 0.1
                    self.List = self.FitsFrame.List
            self.SexSources = True
            self.NativeSources = False
            self.DaoSources = False


    def GetCatSources(self, maxentr=15, catq='N'): 
            radius = 1.5*self.HalfSize*60.0  # arcmin
            #cntfld = (self.RA,self.DEC)
            crd = self.GetFrameCenterCoords()
            inpcoord = AstroCoordInput(crd[0],crd[1])
            cntfld = (inpcoord.RA, inpcoord.DEC)
            #print cntfld
            #
            if catq == 'N':
                cat = TWOMASS(cntfld[0], cntfld[1], math.floor(radius))
            elif catq == 'O':
                cat = USNO(cntfld[0], cntfld[1], math.floor(radius))
            else:
                cat = USNO(cntfld[0], cntfld[1], math.floor(radius))
            #
            if cat != None:
                cat.GetData()
                cat.sort()
                if len(cat.ListEntries) >= maxentr:
                    self.CatList = cat.ListEntries[:maxentr]
                else:
                    self.CatList = cat.ListEntries    


    def DefineWCSHeader(self):
        # remove alternate WCS
        self.FitsFrame.Header.pop(FitsConstants.CD11,None)
        self.FitsFrame.Header.pop(FitsConstants.CD12,None)
        self.FitsFrame.Header.pop(FitsConstants.CD21,None)
        self.FitsFrame.Header.pop(FitsConstants.CD22,None)
        # center ref
        self.FitsFrame.Header.set(FitsConstants.CRVAL1,self.CRVAL1,FitsConstants.CRVALCmt)
        self.FitsFrame.Header.set(FitsConstants.CRVAL2,self.CRVAL2,FitsConstants.CRVALCmt)
        #             
        self.FitsFrame.Header.set(FitsConstants.CTYPE1,self.CTYPE1,FitsConstants.CTYPECmt)
        self.FitsFrame.Header.set(FitsConstants.CTYPE2,self.CTYPE2,FitsConstants.CTYPECmt)
        # rotation
        self.FitsFrame.Header.set(FitsConstants.PC11,self.PC11,FitsConstants.PCCmt)
        self.FitsFrame.Header.set(FitsConstants.PC12,self.PC12,FitsConstants.PCCmt)
        self.FitsFrame.Header.set(FitsConstants.PC21,self.PC21,FitsConstants.PCCmt)
        self.FitsFrame.Header.set(FitsConstants.PC22,self.PC22,FitsConstants.PCCmt)
        # pixel size
        self.FitsFrame.Header.set(FitsConstants.CDELT1,self.CDELT1,FitsConstants.CDELT1Cmt)
        self.FitsFrame.Header.set(FitsConstants.CDELT2,self.CDELT2,FitsConstants.CDELT2Cmt)
        # pixel ref
        self.FitsFrame.Header.set(FitsConstants.CRPIX1,self.CRPIX1,FitsConstants.CRPIX1Cmt)
        self.FitsFrame.Header.set(FitsConstants.CRPIX2,self.CRPIX2,FitsConstants.CRPIX1Cmt)
        # pointing dir
        if not (FitsConstants.RA in list(self.FitsFrame.Header.keys())):
            self.FitsFrame.Header.set(FitsConstants.RA,self.RA)
        if not (FitsConstants.DEC in list(self.FitsFrame.Header.keys())):
            self.FitsFrame.Header.set(FitsConstants.DEC,self.DEC)
        self.FitsFrame.Header.set(FitsConstants.RADECSys,FitsConstants.RADECSysVal)
        self.FitsFrame.Header.set(FitsConstants.AngUnit1,FitsConstants.AngUnit1Val,FitsConstants.AngUnCmt)
        self.FitsFrame.Header.set(FitsConstants.AngUnit2,FitsConstants.AngUnit2Val,FitsConstants.AngUnCmt)
        #
        self.FitsFrame.WCS.updateFromHeader()


    def UpgradeWCSHeader(self, xxx_todo_changeme):
        (crpix1,crpix2,rot,sf) = xxx_todo_changeme
        rotrnd = math.radians(rot)
        self.FitsFrame.Header.set(FitsConstants.PC11,math.cos(rotrnd))
        self.FitsFrame.Header.set(FitsConstants.PC12,-math.sin(rotrnd))
        self.FitsFrame.Header.set(FitsConstants.PC21,math.sin(rotrnd))
        self.FitsFrame.Header.set(FitsConstants.PC22,math.cos(rotrnd))
        if self.CDELT1 >= 0:
            cdelt1 = sf
        else:
            cdelt1 = -sf
        if self.CDELT2 >= 0:
            cdelt2 = sf
        else:
            cdelt2 = -sf
        self.FitsFrame.Header.set(FitsConstants.CDELT1,cdelt1)
        self.FitsFrame.Header.set(FitsConstants.CDELT2,cdelt2)
        self.FitsFrame.Header.set(FitsConstants.CRPIX1,crpix1)
        self.FitsFrame.Header.set(FitsConstants.CRPIX2,crpix2)
        self.FitsFrame.WCS.updateFromHeader()            
  

    
    def ImproveWCS (self, xxx_todo_changeme1):
        (sx,sy,rot,sf) = xxx_todo_changeme1
        self.CDELT1 = self.CDELT1 * sf
        self.CDELT2 = self.CDELT2 * sf
        self.CRPIX1 = self.CRPIX1 + sx
        self.CRPIX2 = self.CRPIX2 + sy
        self.RotationAngle = self.RotationAngle + rot
          
          
                                    
    def FindSolution(self,angtol=0.5,disttol=5.):
        self.GetAAHints()
        #
        reflist = []
        objlist = []
        # centering
        if self.NativeSources:
            for i in self.List:
                objlist.append(((i[0]-self.CRPIX1),(i[1]-self.CRPIX2)))
                #print objlist[-1]
        elif self.DaoSources:
            for i in self.List:
                objlist.append(((i.X-self.CRPIX1),(i.Y-self.CRPIX2)))
                #print objlist[-1]
                #print i.flux
        elif self.SexSources:
            for i in self.List:
                objlist.append(((i.X-self.CRPIX1),(i.Y-self.CRPIX2)))
                #print objlist[-1]
                #print i.flux
        #print 
        radec = []
        for i in self.CatList:
            radec.append((i.RA,i.DEC))
            #print i.H
        RefPos = WCS2Pixel(self.FitsFrame.Header,[(self.CRVAL1,self.CRVAL2)],'astlib')
        #print RefPos
        CatListPix = WCS2Pixel(self.FitsFrame.Header,radec,'astlib')
        for i in CatListPix:
            reflist.append(((i[0]-RefPos[0][0]), (i[1]-RefPos[0][1])))
            #print reflist[-1]
        #
        self.REFCDR = reflist
        self.OBJCDR = objlist
        if self.pixsizeprov and self.rotangleprov:
            mtchres = TriangleMatch(reflist,objlist,angtol,disttol,TrueSize=abs(self.CDELT1)*3600.,TrueAng=self.RotationAngle)
        elif self.pixsizeprov and not self.rotangleprov:
            mtchres = TriangleMatch(reflist,objlist,angtol,disttol,TrueSize=abs(self.CDELT1)*3600.)
        elif not self.pixsizeprov and self.rotangleprov:
            mtchres = TriangleMatch(reflist,objlist,angtol,disttol,TrueAng=self.RotationAngle)
        else:
            mtchres = TriangleMatch(reflist,objlist,angtol,disttol)
        if mtchres != None:
#            print mtchres
            mtchinplist = []
            mtchoutlist = []
            for i in mtchres[0]:
                if self.NativeSources:
                    mtchinplist.append((self.List[i][0],self.List[i][1]))
                elif self.DaoSources:
                    mtchinplist.append((self.List[i].X,self.List[i].Y))
                elif self.SexSources:
                    mtchinplist.append((self.List[i].X,self.List[i].Y))
            for i in mtchres[1]:
                mtchoutlist.append((self.CatList[i].RA,self.CatList[i].DEC))
            self.NStarRes = len(mtchres[0])
            self.ImproveWCS(mtchres[2][0])
            ris = SkyDistSum((self.CRPIX1,self.CRPIX2,self.RotationAngle,math.fabs(self.CDELT1)),mtchinplist,mtchoutlist,self)/len(mtchoutlist)
            dataris = (self.CRPIX1,self.CRPIX2,self.RotationAngle,math.fabs(self.CDELT1))
#            print self.CRVAL1,self.CRVAL2,self.CRPIX1,self.CRPIX2,self.PC11,self.PC12,self.PC21,self.PC22,self.CDELT1,self.CDELT2
            res = SkyDistanceMinimization((self.CRPIX1,self.CRPIX2,self.RotationAngle,math.fabs(self.CDELT1)),mtchinplist,mtchoutlist,self)
            #print (res)
            self.Residuals = res[1]
            self.UpgradeWCSHeader(res[0])
#            print self.CRVAL1,self.CRVAL2,self.CRPIX1,self.CRPIX2,self.PC11,self.PC12,self.PC21,self.PC22,self.CDELT1,self.CDELT2
#            print ris,res[1]
            if self.Residuals <= self.MaxRes:
                self.Astrometrized = True
                self.FitAstro = True
            else:
                # try with no fit
                if ris < self.MaxRes:
                    self.Residuals = ris
                    self.UpgradeWCSHeader(dataris)
                    self.Astrometrized = True
                    self.TriAstro = True 
                else:
                    self.Astrometrized = False
            #
            self.GetPointingShift()
            #
        else:
            self.Astrometrized = False
            self.Residuals = -1
            self.NStarRes = 0
            


    def GetAAHints(self):
        reflist = []
        objlist = []
        #
        if self.NativeSources:
            for i in self.List:
                objlist.append((i[0],i[1]))
                #print objlist[-1]
        elif self.DaoSources:
            for i in self.List:
                objlist.append((i.X,i.Y))
                #print objlist[-1]
                #print i.flux
        elif self.SexSources:
            for i in self.List:
                objlist.append((i.X,i.Y))
                #print objlist[-1]
                #print i.flux
        #
        for i in self.CatList:
            reflist.append((i.RA,i.DEC))
        #
        self.REFCDR = reflist
        self.OBJCDR = objlist
        #
        mtchres = True
        try:
            transf, (source_list, target_list) = aa.find_transform(objlist, reflist)
        except aa.MaxIterError:
            mtchres = False
        #
        if mtchres:
            self.CRVAL1 = transf.translation[0]
            self.CRVAL2 = transf.translation[1]
            self.RotationAngle = transf.rotation
            self.CRPIXEL1 = 0
            self.CRPIXEL2 = 0
            self.UpgradeWCSHeader((0,0,transf.rotation,transf.scale))




    def GetPointingShift(self):
        # derive true center coordinate
        ccord = Pixel2WCS(self.FitsFrame.Header,[(self.FrameCenter[0],self.FrameCenter[1])],'astlib')
#        self.RAshift = ((self.RA - ccord[0][0])*math.cos(math.radians(self.DEC)))*3600.00
#        self.DECshift = (self.DEC - ccord[0][1])*3600.0
        self.RAshift = ((self.FitsFrame.Header[FitsConstants.RA] - ccord[0][0])*math.cos(math.radians(self.FitsFrame.Header[FitsConstants.DEC])))*3600.0
        self.DECshift = (self.FitsFrame.Header[FitsConstants.DEC] - ccord[0][1])*3600.0

            
    def SaveFile(self,nfile):
        frm = aif.open(self.FitsFrame.Name,do_not_scale_image_data=True)
        if self.Astrometrized:
            self.FitsFrame.Header.add_comment('SRPComment: Astrometric solution acceptable.') 
            self.FitsFrame.Header.add_comment('SRPComment: Average residual: %.2g arcsec for %d stars.' % (self.Residuals, self.NStarRes))  
            self.FitsFrame.Header.add_comment('SRPComment: Shift wrt pointing RA and DEC coords: %.1f %.1f arcsec.' % (self.RAshift, self.DECshift)) 
        else:
            self.FitsFrame.Header.add_comment('SRPComment: Astrometric solution not acceptable.') 
            self.FitsFrame.Header.add_comment('Average residual: %.2g arcsec for %d stars.' % (self.Residuals, self.NStarRes))              
        frm[0].header = self.FitsFrame.Header
        #warnings.resetwarnings()
        #warnings.filterwarnings('ignore', category=UserWarning, append=True)
        #warnings.filterwarnings('ignore', category=ResourceWarning, append=True)
        #if self.FitsFrame.BITPIX == 8:
        #    frm[0].scale('uint8')
        #elif self.FitsFrame.BITPIX == 16:
        #    frm[0].scale('uint16')
        #elif self.FitsFrame.BITPIX == 32:
        #    frm[0].scale('uint32')
        #elif self.FitsFrame.BITPIX == -32:
        #    frm[0].scale('float32')
        #elif self.FitsFrame.BITPIX == -64:
        #    frm[0].scale('float64')
        frm.writeto(nfile,overwrite=True,output_verify='ignore')
        #warnings.resetwarnings()
        #warnings.filterwarnings('always', category=UserWarning, append=True)
        #warnings.filterwarnings('always', category=ResourceWarning, append=True)




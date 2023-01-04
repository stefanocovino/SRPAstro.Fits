""" Utility functions and classes for SRP

Context : SRP
Module  : Math.py
Version : 1.1.3
Author  : Stefano Covino
Date    : 06/07/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (27/09/2010) First version.
        : (13/10/2010) Quicker triangle scan.
        : (26/08/2013) More control on accepted matches.
        : (14/03/2017) Possibility to remove not useful iterations 
            and porting to SRP.FITS.
        : (06/07/2021) SREPSTATS.
"""

import math

from SRP.SRPMath.TriangleClass import Triangle
from SRP.SRPMath.CartesianCoordAmpRotoTranslation import CartesianCoordAmpRotoTranslation
from SRPSTATS.CoordDistanceMinimization import CoordDistanceMinimization
from SRP.SRPMath.AngleRange import AngleRange



def TriangleMatch (reflist,objlist,angtol=0.5,disttol=5.,TrueSize=None,TrueAng=None):
#    print objlist
#    print reflist
    if len(reflist) >=3 and len(objlist) >=3:
        #
        if len(reflist) >= 5:
            reqmtch = 5
        else:
            reqmtch = len(reflist)
        #
        for i in range(len(reflist)):
            for ii in range(i,len(reflist)):
                for iii in range(ii,len(reflist)):
                    for l in range(len(objlist)):
                        for ll in range(l,len(objlist)):
                            for lll in range(ll,len(objlist)):
                                if i != ii and ii != iii and i != iii and l != ll and ll != lll and l != lll:
                                    #print i,ii,iii,l,ll,lll
                                    triref = Triangle(reflist[i],reflist[ii],reflist[iii],angtol)
                                    triobj = Triangle(objlist[l],objlist[ll],objlist[lll],angtol)
                                    if not (triobj.IsSymmTriangle() or triref.IsSymmTriangle()):
                                        # cannot work with symmetric triangles
                                        # Select not too small or elongated triangles
                                        if min(triobj.Sizes) >= 10.0 and min(triobj.Angles) > 10.0 and triobj.Commensurabili(triref):
                                            if triobj.Rotable(triref):
                                                # with a good candidate compute pairing data
                                                sf = triobj.SizeFactor(triref)
                                                rf = triobj.RotationFactor(triref)
                                                sh = triobj.ShiftFactor(triref)
                                                if TrueSize != None:
                                                    if TrueSize/sf > 1.1 or TrueSize/sf < 0.9:
                                                        continue
                                                if TrueAng != None:
                                                    if abs(AngleRange(TrueAng)-AngleRange(rf)) > 10.:
                                                        continue
                                                #print i,ii,iii,l,ll,lll
                                                #print ("sf,rf,sz", sf,rf,sh)
                                                #print "objs",triobj.Vertici[0],triobj.Vertici[1],triobj.Vertici[2]
                                                #print "refs",triref.Vertici[0],triref.Vertici[1],triref.Vertici[2]
                                                #print "obj-ang", triobj.Angles
                                                #print "ref-ag", triref.Angles
                                                #print "obj-siz", triobj.Sizes
                                                #print "ref-siz", triref.Sizes
                                                #print "ccns",triobj.ConnVert
                                                # find acceptable error
                                                tempinp = []
                                                tempout = []
                                                for zz in triobj.ConnVert:
                                                    tempinp.append(triobj.Vertici[zz[0]].Coord)
                                                    tempout.append(triref.Vertici[zz[1]].Coord)
                                                res = CoordDistanceMinimization(tempout,tempinp,[sh[0],sh[1],rf,1./sf])   
                                                # acceptable tolerance for match 
                                                minassdist = disttol*res[1]
                                                #print minassdist, res
                                                #print
                                                #    
                                                #minassdist = 0
                                                #for zz in triobj.ConnVert:
                                                #    miecord = CartesianCoordAmpRotoTranslation(triobj.Vertici[zz[0]].Coord,sh,rf,1./sf)
                                                #    miecorddist = math.sqrt((miecord[0]-triref.Vertici[zz[1]].Coord[0])**2 + (miecord[1]-triref.Vertici[zz[1]].Coord[1])**2)
                                                #    minassdist = minassdist + miecorddist 
                                                #    #print miecorddist,triobj.Vertici[zz[0]].Coord, miecord, triref.Vertici[zz[1]].Coord
                                                #minassdist = minassdist*5./3.
                                                # with a basic solution try to find a general one
                                                # pair points
                                                nowinp = []
                                                nowout = []
                                                assobj = [False for c in reflist]
                                                for pp in range(len(objlist)):
                                                    # convert coordinates from obj to ref
                                                    miecord = CartesianCoordAmpRotoTranslation(objlist[pp],(res[0][0],res[0][1]),res[0][2],res[0][3])
                                                    #print "miecord",miecord
                                                    mindist = 1e12
                                                    objmem = 0
                                                    for qq in range(len(reflist)):
                                                        #print pp, qq, assobj[qq]
                                                        #print "obj",objlist[pp],"ref",reflist[qq]
                                                        # if point is not assigned and close enough to a frame point go on
                                                        if not assobj[qq]:
                                                            mdist = math.sqrt((miecord[0]-reflist[qq][0])**2 + (miecord[1]-reflist[qq][1])**2)
                                                            if mdist < mindist:
                                                                mindist = mdist
                                                                objmem = qq
                                                    if mindist <= minassdist:
                                                        assobj[objmem] = True
                                                        nowinp.append(pp)
                                                        nowout.append(objmem)
                                                    #print "dist",mindist,minassdist,objmem
                                                # now we check how many associated points we have
                                                #print len(nowout),reqmtch
                                                if len(nowout) >= reqmtch:
                                                    return nowinp,nowout,res
                                            else:
                                                pass
#                                                print "Not rotable."
    return None
    

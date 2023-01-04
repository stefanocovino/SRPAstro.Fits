""" Utility functions and classes for SRP

Context : SRP
Module  : SRPGW
Version : 2.0.1
Author  : Stefano Covino
Date    : 07/03/2020
E-mail  : stefano.covino@inaf.it
URL:    : http://www.brera.inaf.it/utenti/covino

Usage   :

Remarks :

History : (23/08/2017) First version.
        : (22/05/2019) Better computation of the background and its error
        : (03/07/2020) Minor corrections after library update.
"""


#import SRPGW as GW
import numpy
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.table import hstack
from SRPFITS.Photometry.ApErr import ApErr



#def ApyPhot(x,y,data,rds=(5,10,15),backgr=True,gain=1.,ron=0.):
#    aperture = CircularAperture((x,y), r=rds[0])
#    if backgr:
#        baperture = CircularAnnulus((x,y), r_in=rds[1], r_out=rds[2])
#        baperture_masks = baperture.to_mask(method='center')
#        baperture_data = baperture_masks[0].multiply(data)
#        mask = baperture_masks[0].data
#        baperture_data_1d = baperture_data[mask > 0]
#    else:
#        bg = 0.0
#        ebg = 0.0
#        bpix = 0
#    #
#    flux_table = aperture_photometry(data, aperture, method='subpixel')
#    npix = aperture.area()
#    #
#    if backgr:
#         _, bg,ebg = sigma_clipped_stats(baperture_data_1d)
#        bpix = baperture.area()
#    bkg_sum = bg * npix
#    flux = flux_table['aperture_sum'] - bkg_sum
#    #
#    fluxerr = ApErr(flux_table['aperture_sum']/npix,npix,bpix,bg,ebg**2,k=np.pi/2))
#    #
#    return flux,fluxerr
    


def ApyPhot(x,y,data,rds=(5,10,15),backgr=True,gain=1.,ron=0.):
    if type(x) == numpy.float64 or type(y) == numpy.float64:
        x = [x,]
        y = [y,]
    #
    pos = [(i,l) for i,l in zip(x,y)]
    apertures = CircularAperture(pos, r=rds[0])
    if backgr:
        bapertures = CircularAnnulus(pos, r_in=rds[1], r_out=rds[2])
        bapertures_masks = bapertures.to_mask(method='center')
        bpix = bapertures.area
    else:
        bg = 0.0
        ebg = 0.0
        bpix = 0
    #
    if backgr:
        bg = []
        ebg = []
        for mask in bapertures_masks:
            baperture_data = mask.multiply(data)
            baperture_data_1d = baperture_data[mask.data > 0]
            _, bgmed,stdev = sigma_clipped_stats(baperture_data_1d)
            bg.append(bgmed)
            ebg.append(stdev)
    bkg = numpy.array(bg)
    ebkg = numpy.array(ebg)
    #
    phot = aperture_photometry(data, apertures, method='subpixel')
    npix = apertures.area
    phot['annulus_median'] = bkg
    phot['aper_bkg'] = bkg * npix
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    #
    flux = phot['aper_sum_bkgsub']
    if backgr:
        kp = numpy.pi/2
    else:
        kp = 0.
    fluxerr = numpy.sqrt(ApErr(phot['aperture_sum']/npix,npix,bpix,bkg,ebkg**2,k=kp))
    #
    return flux,fluxerr

""" Init file for Fits

Context : SRP
Module  : SRPFits
Version : 1.2.0
Author  : Stefano Covino
Date    : 29/03/2022
E-mail  : stefano.covino@inaf.it
URL     : http://www.me.oa-brera.inaf.it/utenti/covino


Usage   : to be imported

Remarks :

History : (28/09/2010) First named version.
        : (01/10/2010) New function.
        : (06/01/2010) New function: AddHeaderEntry.
        : (01/04/2015) FitsTabsAppend added.
        : (27/05/2017) New functions added.
        : (14/01/2021) WCSRotationDeg and WCSPixelScale added.
"""



__all__ = ['AddHeaderComment', 'AddHeaderEntry', 'FitsConstant', 'FitsImageClass',
           'FitsTabsAppend', 'GetData', 'GetHeader', 'GetHeaderValue', 'GetSpectrum',
            'GetSpectrumPosition', 'GetWCS', 'IsFits', 'WCSPixelScale', 'WCSRotationDeg']



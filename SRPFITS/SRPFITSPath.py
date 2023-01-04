""" Utility functions and classes for SRP

Context : SRP
Module  : System.py
Version : 1.1.0
Author  : Stefano Covino
Date    : 02/06/2021
E-mail  : stefano.covino@brera.inaf.it
URL:    : http://www.merate.mi.astro.it/utenti/covino

Usage   : to be imported

Remarks :

History : (30/05/2017) First version.
        : (29/08/2017) Bug correction.
        : (02/06/2021) Update.
"""


import SRPFITS

def SRPFITSPath():
    for i in SRPFITS.__path__:
        if i.find('SRPFITS') > 0:
            return i
    return None

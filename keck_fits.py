# -*- coding: utf-8 -*-
"""
filehandling
Created on Wed Apr 09 14:23:52 2014

@author: jholt
"""
from astropy.io import fits

def get_array_and_header(file1):

        if isinstance(file1, str):
            hdulist = fits.open(file1)
            hdr = hdulist[0].header
            data = hdulist[0].data
            # data, hdr = pf.getdata(file1, 0, header=True, ignore_missing_end=True, verify='ignore')
            filename = file1

        else:
            # allowing duck typing after that
            data = file1
            hdr = 'error reading header'
            filename = 'file.fits'

        # formatting for output files  
        filename = filename.rstrip('.fits')
        if len(filename.split('/')) > 1:
            filename = filename.split('/')[-1]

        return data, hdr, filename

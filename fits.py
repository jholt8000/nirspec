# -*- coding: utf-8 -*-
"""
filehandling
Created on Wed Apr 09 14:23:52 2014

@author: jholt
"""
import pyfits as pf

#import astropy.io.fits

class Handle_fits():
    
    @staticmethod
    def get_array_and_header(file1):
        if isinstance(file1,str):
            data, hdr = pf.getdata(file1, 0, header=True, ignore_missing_end=True, verify='ignore')
            filename = file1
        else: 
            # allowing duck typing after that
            data=file1
            hdr = 'error reading header'
            filename = 'sci_reduction'
            
        # formatting for output files  
        filename = filename.rstrip('.fits')
        if len(filename.split('/')) > 1: 
            filename = filename.split('/')[-1]  
            
            
        return data, hdr, filename
    

 
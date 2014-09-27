#!/home/koaadmin/Ureka/variants/common/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 16:08:28 2013

@author: jholt
"""

import numpy as np 

import reduce_order
reload(reduce_order)

import nirspec_wavelength_utils
reload(nirspec_wavelength_utils)
import array_manipulate
reload(array_manipulate)
import fits
reload(fits)

import astro_math
reload(astro_math)
import nirspec_util
reload(nirspec_util)

import twod_lambda_fit
reload(twod_lambda_fit)
'''
reduce NIRSPEC data
'''     

class Main():
    ''' 
    This is the main reduction script, it is mostly a wrapper around the reduce_order 
    program 
    
    Parameters:
    --------------
    sci: str
        FITS file name of science
    flat: str
        FITS file name of flat that matches the science
    dark: str
        (optional) FITS files name of dark that matches the science
    do_extract: Boolean
        if do_extract=True, science is extracted with sky subtraction
    do_darksub: Boolean
        subtract the supplied dark, if do_darksub=True, a dark must be supplied
    show_plot: Boolean
        display the plots with pylab, usually used for troubleshooting
    cosmic_clean: Boolean
        clean the cosmic rays -- this is recommended for order location accuracy
    max_iter: int
        input variable to cosmic ray 
    sig_clip: float
        input variable to cosmic ray 
    sig_frac: float
        input variable to cosmic ray 
    obj_lim: float
        input variable to cosmic ray 
    ext_height: int
        height in pixels of continuum to be extracted
    sky_distance: int
        length in pixels away from center of continuum to extract sky
    write_fits: Boolean
        write out fits output of each individual rectified order
    sky_height: int
        height in pixels of sky to be extracted
    rawpath: str
        location of input files
    outpath: str
        desired location of output files
    write_plots: Boolean
        write PNG snapshots of output data
    sky_sigma: float
        sky line threshold for signal/noise
        
    Returns:
    -------------
    FITS files for each cut-out and rectified order
    PNG files for each plot
    
    Examples:
        Instansiate reduction object
        > reduction_object=nirspecOO.Main('NS.20120106.16518.fits','NS.20120106.05719.fits',do_extract=True,cosmic_clean=True,write_fits=False,write_plots=False,show_plot=True)
        Reduce array         
        > reduction_object.reduce_nirspec()
        
    '''
    
    def __init__(self, sci, flat, dark = '', do_extract = False, do_darksub = True, 
                 show_plot = True, cosmic_clean = False, 
                 max_iter = 3, sig_clip = 5.0, sig_frac = 0.3, obj_lim = 5.0, 
                 ext_height = 3, sky_distance = 5, write_fits = False, 
                 sky_height = 5, rawpath = '', outpath = '', write_plots = False,  
                 sky_sigma = 2.25):
             
            ''' Initialize reduction'''
             
            self.sci, self.sciheader, self.sciname = fits.Handle_fits.get_array_and_header(sci)

            if flat:             
                self.flat, self.flatheader, self.flatname = fits.Handle_fits.get_array_and_header(flat) 
            
            else:
                print 'I need a flat to reduce science frame'
                return 
                 
            self.dark = dark
            self.do_extract = do_extract
            self.do_darksub = do_darksub
            self.show_plot = show_plot
            self.cosmic_clean = cosmic_clean
            self.ext_height = ext_height
            self.sky_distance = sky_distance
            self.write_fits = write_fits
            self.sky_height = sky_height
            self.rawpath = rawpath
            self.outpath = outpath
            self.write_plots = write_plots
            self.sky_sigma = sky_sigma
            self.max_iter = max_iter
            self.sig_clip = sig_clip
            self.sig_frac = sig_frac
            self.obj_lim = obj_lim 

            #self.reduce_nirspec()
    
    def reduce_nirspec(self):
             
        ''' takes input raw science and raw flat and produces a rectified, 
        extracted, wavelength calibrated spectrum for quick-look purposes'''
     
        # instantiate nirspec specific bookkeeping and header utils
        self.nh = nirspec_util.Nirspec_header(self.sciheader)       
        
        nb = nirspec_util.Nirspec_bookkeeping(self.sciname,
                                              outpath = self.outpath) 
                                                        
        self.logger = nb.setup_logger()
        
        self.logger.info(str('Starting Data Reduction for science file '+self.sciname)) 

        # Get precanned starting self.order_nums & self.order_num sizes for each filter
        self.data_dict = self.nh.get_data_dict()
        
        if self.data_dict['startord'] == 9999: 
            self.logger.error('Filter name '+str(self.header['filname'])+' is not ready')
            return 
            
        # instanciate sciobj and flatobj objects
        sciobj = array_manipulate.SciArray(self.sci)
        flatobj = array_manipulate.FlatArray(self.flat)
        
        # this should be in a method
        if self.dark:
            sciobj.data = sciobj.data - self.dark
            flatobj.data = flatobj.data - self.dark
            
        #inherit self (current reduction) attributes to sciobj and flatobj
        sciobj.__dict__.update(self.__dict__)
        flatobj.__dict__.update(self.__dict__)
             
        if self.cosmic_clean: 
            self.logger.info(str('cosmic ray cleaning using LA Cosmic'))

            # cosmic clean sciobj.data and flatobj.data #
            sciobj.cosmic()
            flatobj.cosmic()
                    
        if self.sciheader['dispers'] == 'low': 
            self.logger.info('low dispersion data')
            self.low_disp = True
            
        else: self.low_disp = False
   
        # this starting self.order_num number is used in getpos function to find possible 
        # starting point for flatobj.centroid calculation    
        self.order_num = self.data_dict['startord']     
            
        # shift and subtract flat from itself 
        # to get top and bottom starting locations
        #sets flatobj.tops, flatobj.bots 
        flatobj.make_tops_bots()
         
        sciobj.do_extract_orig = sciobj.do_extract
 
        #initialize variables and lists
        lhs_top=0
        all_bad_sol=[]
        all_order_objects=[]
        orig_pix_x=[]
        order_number=[]
        matched_sky_line=[]
        
        # go through each of the orders and reduce, if edges found
        # reduce => find edges, trace edges, cut out order, rectify order
        #           flatfield correct, extract 1d 
        while self.order_num > 0 and lhs_top < 1034:
            all_nirspec_data = reduce_order.reduce_order(self, sciobj, flatobj)
            if len(all_nirspec_data[0]) > 0:
                order, sciorder, lineobj, flatobj, traceobj, found_wavelength = all_nirspec_data[0]
                all_order_objects.append(all_nirspec_data[0])

                if found_wavelength:       
                    # store wavlength solution
                    orig_pix_x.append((np.array(lineobj.matchesidx)))
                    order_number.append((np.zeros(len(lineobj.matchesidx))+order))
                    matched_sky_line.append((np.array(lineobj.matchesohx)))
                    
                else:
                    all_bad_sol.append((lineobj.matchesdx, lineobj.matchesohx, lineobj.matchesidx, order))
                    
            if all_nirspec_data[1] > 1034:
                self.order_num=self.order_num-100 #break out early, don't go trhough all orders to order=1
            else:
                self.order_num=self.order_num-1
        if len(orig_pix_x) > 0:    
            orig_pix_x=np.hstack(orig_pix_x)
            order_number=1./np.hstack(order_number)
            matched_sky_line=np.hstack(matched_sky_line)

        p1 = twod_lambda_fit.twodfit(np.array(orig_pix_x), 
                                     np.array(order_number), 
                                     np.array(matched_sky_line), logger=self.logger,
                                     lower_len_points=10., sigma_max=0.5)
        
        
        for (order_num, sciorder, lineobj, flatobj, traceobj, found_wavelength) in all_order_objects:  
            ### all this belongs elsewhere ###
            newdx = np.arange(1024)
            newy=1./order_num
            newoh=np.ravel(p1[0] + p1[1]*newdx + p1[2]*newdx**2 + p1[3]*newy + p1[4]*newdx*newy + p1[5]*(newdx**2)*newy )
           
            sciorder.dx_2dfit = astro_math.conv_Ang_to_mu(newoh) 
            sciorder.dx = astro_math.conv_Ang_to_mu(sciorder.dx)
            lineobj.matchesohx = astro_math.conv_Ang_to_mu(lineobj.matchesohx)
            lineobj.bigohx = astro_math.conv_Ang_to_mu(lineobj.bigohx)

            #### fit a line to final x-axis
            
            nb.make_nirspec_final_FITS_and_plots(self, order_num, sciorder, lineobj, traceobj, flatobj, sciobj)
                              
        
        self.logger.info('Finished quicklook reduction for science file '+self.sciname)
        
        nb.close_logger()

        
    
    
        
         
         
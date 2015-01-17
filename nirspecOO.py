#!/home/koaadmin/Ureka/variants/common/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 16:08:28 2013

@author: jholt
"""

import numpy as np

import reduce_order
reload(reduce_order)

import wavelength_utils
reload(wavelength_utils)
import array_manipulate
reload(array_manipulate)
import keck_fits
reload(keck_fits)

import astro_math
reload(astro_math)
import nirspec_util
reload(nirspec_util)

import twod_lambda_fit
reload(twod_lambda_fit)
"""
reduce NIRSPEC data
"""


class Main():
    """ 
    This is the main reduction script, it is mostly a wrapper around the reduce_order 
    program 
    
    Parameters:
    --------------
    :param sci_name: str
        FITS file name of science
    :param flat_name: str
        FITS file name of flat that matches the science
    :param dark_name: str
        (optional) FITS files name of dark that matches the science
    :param do_extract: Boolean
        if do_extract=True, science is extracted with sky subtraction
    :param do_darksub: Boolean
        subtract the supplied dark, if do_darksub=True, a dark must be supplied
    :param show_plot: Boolean
        display the plots with pylab, usually used for troubleshooting
    :param cosmic_clean: Boolean
        clean the cosmic rays -- this is recommended for order location accuracy
    :param max_iter: int
        input variable to cosmic ray 
    :param sig_clip: float
        input variable to cosmic ray 
    :param sig_frac: float
        input variable to cosmic ray 
    :param obj_lim: float
        input variable to cosmic ray 
    :param ext_height: int
        height in pixels of continuum to be extracted
    :param sky_distance: int
        length in pixels away from center of continuum to extract sky
    :param sky_height: int
        height in pixels of sky to be extracted
    :param rawpath: str
        location of input files
    :param outpath: str
        desired location of output files
    :param write_fits: Boolean
        write out fits output of each individual rectified order
    :param write_plots: Boolean
        write PNG snapshots of output data
    :param sky_sigma: float
        sky line threshold for signal/noise
        
    Returns:
    -------------
    FITS files for each cut-out and rectified order - if write_fits=True
    PNG files for each plot - if write_plots=True
    
    Examples:
        Instantiate reduction object
        > reduction_1=nirspecOO.Main('NS.20120106.16518.fits','NS.20120106.05719.fits',do_extract=True,cosmic_clean=True,write_fits=False,write_plots=False,show_plot=True)
        Reduce array         
        > reduction_1.reduce_nirspec()
        
    """

    def __init__(self, sci_name, flat_name='', flat_array=[], dark_name='', dark_array=[], do_extract=True,\
                 do_darksub=True, show_plot=True, cosmic_clean=True, max_iter=3, sig_clip=5.0, sig_frac=0.3,\
                 obj_lim=5.0, ext_height=3, sky_distance=5, sky_height=5, rawpath='', outpath='', write_plots=False, \
                 write_fits=False, sky_sigma=2.25, traceWidth=10, backgroundWidth=30, traceMean=True, traceLast=False,\
                 traceDelta=1.9, verbose=True):

        """ Initialize reduction"""

        self.flat_array = flat_array
        self.dark_array = dark_array
        self.dark_name = dark_name
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
        self.backgroundWidth = backgroundWidth
        self.traceMean = traceMean
        self.traceLast = traceLast
        self.traceDelta = traceDelta
        self.traceWidth = traceWidth

        self.sci, self.sciheader, self.sciname = keck_fits.get_array_and_header(sci_name)

        if flat_name:
            self.flat, self.flatheader, self.flatname = keck_fits.get_array_and_header(flat_name)

        else:
            if flat_array.any():
                print 'using flat array '
                self.flat = flat_array
                self.flatname = 'unknown_flat_name'
            else:
                print 'I need a flat to reduce science frame'
                return

        # instantiate nirspec specific bookkeeping and header utils
        self.nh = nirspec_util.NirspecHeader(self.sciheader)
        self.data_dict = self.nh.get_data_dict()

        # this starting self.order_num number is used in getpos function to find possible
        # starting point for flatObj.centroid calculation
        self.order_num = self.data_dict['startord']

        self.nb = nirspec_util.NirspecBookkeeping(self.sciname,
                                                  outpath=self.outpath)
        self.logger = self.nb.setup_logger(self.sciname, self.outpath, verbose)

        if self.sciheader['dispers'] == 'low':
            self.logger.info('low dispersion data')
            self.low_disp = True

        else:
            self.low_disp = False

    def reduce_nirspec(self):

        """ takes input raw science and raw flat and produces a rectified, 
        extracted, wavelength calibrated spectrum for quick-look purposes

        """

        self.logger.info(str('Starting Data Reduction for science file ' + self.sciname))

        # Get pre-canned starting self.order_nums & self.order_num sizes for each filter

        if self.data_dict['startord'] == 9999:
            self.logger.error('Filter name ' + str(self.header['filname']) + ' is not ready')
            return

        # instantiate sciObj and flatObj objects
        sciObj = array_manipulate.SciArray(self.sci)
        flatObj = array_manipulate.FlatArray(self.flat)

        # this should be in a method
        if self.dark_name:
            sciObj.data -= self.dark_name
            flatObj.data -= self.dark_name

        if self.cosmic_clean:
            self.logger.info(str('cosmic ray cleaning using LA Cosmic'))

            # cosmic clean sciObj.data and flatObj.data #
            sciObj.cosmic(sig_clip=self.sig_clip, sig_frac=self.sig_frac,
                                 obj_lim=self.obj_lim)
            flatObj.cosmic(sig_clip=self.sig_clip, sig_frac=self.sig_frac,
                                 obj_lim=self.obj_lim)

        import pylab as pl
        pl.figure(1)
        pl.imshow(flatObj.data)
        pl.figure(2)
        pl.imshow(sciObj.data)
        pl.figure(3)
        pl.imshow(flatObj.tops)
        pl.figure(4)
        pl.imshow(flatObj.bots)

        # initialize variables and lists
        lhs_top = 0
        all_order_objects = []
        orig_pix_x = []
        order_number_array = []
        matched_sky_line = []

        # go through each of the orders and reduce, if edges found
        # reduce => find edges, trace edges, cut out order, rectify order
        # flatfield correct, extract 1d 
        while self.order_num > 0 and lhs_top < 1034:

            # ## This is the main reduction of each order ###
            reduced_order_object = reduce_order.Reduce_order(self.order_num, logger = self.logger,
                                                             ext_height=self.ext_height, sky_distance=self.sky_distance,
                                                             sky_height=self.sky_height, sky_sigma=self.sky_sigma, hdr_obj=self.nh,
                                                             sciobj=sciObj, flatobj=flatObj, sci_data=[], flat_data=[],
                                                             order_threshold=self.data_dict['order_threshold'],
                                                             do_extract=self.do_extract, padding=self.data_dict['padding'],
                                                             traceWidth=self.traceWidth, backgroundWidth=self.backgroundWidth,
                                                             traceMean=self.traceMean, traceLast=False, traceDelta=self.traceDelta)
            reduced_order_object.reduce_order()

            lhs_top = reduced_order_object.lhs_top

            if reduced_order_object.lineobj is not None:
                # ## add the reduced_order_object with all arrays and info about reduction to master tuple
                all_order_objects.append(reduced_order_object)
                if reduced_order_object.found_wavelength:
                    # store wavelength solution
                    orig_pix_x.append((np.array(reduced_order_object.lineobj.matchesidx)))
                    matched_sky_line.append((np.array(reduced_order_object.lineobj.matchesohx)))
                    # ## create an array with the same length as the number of sky lines identified and matched
                    order_number_array.append(np.zeros(len(reduced_order_object.lineobj.matchesdx)) + self.order_num)

            if lhs_top > 1034:
                # don't go through all orders to order_num=1 once we are off detector
                self.order_num -= 100
            else:
                self.order_num -= 1

        if len(orig_pix_x) > 0:
            # reduce dimensions of each array of matched sky lines from each order into
            # single 1d array, "flatten" array

            orig_pix_x_stack = np.hstack(orig_pix_x)
            # the 2d solution fits the inverse order number better
            order_number_array_stack = 1. / np.hstack(order_number_array)
            matched_sky_line_stack = np.hstack(matched_sky_line)

            p1, newoh, dataZZ = twod_lambda_fit.twodfit(np.array(orig_pix_x_stack),
                                         np.array(order_number_array_stack),
                                         np.array(matched_sky_line_stack), logger=self.logger,
                                         lower_len_points=10., sigma_max=0.3)
        else:
            p1 = []

        # Go back through each order and apply the 2d wavelength fit found
        for order_object in all_order_objects:

            newoh = twod_lambda_fit.applySolution(order_object, p1)

            # ## make plots and output FITS files
            self.nb.make_nirspec_final_FITS_and_plots(self, order_object, order_object.sciorder, order_object.lineobj,
                                                     order_object.flatobj, order_object.sciobj,
                                                      newoh)

        self.logger.info('Finished quicklook reduction for science file ' + self.sciname)

        self.nb.close_logger()
        #return np.array(orig_pix_x_stack), np.array(order_number_array_stack),  np.array(matched_sky_line_stack), self.logger

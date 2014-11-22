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
    :param write_fits: Boolean
        write out fits output of each individual rectified order
    :param sky_height: int
        height in pixels of sky to be extracted
    :param rawpath: str
        location of input files
    :param outpath: str
        desired location of output files
    :param write_plots: Boolean
        write PNG snapshots of output data
    :param sky_sigma: float
        sky line threshold for signal/noise
        
    Returns:
    -------------
    FITS files for each cut-out and rectified order
    PNG files for each plot
    
    Examples:
        Instantiate reduction object
        > reduction_1=nirspecOO.Main('NS.20120106.16518.fits','NS.20120106.05719.fits',do_extract=True,cosmic_clean=True,write_fits=False,write_plots=False,show_plot=True)
        Reduce array         
        > reduction_1.reduce_nirspec()
        
    """

    def __init__(self, sci_name, flat_name, dark_name='', do_extract=True, do_darksub=True,
                 show_plot=True, cosmic_clean=True,
                 max_iter=3, sig_clip=5.0, sig_frac=0.3, obj_lim=5.0,
                 ext_height=3, sky_distance=5, write_fits=False,
                 sky_height=5, rawpath='', outpath='', write_plots=False,
                 sky_sigma=2.25):

        """ Initialize reduction"""

        self.sci, self.sciheader, self.sciname = fits.Handle_fits.get_array_and_header(sci_name)

        if flat_name:
            self.flat, self.flatheader, self.flatname = fits.Handle_fits.get_array_and_header(flat_name)

        else:
            print 'I need a flat to reduce science frame'
            return

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

        # self.reduce_nirspec()

    def reduce_nirspec(self):

        """ takes input raw science and raw flat and produces a rectified, 
        extracted, wavelength calibrated spectrum for quick-look purposes
        :type self: object
        """

        # instantiate nirspec specific bookkeeping and header utils
        self.nh = nirspec_util.NirspecHeader(self.sciheader)

        nb = nirspec_util.NirspecBookkeeping(self.sciname,
                                             outpath=self.outpath)

        self.logger = nb.setup_logger()

        self.logger.info(str('Starting Data Reduction for science file ' + self.sciname))

        # Get pre-canned starting self.order_nums & self.order_num sizes for each filter
        self.data_dict = self.nh.get_data_dict()

        if self.data_dict['startord'] == 9999:
            self.logger.error('Filter name ' + str(self.header['filname']) + ' is not ready')
            return

            # instanciate sciobj and flatobj objects
        sciobj = array_manipulate.SciArray(self.sci)
        flatobj = array_manipulate.FlatArray(self.flat)

        # this should be in a method
        if self.dark_name:
            sciobj.data -= self.dark_name
            flatobj.data -= self.dark_name

        # inherit self (current reduction) attributes to sciobj and flatobj
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

        else:
            self.low_disp = False

        # this starting self.order_num number is used in getpos function to find possible 
        # starting point for flatobj.centroid calculation    
        self.order_num = self.data_dict['startord']

        # shift and subtract flat from itself 
        # to get top and bottom starting locations
        # sets flatobj.tops, flatobj.bots 
        flatobj.make_tops_bots()

        sciobj.do_extract_orig = sciobj.do_extract

        # initialize variables and lists
        lhs_top = 0
        all_bad_sol = []
        all_order_objects = []
        all_order_matches = []
        orig_pix_x = []
        order_number_array = []
        matched_sky_line = []

        # go through each of the orders and reduce, if edges found
        # reduce => find edges, trace edges, cut out order, rectify order
        # flatfield correct, extract 1d 
        while self.order_num > 0 and lhs_top < 1034:

            reduced_order_object = reduce_order.Reduce_order(self, sciobj, flatobj)
            reduced_order_object.reduce_order()

            lhs_top = reduced_order_object.traceobj.lhs_top

            if reduced_order_object.lineobj is not None and self.do_extract:

                all_order_objects.append(reduced_order_object)

                order_number_array.append(np.zeros(len(reduced_order_object.lineobj.matchesdx)) +self.order_num)

                if reduced_order_object.found_wavelength:
                    # store wavelength solution
                    orig_pix_x.append((np.array(reduced_order_object.lineobj.matchesidx)))
                    all_order_matches.append(
                        (np.zeros(len(reduced_order_object.lineobj.matchesidx)) + self.order_num))
                    matched_sky_line.append((np.array(reduced_order_object.lineobj.matchesohx)))

                else:
                    all_bad_sol.append((reduced_order_object.lineobj.matchesdx,
                                        reduced_order_object.lineobj.matchesohx,
                                        reduced_order_object.lineobj.matchesidx, self.order_num))

            if lhs_top > 1034:
                # break out early, don't go through all orders to order_num=1
                self.order_num -= 100
            else:
                self.order_num -= 1

                # need to add a sanity checker to first fit some function to each order's OH and IDs
                # opx, ona, msl = sanity_check(orig_pix_x, order_number_array, matched_sky_line)

        if len(orig_pix_x) > 0:
            # reduce dimensions of each array into single 1d array, "flatten" array
            orig_pix_x_stack = np.hstack(orig_pix_x)
            #try:
            order_number_array_stack = 1. / np.hstack(order_number_array)
            print 'onas=',order_number_array_stack
            matched_sky_line_stack = np.hstack(matched_sky_line)
            print 'msls=',matched_sky_line_stack
            print "orig_pix_x_stack=",orig_pix_x_stack
            p1 = twod_lambda_fit.twodfit(np.array(orig_pix_x_stack),
                                         np.array(order_number_array_stack),
                                         np.array(matched_sky_line_stack), logger=self.logger,
                                         lower_len_points=10., sigma_max=0.5)


            print '2d lambda not working'

        for order_object in all_order_objects:

            # ## all this belongs elsewhere ###
            newdx = np.arange(1024)
            newy = 1. / order_object.sciorder.order_num
            newoh = np.ravel(
                p1[0] + p1[1] * newdx + p1[2] * newdx ** 2 + p1[3] * newy + p1[4] * newdx * newy + p1[5] * (
                    newdx ** 2) * newy)

            #setattr(reduced_order_object.sciorder, "dx_2dfit", astro_math.conv_ang_to_mu(newoh))
            #reduced_order_object.sciorder.dx = astro_math.conv_ang_to_mu(reduced_order_object.sciorder.dx)
            #reduced_order_object.lineobj.matchesohx = astro_math.conv_ang_to_mu(reduced_order_object.lineobj.matchesohx)
            #reduced_order_object.lineobj.bigohx = astro_math.conv_ang_to_mu(reduced_order_object.lineobj.bigohx)

                # ### fit a line to final x-axis
            nb.make_nirspec_final_FITS_and_plots(self, order_object.sciorder, order_object.lineobj, order_object.traceobj, order_object.flatobj, order_object.sciobj, newoh)

        self.logger.info('Finished quicklook reduction for science file ' + self.sciname)

        nb.close_logger()

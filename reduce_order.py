# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:41:32 2014

@author: jholt
"""
import sys
import nirspec_util

reload(nirspec_util)
import astro_math

reload(astro_math)
import array_manipulate

reload(array_manipulate)
import wavelength_utils

reload(wavelength_utils)

import trace_order
import logging
reload(trace_order)
from fudge_constants import NirspecFudgeConstants

try:
    import numpy as np
except:
    print "ERROR: you do not have numpy installed!"
    print "       exiting...."
    sys.exit()

# This should be more stand-alone than it is, the input args should be spelled
# out better

class Reduce_order(object):

    #reduce each order found

    def __init__(self, order_num, logger = '', ext_height=3, sky_distance=5, sky_height=5, sky_sigma=2.25, hdr_obj='',\
                 header = '', order_threshold=1000,\
                 sciobj='', flatobj='', sci_data=[], flat_data=[], do_extract=True, padding=10,\
                 traceWidth=10, backgroundWidth=30, traceMean=True,\
                 traceLast=False, traceDelta=1.9):

        """
        :param order_num: number
            order number, used in finding starting locations for top and bottom of order
        :param logger: instance
            logger object created with utils script, if left blank will be created as reduction_order_<order_num>.log
        :param ext_height: number
            how many pixels above center of continuum to extract
        :param sky_distance: number
            how many pixels away from the center of the continuum to use as the sky lines
        :param sky_height: number
            range of pixels to use for finding the sky lines
        :param sky_sigma: number
            sigma above noise to use for sky line detection
        :param sciobj: instance
            science object using array_manipulate.SciArray, containing necessary reduction methods.
            if left empty, the sci_data array will be used
        :param flatobj: instance
            flat object using array_manipulate.FlatArray, containing necessary reduction methods
            if left empty, the flat_data array will be used
        :param sci_data: numpy array
            numpy array of science data to reduce
            if left empty, the sciobj parameter will be used
        :param flat_data: numpy array
            numpy array of science data used in reduction
            if left empty, the flatobj parameter will be used

        :param do_extract: Boolean
            extract?
        :raise AttributeError:
        """
        self.padding = padding
        self.order_num = order_num
        self.do_extract = do_extract
        self.ext_height = ext_height
        self.sky_distance = sky_distance
        self.sky_height = sky_height
        self.sky_sigma = sky_sigma

        self.backgroundWidth = backgroundWidth
        self.traceMean = traceMean
        self.traceLast = traceLast
        self.traceDelta = traceDelta
        self.traceWidth = traceWidth
        self.order_threshold = order_threshold

        if sciobj == '' and flatobj == '' and sci_data == []:
            print "Need sciobj and flatobj instances of array_manipulate.SciArray OR"
            print " need data numpy array (from astropy.io.fits) "
            raise AttributeError

        if header == '' and hdr_obj =='':
            print "Need header object instances of nirspec_util OR"
            print " need header (from astropy.io.fits) "
            raise AttributeError

        if not isinstance(hdr_obj, nirspec_util.NirspecHeader):
            self.hdr_obj = nirspec_util.NirspecHeader(header)
        else:
            self.hdr_obj = hdr_obj

        if not isinstance(logger, logging.Logger):
            self.logger = nirspec_util.NirspecBookkeeping.setup_logger('reduction_order_'+order_num+'.log', '.', verbose=False)
        else:
            self.logger = logger

        if isinstance(sciobj, array_manipulate.SciArray):
            self.sciobj = sciobj
        else:
            self.sciobj = array_manipulate.SciArray(sci_data)

        if isinstance(flatobj, array_manipulate.FlatArray):
            self.flatobj = flatobj
        else:
            self.flatobj = array_manipulate.FlatArray(flat_data)

        self.detector_height = self.flatobj.data.shape[0]


        #returns:
        self.found_wavelength = False
        self.lineobj = None
        self.sciorder = None

    def reduce_order(self):

        # make order-specific extraction bool
        #global sky_line_fit

        assert isinstance(self.sciobj, array_manipulate.SciArray), "need to instantiate science array"

        self.logger.info('-------Starting on Order ' + str(self.order_num) + '---------\n \n')

        # find order position on detector ####

        # Use grating eqn empirical fitting to determine the theoretical
        # starting location of the top and bottom of order as well as
        # starting wavelengths for that order.

        # location (on lhs of detector) of top of order, bottom of order, and theoretical wavelength array
        # this is a  NIRSPEC centric routine, uses header (filter, echelle and disp angles) and order_num
        self.lhs_top_theory, self.lhs_bot_theory, theory_x = self.hdr_obj.get_theory_order_pos(self.order_num)
        self.logger.info('Theory -- left bot = ' + str(int(self.lhs_bot_theory)) +
                                      ' left top = ' + str(int(self.lhs_top_theory)))


        if (self.lhs_top_theory < self.detector_height + NirspecFudgeConstants.total_chip_padding
            and self.lhs_bot_theory > -NirspecFudgeConstants.total_chip_padding):

            traceobj = trace_order.Trace(self.lhs_top_theory, self.lhs_bot_theory,
                                                     order_threshold=self.order_threshold, order_num=self.order_num,
                                                     logger=self.logger, flatobj=self.flatobj, flat_data=[],
                                                     traceWidth=self.traceWidth,
                                                     backgroundWidth=30, traceMean=self.traceMean,
                                                     traceLast=False, traceDelta=1.9)

            # instance attributes defined outside __init__, need to refactor
            self.avg_spectroid, self.top_spectroid, self.bot_spectroid, self.lhs_top, self.lhs_bot = traceobj.trace_order()

        else:
            self.logger.info('order: ' + str(self.order_num) + ' not on detector  ')
            self.lhs_top = self.lhs_top_theory
            # return to main and skip this order
            return

        # # cut out order and spatially rectify ####
        if traceobj.trace_success:

            self.logger.info('cutting out order' + str(self.order_num) + ' from science and flat')

            # ## cut out the order from the science and the flat ###
            # ## include padding on the top and bottom to ensure order is on the cutout
            # ## array and to avoid cutting into the science when order is straightened

            self.padding = astro_math.fudge_padding(self.avg_spectroid, self.padding)

            # determine highest point to decide where to cut out the array to isolate order
            highest_top = max(self.top_spectroid[0], self.top_spectroid[-10])

            # order_slice that contains only the cutout around the order
            sci_order_slice = self.sciobj.cut_out(padding=self.padding, lower=self.lhs_bot,
                                upper=highest_top)

            # order_slice that contains only the cutout around the order
            flat_order_slice = self.flatobj.cut_out(padding=self.padding, lower=self.lhs_bot,
                                 upper=highest_top)

            self.top_spectroid, self.avg_spectroid, self.bot_spectroid, order_shifted = astro_math.shift_order(self.top_spectroid,
                                                                                                self.avg_spectroid,
                                                                                                self.bot_spectroid,
                                                                                                self.padding)

            # make instances of array manip class using just the order slice
            self.sciorder = array_manipulate.SciArray(sci_order_slice)

            setattr(self.sciorder, 'order_num', self.order_num)

            flatorder = array_manipulate.FlatArray(flat_order_slice)

            self.logger.info('masking out off-order locations')

            # ## normalize flatorder data ###

            # specify the locations of actual order and off order

            on_order, off_order = flatorder.mask_order(self.top_spectroid, self.bot_spectroid)

            self.logger.info('normalizing flat order ' + str(self.order_num))

            # sets flatorder.normalized and flatorder.flat_mean

            #try:
            try:
                normalized_flat, flat_mean = flatorder.normalize(on_order, off_order, mask=True, instr="NIRSPEC")

                self.logger.info('flatfielding science')

                # sets self.sciorder.masked
                masked = self.sciorder.mask_off_order(on_order)

                # Where is this flat corrected science array used? should norm_data be data?
                #self.sciorder.flat_corrected_data = self.sciorder.data / normalized_flat
                self.sciorder.data = self.sciorder.data / normalized_flat


            except:
                self.logger.info('could not flatfield')

            # rectify the order slice, sets self.sciorder.rectified
            self.logger.info('shifting order using order edges')
            # fit a 3d poly to center-of-mass fit to smooth interpolation shift and get rid of any bumps due to uneven
            # illumination on flat order edges

            fit_spectroid, foo = astro_math.fit_poly(self.avg_spectroid[0:-100], xes=np.arange(self.sciorder.data.shape[1]), deg=2)

            rectified = self.sciorder.interp_shift(fit_spectroid, orientation='vertical', pivot='middle')
            rectified_flat = flatorder.interp_shift(fit_spectroid, orientation='vertical', pivot='middle')

            # remove the padding and start at lhs_bot to show plots in correct place
            # instance attributes defined outside __init__, need to refactor
            self.top_spectroid, self.avg_spectroid, self.bot_spectroid = astro_math.shift_order_back(self.top_spectroid, self.avg_spectroid, self.bot_spectroid, self.padding, order_shifted, self.lhs_bot)

            # overwrite data array as rectified array
            self.sciorder.data = rectified

            order_rectified = True

        else:  # could not find fit along flat order edge
            self.logger.info('WARNING: did not find a good order fit, not rectifying spatially ')
            self.logger.info('         and therefore not extracting continuum ')
            self.logger.error('Finished reducing order = ' + str(self.order_num))

            return

        # # determine sky and continuum locations  ####
        # take a the crosscut of the order and determine place off-peak for sky lines
        peak, crosscut = self.sciorder.find_peak()
        self.sciorder.__setattr__('peak', peak)
        self.sciorder.__setattr__('crosscut',crosscut)

        try:
            # i dont like this
            ext_range, sky_range_top, sky_range_bot = self.sciorder.setup_extraction_ranges(self.ext_height,
                                                                                            self.sky_distance,
                                                                                            self.sky_height,
                                                                                            peak,
                                                                                            self.logger)
            self.sciorder.__setattr__('ext_range',ext_range)
            self.sciorder.__setattr__('sky_range_top',sky_range_top)
            self.sciorder.__setattr__('sky_range_bot',sky_range_bot)

        except:
            self.logger.info('WARNING: could not find correct locations for continuum and sky')
            self.logger.info('         and therefore not extracting continuum ')
            return
            # ## Horizontal order rectification ###

        # sky_line to use in rectification using sky lines
        sky_line, sky_line_fit_success = self.sciorder.find_skyline_trace(sky_sigma=self.sky_sigma,
                                         padding = self.padding)

        try:
            sky_line_fit, foo = astro_math.fit_poly(sky_line,
                                                    xes=np.arange(self.sciorder.data.shape[0]),
                                                    deg=1)
            sky_line_fit_success = True

        except:
            sky_line_fit_success = False

        if sky_line_fit_success:  # or len(sky_line_fit[20:]) > 20:
            # # rectify spectral direction using sky line centroid fits ####
            rectified = self.sciorder.interp_shift(sky_line_fit, orientation='horizonal',
                                       pivot='peak')

            self.sciorder.data = rectified

        else:
            text = 'WARNING could not fit sky lines, not rectifying spectral dimension '
            self.do_extract = False
            self.logger.info(text)

        # Extract spectrum and subtract sky ##
        # Gain = 4 e- / ADU  and readnoise = 625 e- = 156 ADU
        if self.do_extract:

            # extract
            cont, sky, extract_status = self.sciorder.sum_extract(ext_range, sky_range_bot,
                                      sky_range_top)

            self.sciorder.__setattr__('cont',cont)

            if extract_status == 0:

                self.logger.error('could not extract order ' + str(self.order_num))
                self.logger.error('Finished reducing order = ' + str(self.order_num))
                return

            # identify sky lines with catalogue sky line locations

            # ### Wavelength skyline identification ########
            # sets self.lineobj.identify_status

            self.lineobj = wavelength_utils.LineId(theory_x, sky, False, self.logger)

            self.logger.info(
                'reduce_order: shift between sky list and real sky lines = ' + str(int(self.lineobj.lambda_shift)))

            if self.lineobj.identify_status < 1:
                self.found_wavelength = False
                self.logger.info('problem with sky line identification')
                self.logger.error(
                        'reduce_order: Could not find sky lines: only doing zeroith order wavelength calibration ')

            else:
                # find the solution between current zeroith order solution and real lambda
                p0 = np.polyfit(self.lineobj.matchesohx, self.lineobj.matchesdx, deg=1)
                disp = p0[0]
                offset = p0[1]

                self.logger.info('linear solution bw matches disp=' + str(disp) + ' offset =' + str(offset))

                # if the matches between theory and sky line list are too far off from a linear fit, dont use matches
                if NirspecFudgeConstants.disp_upper_limit > abs(disp) > NirspecFudgeConstants.disp_lower_limit:
                    self.found_wavelength = True

                else:
                    self.logger.error('bad fit: only doing zeroith order wavelength calibration ')
                    self.found_wavelength = False

                    # need to store all the pre-wavelength fixed data to make out files

        self.logger.error('Finished reducing order = ' + str(self.order_num))


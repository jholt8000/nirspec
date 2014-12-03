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
import nirspec_wavelength_utils

reload(nirspec_wavelength_utils)

import trace_order

reload(trace_order)
from fudge_constants import NirspecFudgeConstants

try:
    import numpy as np
except:
    print "ERROR: you do not have pyfits, pylab, numpy, and/or scipy installed!"
    print "       exiting...."
    sys.exit()

# This should be more stand-alone than it is, the input args should be spelled
# out better

class Reduce_order(object):
    """
    reduce each order found
    :rtype : object
    Parameters:
    --------------
    self.reduction: object
        The reduction object from the wrapper nirspecOO, has attributes: 
            self.reduction.do_extract
            self.reduction.order_num (should be order_num)
            self.reduction.logger

    self.sciobj: object
        The science array object, contains self.sciobj.data array
    self.flatobj: object
        The flat array object, contains flatoj.data
    """
    # ## attributes could be padding and dx as those change throughout ###

    def __init__(self, reduction, sciobj, flatobj):
        self.reduction = reduction
        self.sciobj = sciobj
        self.flatobj = flatobj
        self.found_wavelength = False
        self.lineobj = None
        self.traceobj = None
        self.sciorder = None

    def reduce_order(self):

        # make order-specific extraction bool
        global sky_line_fit
        do_extract = self.reduction.do_extract

        # keep this order number (self.reduction.order_num will change)
        self.reduction.order_num = self.reduction.order_num

        self.reduction.logger.info('-------Starting on Order ' + str(self.reduction.order_num) + '---------\n \n')

        # find order position on detector ####

        # Use grating eqn empirical fitting to determine the theoretical
        # starting location of the top and bottom of order as well as
        # starting wavelengths for that order.
        # lhs_bot_theory,lhs_top_theory, dx, lhs_bot, lhs_top, bot_spectroid,
        # top_spectroid, avg_spectroid, highest_top
        self.traceobj = trace_order.Trace_order_utils(self.reduction, self.flatobj, self.sciobj)

        self.traceobj.trace_order()

        # Ensure that order was traced successfully and is on detector
        if not self.traceobj.trace_success:
            self.reduction.logger.info('order ' + str(self.reduction.order_num) + ' not on array')
            # return to main self.reduction and skip this order
            return

            # # cut out order and spatially rectify ####
        if self.traceobj.trace_success:

            self.reduction.logger.info('cutting out order' + str(self.reduction.order_num) + ' from science and flat')

            # ## cut out the order from the science and the flat ###
            # ## include padding on the top and bottom to ensure order is on the cutout
            # ## array and to avoid cutting into the science when order is straightened

            self.traceobj.fudge_padding()

            # sets self.sciobj.order_slice that contains only the cutout around the order
            self.sciobj.cut_out(padding=self.traceobj.padding, lower=self.traceobj.lhs_bot,
                                upper=self.traceobj.highest_top)

            # sets self.flatobj.order_slice that contains only the cutout around the order
            self.flatobj.cut_out(padding=self.traceobj.padding, lower=self.traceobj.lhs_bot,
                                 upper=self.traceobj.highest_top)

            self.traceobj.shift_order()

            # make instances of array manip class using just the order slice
            self.sciorder = array_manipulate.SciArray(self.sciobj.order_slice)

            setattr(self.sciorder, 'order_num', self.reduction.order_num)

            flatorder = array_manipulate.FlatArray(self.flatobj.order_slice)

            # copy new order-specific dx array to be attribute of self.sciorder
            self.sciorder.dx = self.traceobj.dx

            self.reduction.logger.info('masking out off-order locations')

            # ## normalize flatorder data ###

            # specify the locations of actual order and off order
            # sets flatorder.on_order and flatorder.off_order

            flatorder.mask_order(self.traceobj.top_spectroid, self.traceobj.bot_spectroid)

            self.reduction.logger.info('normalizing flat order ' + str(self.reduction.order_num))

            # sets flatorder.normalized and flatorder.flat_mean

            try:
                flatorder.normalize(flatorder.on_order, flatorder.off_order, mask=True, instr="NIRSPEC")

                self.reduction.logger.info('flatfielding science')

                # sets self.sciorder.masked
                self.sciorder.mask_off_order(flatorder.on_order)

                # Where is this flat corrected science array used? should norm_data be data?
                self.sciorder.norm_data = self.sciorder.data / flatorder.normalized

            except:
                self.reduction.logger.info('could not flatfield')

            # rectify the order slice, sets self.sciorder.rectified
            self.reduction.logger.info('shifting order using order edges')
            # fit a 3d poly to center-of-mass fit to smooth interpolation shift and get rid of any bumps due to uneven
            # illumination on flat order edges
            fit_spectroid, foo = astro_math.fit_poly(self.traceobj.avg_spectroid,
                                                     xes=np.arange(self.sciorder.data.shape[1]),
                                                     deg=3)
            self.sciorder.interp_shift(fit_spectroid, orientation='vertical', pivot='middle')

            # remove the padding and start at lhs_bot to show plots in correct place
            self.traceobj.shift_order_back()

            # overwrite data array as rectified array
            self.sciorder.data = self.sciorder.rectified

            order_rectified = True

        else:  # could not find fit along flat order edge
            self.reduction.logger.info('WARNING: did not find a good order fit, not rectifying spatially ')
            self.reduction.logger.info('         and therefore not extracting continuum ')
            return

        # # determine sky and continuum locations  ####
        # take a the crosscut of the order and determine place off-peak for sky lines
        # sets self.sciorder.crosscut and self.sciorder.peak
        self.sciorder.find_peak(order_rectified)


        # sets self.sciorder.ext_range, self.sciorder.sky_range_bot, self.sciorder.sky_range_top
        try:
            self.sciorder.setup_extraction_ranges(self.reduction.ext_height, self.reduction.sky_distance,
                                                  self.reduction.sky_height, self.sciorder.peak,
                                                  self.reduction.order_num,
                                                  self.reduction.logger)
        except:
            self.reduction.logger.info('WARNING: could not find correct locations for continuum and sky')
            self.reduction.logger.info('         and therefore not extracting continuum ')
            return
            # ## Horizontal order rectification ###

        # sets self.sciorder.sky_line to use in rectification using sky lines
        self.sciorder.find_skyline_trace(sky_sigma=self.reduction.sky_sigma,
                                         padding=self.traceobj.padding)

        try:
            sky_line_fit, foo = astro_math.fit_poly(self.sciorder.sky_line,
                                                    xes=np.arange(self.sciorder.data.shape[0]),
                                                    deg=1)
            sky_line_fit_success = True

        except:
            sky_line_fit_success = False

        if sky_line_fit_success:  # or len(sky_line_fit[20:]) > 20:
            # # rectify spectral direction using sky line centroid fits ####
            # sets self.sciorder.rectified
            self.sciorder.interp_shift(sky_line_fit, orientation='horizonal',
                                       pivot='peak')

            # import pylab as pl
            # pl.clf()
            # pl.imshow(self.sciorder.rectified)
            self.sciorder.data = self.sciorder.rectified

        else:
            text = 'WARNING could not fit sky lines, not rectifying spectral dimension '
            self.reduction.logger.info(text)

        # Extract spectrum and subtract sky ##
        # Gain = 4 e- / ADU  and readnoise = 625 e- = 156 ADU
        if do_extract:

            # sets self.sciorder.cont, self.sciorder.skys, self.sciorder.extract_status
            self.sciorder.sum_extract(self.sciorder.ext_range, self.sciorder.sky_range_bot,
                                      self.sciorder.sky_range_top)

            # import pylab as pl
            # pl.figure(12)
            # pl.clf()
            # pl.plot(self.sciorder.cont)
            # pl.plot(self.sciorder.skys)
            # pl.show()
            if self.sciorder.extract_status == 0:
                self.reduction.logger.error('could not extract order ' + str(self.reduction.order_num))

                return


            # identify sky lines with catalogue sky line locations

            # ### Wavelength skyline identification ########
            self.lineobj = nirspec_wavelength_utils.LineId(self.reduction.low_disp,
                                                           self.sciorder.dx,
                                                           self.sciorder.skys)
            # Test commit
            # # find and apply wavelength shift ###

            # Read in the sky line list. sets self.lineobj.ohx, self.lineobj.ohy
            # skyline
            # list is determined using low_disp and fudge_constants ohlinelist
            self.lineobj.read_OH()

            # make a synthetic sky spectrum using line list information with width
            # the size of the data and with sigma = 0.2 (found empirically)
            # sets self.lineobj.fake_sky
            self.lineobj.gauss_sky(self.lineobj.ohx, self.lineobj.ohy, 0.2)

            # Cross correlate data: self.lineobj.dx and self.lineobj.fake_sky
            # with the synthetic sky, makes self.lineobj.lambda_shift
            self.lineobj.find_xcorr_shift(self.lineobj.fake_sky)

            self.reduction.logger.info(
                'shift between sky list and real sky lines = ' + str(int(self.lineobj.lambda_shift)))

            if abs(self.lineobj.lambda_shift) < NirspecFudgeConstants.max_shift_from_theory:
                self.sciorder.dx = self.sciorder.dx + self.lineobj.lambda_shift
                self.lineobj.dx = self.sciorder.dx  # i hate this, identify should have dx as a parameter -JH
                self.reduction.logger.info(
                    'applied the shift = ' + str(self.lineobj.lambda_shift) + ' to aid in sky line identification')

                # match sky lines
            # sets self.lineobj.matchesdx, self.lineobj.matchesohx, self.lineobj.matchesohy,
            # self.lineobj.bigohx, self.lineobj.bigohy, self.lineobj.identify_status
            self.lineobj.identify(self.lineobj.ohx, self.lineobj.ohy)

            if self.lineobj.identify_status < 1:
                self.reduction.logger.info('problem with sky line identification')

                if abs(self.lineobj.lambda_shift) < NirspecFudgeConstants.max_shift_from_theory:
                    self.sciorder.dx = self.sciorder.dx - self.lineobj.lambda_shift
                    self.lineobj.dx = self.sciorder.dx
                    self.reduction.logger.info(
                        'removed the shift since sky is unreliable : ' + str(-self.lineobj.lambda_shift))
                    self.reduction.logger.error(
                        'Could not find sky lines: only doing zeroith order wavelength calibration ')
                    self.found_wavelength = False

            else:
                # find the solution between current zeroith order solution and real lambda
                p0 = np.polyfit(self.lineobj.matchesohx, self.lineobj.matchesdx, deg=1)

                disp = p0[0]
                offset = p0[1]

                self.reduction.logger.info('linear solution bw matches disp=' + str(disp) + ' offset =' + str(offset))

                # if the matches between theory and sky line list are too far off from a linear fit, dont use matches
                if NirspecFudgeConstants.disp_upper_limit > abs(disp) > NirspecFudgeConstants.disp_lower_limit:
                    self.found_wavelength = True
                    # astro_math.order_wavelength_solution(self.lineobj.matchesdx, self.lineobj.matchesohx, self.sciorder.dx, order)

                else:
                    self.reduction.logger.error('bad fit: only doing zeroith order wavelength calibration ')
                    self.found_wavelength = False

                    # need to store all the pre-wavelength fixed data to make out files



# -*- coding: utf-8 -*-
"""
Created on Tue Jul 01 15:11:52 2014

@author: jholt
"""


from spectroid import spectroid
import numpy as np
import astro_math
import nirspec_util.NirspecHeader
import array_manipulate
import fudge_constants

# need to deal with how to find starting locations for top and bottom that is maybe a separate function that takes input
#  order_num, header (for filter, echelle, disp), took out and is in reduce_order


class Trace(object):
    
    def __init__(self, lhs_top_theory, lhs_bot_theory, padding=25,
                 order_threshold=1000, order_num=0,  logger='', header_util='', flatobj='', flat_data=[],
                 traceWidth=10, backgroundWidth=30, startingLocation=926, traceMean=True,
                 traceLast=False, traceDelta=1.9):

        """
        
        :param traceWidth:
            input parameter to spectroid call, traceWidth=10: is distance away, up and down, from center to search for centroid
        :param startingLocation:
            input parameter to spectroid call, starting location for centroid at x=0
        :param padding:
            how much to pad before cutting out order, this shouldn't be a param of trace_order
        :param order_threshold:
            used in zeroing in on actual order edges
        :param order_num:
            this might not need to be used here
        :param logger: instance
                    logger object created with utils script, if left blank will be created as reduction_order_<order_num>.log
        :param flatobj: instance
            flat object using array_manipulate.FlatArray, containing necessary reduction methods
            if left empty, the flat_data array will be used
        :param flat_data: numpy array
            numpy array of science data used in reduction
            if left empty, the flatobj parameter will be used

        :return:
        """
        self.backgroundWidth = backgroundWidth
        self.traceMean = traceMean
        self.traceLast = traceLast
        self.traceDelta = traceDelta
        self.lhs_top_theory = lhs_top_theory
        self.lhs_bot_theory = lhs_bot_theory
        self.order_threshold = order_threshold
        self.padding = padding
        self.startingLocation = startingLocation
        self.traceWidth = traceWidth

        if flatobj == '' and flat_data == []:
            print "Need flatobj instances of array_manipulate.SciArray OR"
            print " need data numpy array (from astropy.io.fits) "
            raise AttributeError

        if not isinstance(logger, nirspec_util.NirspecBookkeeping.setup_logger):
            self.logger = nirspec_util.NirspecBookkeeping.setup_logger('reduction_order_'+order_num+'.log', '.', verbose=False)

        if isinstance(flatobj, array_manipulate.FlatArray):
            self.flatobj = flatobj
        else:
            self.flatobj = array_manipulate.FlatArray(flat_data)

        self.trace_success = False
        self.lhs_bot_theory = lhs_bot_theory
        self.lhs_top_theory = lhs_top_theory
        self.theory_x = []


    @property
    def trace_order(self, padding=None):
        ''' finds actual edge positions using lhs_top_theory and lhs_bot_theory as starting points
            runs spectroid on those peaks, theory_x is theoretical wavelength array '''
        # this should be returning parameters

        lhs_top = self.determine_lhs_edge_pos(self.flatobj.tops, self.lhs_top_theory,
                                                   self.order_threshold)

        lhs_bot = self.determine_lhs_edge_pos(self.flatobj.bots, self.lhs_bot_theory,
                                                   self.order_threshold)

        self.logger.info('Measured -- left bot = ' + str(int(lhs_bot)))
        self.logger.info('            left top = ' + str(int(lhs_top)))

        if lhs_top > lhs_bot:

            # # spectroid on top and bottom of order location #####
            # call centroiding task, using the location of the peak and zero
            # background subtraction
            # bft = bad fit top = number of times spectroid had to self-correct

            top_spectroid, bft = spectroid(self.flatobj.tops, traceWidth=self.traceWidth,
                                                backgroundWidth=self.backgroundWidth,
                                                startingLocation=self.startingLocation, traceMean=self.traceMean,
                                                traceLast=self.traceLast, traceDelta=self.traceDelta)


            self.logger.info('had to self correct on top = ' + str(bft) + ' times ')

            if bft < fudge_constants.NirspecFudgeConstants.badfit_limit:
                traced_top = True
            try:
                # determine highest point to decide where to cut out the array to isolate order
                highest_top = max(top_spectroid[0], top_spectroid[1010])
            except:
                traced_top = False

            bot_spectroid, bfb = spectroid(self.flatobj.bots, traceWidth=self.traceWidth,
                                                backgroundWidth=self.backgroundWidth,
                                                startingLocation=self.startingLocation, traceMean=self.traceMean,
                                                traceLast=self.traceLast, traceDelta=self.traceDelta)

            self.logger.info('had to self correct on bottom = ' + str(bfb) + ' times ')

            if bfb < fudge_constants.NirspecFudgeConstants.badfit_limit:
                traced_bot = True

            try:
                lhs_bot = bot_spectroid[1]
            except:
                traced_bot = False

            if traced_top and traced_bot:
                #find the mean of the centroid along the top of the order and
                # along the bottom of the order
               avg_spectroid = (top_spectroid + bot_spectroid) / 2.
            elif traced_top:
                self.logger.info('using only top curve ')
                avg_spectroid = top_spectroid - ((lhs_top - lhs_bot) / 2) + 1.
            elif traced_bot:
                self.logger.info('using only bottom curve ')
                avg_spectroid = self.bot_spectroid + ((lhs_top - lhs_bot) / 2) + 1.

            else:
                self.logger.info('could not trace order ' + str(self.order_num))
                self.trace_success = False
                lhs_top = self.lhs_top_theory

        else:
            self.logger.info('top of order cannot be below bottom')
            self.logger.info('skipping order: ' + str(self.order_num))
            self.trace_success = False
            lhs_top = self.lhs_top_theory

        if self.trace_success and self.avg_spectroid.any():
            # update self.padding
            padding = astro_math.fudge_padding()
            # make self.avg_spectroid
            avg_spectroid = astro_math.smooth_spectroid(avg_spectroid.any())

            self.trace_success = True
            return avg_spectroid, top_spectroid, bot_spectroid, lhs_top, lhs_bot

        else:
            return [], [], [], lhs_top, lhs_bot



    def determine_lhs_edge_pos(self, edges, theory_lhs):
        """  find location of either top or bottom of the order using theoretical position
        as a starting point """
        # Find the peaks in the shifted/subracted flat file near the theoretical peaks
        lhs_edge_pos = self.nh.get_actual_order_pos(edges, theory_lhs, self.order_threshold)

        self.logger.info('found a peak at ' + str(lhs_edge_pos))

        # Ensure the theoretical location and the actual location of 
        # the top of the order are close enough
        found_edge = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs, self.order_threshold)

        if not found_edge:
            ''' lower threshold and try again '''
            self.logger.info('searching for top edge on a fainter flat')
            self.logger.info('  ---> this might affect the rectification')

            lhs_edge_pos = self.nh.get_actual_order_pos(edges, theory_lhs,
                                                                     self.order_threshold / 3)
            found_top = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs,
                                                    self.order_threshold)

            if not found_top:
                self.logger.info('Flat is too faint in this order ')
                self.logger.info('Cannot find order location')

        # return variable because it could be for top or bottom location
        return lhs_edge_pos
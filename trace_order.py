# -*- coding: utf-8 -*-
"""
Created on Tue Jul 01 15:11:52 2014

@author: jholt
"""


from spectroid import spectroid
import astro_math
import nirspec_util
import logging
import array_manipulate
import fudge_constants

# need to deal with how to find starting locations for top and bottom that is maybe a separate function that takes input
#  order_num, header (for filter, echelle, disp), took out and is in reduce_order

class Trace(object):
    
    def __init__(self, lhs_top_theory, lhs_bot_theory,
                 order_sigma=1000, order_threshold=40, order_num=0,  logger='', flatobj='', flat_data=[],
                 traceWidth=10, backgroundWidth=30, traceMean = True,
                 traceLast=False, traceDelta=1.9):

        """
        
        :param traceWidth:
            input parameter to spectroid call, traceWidth=10: is distance away, up and down, from center to search for centroid

        :param order_threshold:
            acceptable distance from theoretical order location
        :param order_sigma
            above noise for order ege
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
        self.order_sigma = order_sigma
        self.backgroundWidth = backgroundWidth
        self.traceMean = traceMean
        self.traceLast = traceLast
        self.traceDelta = traceDelta
        self.lhs_top_theory = lhs_top_theory
        self.lhs_bot_theory = lhs_bot_theory
        self.order_threshold = order_threshold
        self.traceWidth = traceWidth
        self.order_num = order_num

        if flatobj == '' and flat_data == []:
            print "Need flatobj instances of array_manipulate.SciArray OR"
            print " need data numpy array (from astropy.io.fits) "
            raise AttributeError("Need flatobj instances of array_manipulate.SciArray OR need data numpy array (from astropy.io.fits) ")

        if not isinstance(logger, logging.Logger):
            self.logger = nirspec_util.NirspecBookkeeping.setup_logger('reduction_order_'+order_num+'.log', '.', verbose=False)
        else:
            self.logger = logger

        if isinstance(flatobj, array_manipulate.FlatArray):
            self.flatobj = flatobj
        else:
            self.flatobj = array_manipulate.FlatArray(flat_data)

        self.trace_success = False
        self.lhs_bot_theory = lhs_bot_theory
        self.lhs_top_theory = lhs_top_theory

    def trace_order(self):
        ''' finds actual edge positions using lhs_top_theory and lhs_bot_theory as starting points
            runs spectroid on those peaks'''
        # this should be returning parameters
        self.logger.info("searching for lhs tops")
        lhs_top = self.determine_lhs_edge_pos(self.flatobj.tops, self.lhs_top_theory)
        self.logger.info("searching for lhs bots")
        lhs_bot = self.determine_lhs_edge_pos(self.flatobj.bots, self.lhs_bot_theory)

        self.logger.info('Measured -- left bot = ' + str(int(lhs_bot)))
        self.logger.info('            left top = ' + str(int(lhs_top)))

        if lhs_top > lhs_bot:

            # # spectroid on top and bottom of order location #####
            # call centroiding task, using the location of the peak and zero
            # background subtraction
            # bft = bad fit top = number of times spectroid had to self-correct
            try:
                top_spectroid, bft = spectroid(self.flatobj.tops, traceWidth=self.traceWidth,
                                                backgroundWidth=self.backgroundWidth,
                                                startingLocation=lhs_top, traceMean=self.traceMean,
                                                traceLast=self.traceLast, traceDelta=self.traceDelta)

            except:
                self.logger.error('Could not trace starting at '+str(lhs_top))
                top_spectroid, bft = (0, 9999)

            self.logger.info('had to self correct on top = ' + str(bft) + ' times ')

            if bft < fudge_constants.NirspecFudgeConstants.badfit_limit:
                try:
                    traced_top = True
                except:
                    traced_top = False
            else:
                traced_top = False
            try:
                bot_spectroid, bfb = spectroid(self.flatobj.bots, traceWidth=self.traceWidth,
                                                backgroundWidth=self.backgroundWidth,
                                                startingLocation=lhs_bot, traceMean=self.traceMean,
                                                traceLast=self.traceLast, traceDelta=self.traceDelta)
            except:
                self.logger.error('Could not trace starting at '+str(lhs_bot))
                bot_spectroid, bfb = (0, 9999)

            self.logger.info('had to self correct on bottom = ' + str(bfb) + ' times ')

            if bfb < fudge_constants.NirspecFudgeConstants.badfit_limit:

                try:
                    lhs_bot = bot_spectroid[1]
                    traced_bot = True
                except:
                    traced_bot = False
            else:
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
                avg_spectroid = bot_spectroid + ((lhs_top - lhs_bot) / 2) + 1.

            else:
                self.logger.info('could not trace order ' + str(self.order_num))
                self.trace_success = False
                lhs_top = self.lhs_top_theory

        else:
            self.logger.info('top of order cannot be below bottom')
            self.logger.info('skipping order: ' + str(self.order_num))
            self.trace_success = False
            lhs_top = self.lhs_top_theory
            traced_top = False
            traced_bot = False

        if traced_top or traced_bot:
            # avg_spectroid will be a numpy array
            self.trace_success = True
            return avg_spectroid, top_spectroid, bot_spectroid, lhs_top, lhs_bot

        else:
            return [], [], [], lhs_top, lhs_bot

    def determine_lhs_edge_pos(self, edges, theory_lhs):
        """  find location of either top or bottom of the order using theoretical position
        as a starting point """
        # Find the peaks in the shifted/subracted flat file near the theoretical peaks
        lhs_edge_pos = astro_math.get_actual_order_pos(edges, theory_lhs, self.order_sigma)

        self.logger.info('found a peak at ' + str(lhs_edge_pos))

        # Ensure the theoretical location and the actual location of 
        # the top of the order are close enough
        found_edge = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs, self.order_threshold)

        if not found_edge:
            ''' lower threshold and try again '''
            self.logger.info('searching for top edge on a fainter flat')
            self.logger.info('  ---> this might affect the rectification')

            lhs_edge_pos = astro_math.get_actual_order_pos(edges, theory_lhs,
                                                                     self.order_sigma / 3)
            found_top = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs,
                                                    self.order_threshold / 3)

            if not found_top:
                self.logger.info('Flat is too faint in this order ')
                self.logger.info('Cannot find order location')
                return 9999

        # return variable because it could be for top or bottom location
        return lhs_edge_pos
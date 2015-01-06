# -*- coding: utf-8 -*-
"""
array manipulation utilities

@author: jholt 20140315
"""
import numpy as np
import scipy as sp
import astro_math

reload(astro_math)
try:
    from scipy.signal import argrelextrema
except:
    print 'need to update scipy to get argrelextrema'
from spectroid import spectroid
from fudge_constants import NirspecFudgeConstants


class BaseArray(object):
    """
    Base class for all science and flat array objects

    """

    def __init__(self, data):
        self.data = data
        self.interp_shifted = False
        self.cosmic_cleaned = False

    def __multiply__(self, mult):
        """ multiply array by mult = constant or mult = array """
        return self.data * mult

    def __div__(self, denominator):
        """ divide array by denominator = constant or array"""
        return self.data / denominator

    def __sub__(self, sub):
        """ subtract sub from array"""
        #data, dname = fits.Handle_fits.ensure_array(sub)
        return self.data - sub

    def mask_off_order(self, on_order):
        """
       mask out off-order to avoid big numbers there messing up the scaling
       Parameters:
       on_order: array
           array of same size as self.data with off-order masked
       """
        return self.data * on_order

    def interp_shift(self, curve_array, orientation='vertical', pivot='middle'):
        """ shift 2d image array using scipy.interpolate
            this is used to straighten orders using the edge array in curve_array
            usually curve array is found using spectroid
            Parameters:
            --------------
            curve_array : array
                a trace of the order curvature
            orientation: str
                which way to apply shift
            pivot
                origin of the shift, if using sky lines use pivot='peak'
                assuming the peak is the continuum
        """
        orientation = orientation.strip()
        mid = self.data.shape[0] / 2

        shifted = []
        if pivot == 'middle':
            # shift the curve_array to be centered at the middle of the order
            # multiply by -1 to ensure we interpolate.shift the opposite of 
            # the function found 
            shift_curve_array = -1. * (curve_array - mid)

        elif pivot == 'peak':
            # first find the peak 
            if orientation == 'horizontal':
                cross_cut = self.data.sum(axis=0)
            else:
                cross_cut = self.data.sum(axis=1)
            # peak0 = continuum, this will be pivot point -- no change in interp shift
            peak = np.argmax(cross_cut)
            # need final curve_array at peak to be zero
            shift_curve_array = -1 * (curve_array - curve_array[peak])

        else:
            print "ERROR: pivot is not either 'middle' or 'peak' "
            return

        # for each column (each x),  if vertical
        for i in range(0, len(shift_curve_array)):

            if orientation == 'vertical':
                order_slice = self.data[:, i]  # i is x
            else:
                order_slice = self.data[i, :]  # i is y

            # shift order slice to rectify in spatial dimension
            shifted.append(sp.ndimage.interpolation.shift(order_slice,
                                                          shift_curve_array[i],
                                                          order=3,
                                                          mode='nearest',
                                                          prefilter=True))

        rectified = np.array(shifted)
        if orientation == 'vertical': rectified = rectified.transpose()
        # self.rectified = rectified
        self.interp_shifted = True
        return rectified

    def cosmic(self, max_iter=3, sig_clip=5.0, sig_frac=0.3, obj_lim=5.0):
        """ call LA cosmic routine by Malte Tewes"""
        import cosmics  # LA cosmic routine

        reload(cosmics)
        c = cosmics.cosmicsImage(self.data, sigclip=sig_clip, sigfrac=sig_frac,
                                 objlim=obj_lim, verbose=False)
        c.run(max_iter)
        self.data = c.cleanarray
        self.cosmic_cleaned = True

    def find_peak(self, order_shifted = True, use_range="[:]"):
        """ this should be in SciArray
        ---------------
        sciorder: sciArray object with data attribute
        order_shifted: Boolean
            if the order has been shifted or rectified use whole axis

        """
        if order_shifted:
            crosscut = self.data.sum(axis=1)
        else:
            #this isn't right
            crosscut = self.data[:30].sum(axis=1)

        if crosscut.any():
            peak = np.argmax(crosscut)
            location_of_peak = crosscut[peak]
        else:
            peak = 0.0

        return peak, crosscut

    def cut_out(self, padding=30., lower=10., upper=30., orientation='landscape'):
        """ cut out orders -
        adds padding to top and bottom of array
        ensures that bottom-lower is on array
        Parameters:
        -------------
        padding: float
            amount to add to upper location and subract from lower location before
            cutting
        lower: float
            lower location to cut
        upper: float
            upper location to cut
        orientation: str
            if orientation = 'landscape' then upper and lower are in y, otherwise
            they are in x

        Returns:
        -------------
        cut out array
        """

        if orientation != 'landscape':
            data_trans = self.data.transpose()
        else:
            data_trans = self.data

        if lower > padding:
            order_slice = data_trans[lower - padding:upper + padding, :]
        else:
            order_slice = data_trans[0:upper + padding, :]

        if orientation != 'landscape':
            order_slice = order_slice.transpose()

        return order_slice

    def mask_order(self, upper_edge, lower_edge):
        """  use the traces of the orders to make on and off order arrays
         Used primarily to mask out the array data before normalizing or dividing
         Parameters
         ----------
         lower_edge: array that traces the lower edge of the order
         upper_edge: array that traces the upper edge of the order
         Returns
         ----------
         on_order: array with the off-order pixels set to 0.
         off_order: array with the on-order pixels set to 0.
        """

        y, x = np.indices(self.data.shape, dtype=np.float32)

        off_top = y > upper_edge
        off_bot = y < lower_edge
        off_order = off_top | off_bot

        belowtop = y < upper_edge
        abovebot = y > lower_edge
        on_order = belowtop & abovebot

        return on_order, off_order

class SciArray(BaseArray):
    """
    SciArray has all the methods of BaseArray and additional routines that do
    sky line fitting, sky vs continuum ranges, extraction
    """

    def __init__(self, data):
        self.data = data

    def find_skyline_trace(self, sky_sigma=4.0, padding=30.):
        """ this both finds peaks and measures centroids of lines -- should be separate methods"""

        # transpose the array because spectroid can only read horizontal peaks for now
        npsts = self.data.transpose()

        # The order cutout has padding on each side. In order to find the sky lines we should 
        # only look at the central section of the cut out array
        npsts = npsts[:, padding + 5:npsts.shape[1] - 5 - padding]

        cc = np.sum(npsts[:, 0:5], axis=1)

        locpeaks = argrelextrema(cc, np.greater)
        locmaxes = np.where(cc[locpeaks[0]] > sky_sigma * cc.mean())
        maxes = np.array(locpeaks[0][locmaxes[0]])

        deletelist = []
        # remove 'close' real sky line values
        for i in range(1, len(maxes)):
            if abs(maxes[i] - maxes[i - 1]) < 5:
                deletelist.append(i)
        maxes = np.delete(maxes, deletelist, None)

        peaks = cc[maxes]
        sortorder = np.argsort(peaks)
        maxes = maxes[sortorder]
        maxes = maxes[::-1]

        skydict = {}
        centroid_sky_sum = np.array([])
        fitnumber = 0

        for maxskyloc in maxes:
            if 10 < maxskyloc < 1010:

                centroid_sky, badfit = spectroid(npsts, traceWidth=3, backgroundWidth=0, startingLocation=maxskyloc,
                                                 traceMean=True, traceLast=False, traceDelta=0.8)

                if badfit == 9999:
                    continue  # skip this skyline

                # average up the good ones
                if badfit < 10:
                    skydict[fitnumber] = centroid_sky
                    fitnumber += 1
                    # pl.plot(centroid_sky,'k')
                    if centroid_sky_sum.any():
                        centroid_sky_sum = centroid_sky_sum + centroid_sky - centroid_sky[0]
                    else:
                        centroid_sky_sum = centroid_sky - centroid_sky[0]
                if fitnumber > 2:
                    break

        if centroid_sky_sum.any():
            sky_line = centroid_sky_sum / fitnumber
            sky_line_success = True
        else:
            sky_line = NirspecFudgeConstants.badval
            sky_line_success = False

        return sky_line, sky_line_success

    def setup_extraction_ranges(self, ext_height, sky_distance, sky_height,
                                peak, order, logger):

        print "ext_height, sky_distance, sky_height, peak, order=",ext_height, sky_distance, sky_height, peak, order

        if ext_height % 2:  # odd values typically across the continuum - recall range(-1, 2) is (-1, 0, 1)
            ext_range = range(int((1 - ext_height) / 2.), int((ext_height + 1) / 2.))
        else:  # even
            ext_range = range((-ext_height) / 2, ext_height / 2)  # have to make assumption to keep integers
            # range(-2, 2) = (-2, -1, 0, 1) so extra part below continuum

        # sky_distance
        sky_range_top = range(ext_range[-1] + sky_distance,
                              ext_range[-1] + sky_distance + sky_height)  # [6, 7, 8, 9, 10]
        sky_range_bot = range(ext_range[0] - sky_distance - sky_height + 1,
                              ext_range[0] - sky_distance + 1)  # [-10, -9, -8, -7, -6]

        if ((peak + sky_range_bot[-1] + 1) > self.data.shape[0]) or ((peak + sky_range_bot[0]) < 0):
            logger.error(' array_manipulate.setup_extraction ranges: bottom sky range ' + str(sky_range_bot) + ' outside of order')
            logger.error('    trying smaller sky range and distance')
            sky_distance = min(peak - sky_range_bot, self.data.shape[0] - peak)
            sky_height = 2

            if ext_range[0] - sky_distance + 1 < self.data.shape[0]:
                sky_range_bot = range(ext_range[0] - sky_distance - sky_height + 1, ext_range[0] - sky_distance + 1)
            else:
                logger.error(' array_manipulate.setup_extraction ranges: not using a bottom sky subtraction, peak is too close to edge')
                sky_range_bot = []

        if (peak + sky_range_top[-1] + 1) > self.data.shape[0]:
            logger.error(' array_manipulate.setup_extraction ranges: top sky range ' + str(sky_range_top) + ' outside of order')
            logger.error('    trying smaller sky range and distance')
            sky_distance = 1
            sky_height = 4
            if peak + sky_distance + sky_height < self.data.shape[0]:
                sky_range_top = range(ext_range[-1] + sky_distance, ext_range[-1] + sky_distance + sky_height)
            else:
                logger.error(' array_manipulate.setup_extraction ranges: not using a top sky subtraction, peak is too close to edge')
                sky_range_top = []

        if self.data.shape[0] - peak < 2 or peak < 3:
            logger.error('WARNING  array_manipulate.setup_extraction ranges: cannot find a peak in order ' + str(order) + ' skipping extraction ')
            return 'ERROR', 'ERROR', 'ERROR'

        return ext_range, sky_range_top, sky_range_bot


    def get_sky(self, peak, sky_range_bot, sky_range_top):

        sky_height_bot = len(sky_range_bot)
        sky_height_top = len(sky_range_top)

        try:
            sky1 = np.sum(self.data[peak + i, :] for i in sky_range_bot)
        except:
            sky1 = [0]
            sky_height_bot = 0.

        try:
            sky2 = np.sum(self.data[peak + i, :] for i in sky_range_top)
        except:
            sky2 = [0]
            sky_height_top = 0.

        if sky_height_top > 0. and sky_height_bot > 0.:
            sky = (sky1 + sky2) / (sky_height_top + sky_height_bot)
            sky = np.array(sky)
            sky = sky - np.median(sky)
            sky = sky

        else:
            sky = NirspecFudgeConstants.badval

        return sky

    def sum_extract(self, ext_range, sky_range_bot, sky_range_top):
        """ extract peak

        Parameters
        ----------
        ext_range: tuple
          range for extraction of continuum
        sky_range_bot: tuple
          range of where to
        sky_range_top: tuple

        Returns
        -----------
        extracted spectrum
        extracted off peak "sky" spectra
        status
        """
        cross_cut = self.data.sum(axis=1)

        peak = np.argmax(cross_cut)
        data_array = self.data + 0.00000000001  # avoid div by 0
        try:
            extracted = np.sum(data_array[peak - i, :] for i in ext_range)
            extracted /= len(ext_range)
            # sets self.sky
            sky = self.get_sky(peak, sky_range_bot, sky_range_top)
            cont = extracted - sky  # subtract sky from peak
            sky = sky - sky.mean()
            extract_status = 1
        except:
            cont = []
            sky = 'bad'
            extract_status = 0

        return cont, sky, extract_status

class FlatArray(BaseArray):
    """
    FlatArray has all the methods of BaseArray and additional routines that do
    normalization of an order with on/off order masking, shifting and subtracting
    flat to create array of just tops and just bottoms of orders, and traces
    along each top and bottom of order to find the order curvature and location
    """

    def __init__(self, data):
        """


        :type self: object
        :rtype : object
        """
        self.data = data
        # shift and subtract flat from itself 
        # to get top and bottom starting locations
        # sets self.tops, self.bots

        self.make_tops_bots()

    def normalize(self, on_order=1.0, off_order=1.0, mask=True, instr="NIRSPEC"):
        """
        Normalize self.data
        Parameters:
        --------------
        on_order: array
           pixel array of same size of self.data with off-order pixels masked
        off_order: array
           pixel array of same size of self.data with off-order pixels masked
        mask: bool
           when mask=True, off-order pixels are not used in computation of mean
        instr: str
           if instr="NIRSPEC" certain pixels in the flat are not used

        """

        if instr == "NIRSPEC":
            # flat blows up past 1000
            self.data[:, 1000:] = 1.0

        if mask:
            # turn off-order pixels to zero before finding mean
            data_on = self.data * on_order

            # take mean of the only-on order
            mx = np.ma.masked_array(self.data, mask=off_order)

            flat_mean = mx.mean()

        else:
            flat_mean = self.data.mean()
            data_on = self.data

        # create normalized flat
        normalized = data_on / flat_mean

        # around the edges of the order can blow up when div by mean, set those to one
        normalized[np.where(normalized > 10.0)] = 1.0

        # avoid zeroes
        normalized[np.where(normalized == 0.0)] = 1.0
        normalized[np.where(normalized < 0.2)] = 1.0

        return normalized, flat_mean
        #self.normalized = normalized
        #self.flat_mean = flat_mean

    def make_tops_bots(self):
        """
        make arrays with only the top edges of the orders and only the bottom
        edges of the orders. typically a flat is just shifted and subtracted
        from itself.

        """
        f2 = np.roll(self.data, 5, axis=0)
        tops = f2 - self.data
        bots = self.data - f2
        self.tops = tops
        self.bots = bots

    @staticmethod
    def trace_order(tops, bots, nh, data_dict, order, logger, sciobj):
        # moved to trace_order.py, this version is not used
        # this finds theoretical order position, finds the peaks in tops,bots,
        # runs spectroid on those peaks, dx is theoretical wavelength array

        lhs_top_theory, lhs_bot_theory, dx = nh.get_theory_order_pos(order)

        lhs_bot = -1
        lhs_top = -1

        cb = np.array(0)
        ct = np.array(0)
        cm = np.array(0)
        highest_top = -1

        logger.info('Theory -- left bot = ' + str(int(lhs_bot_theory)) + \
                    ' left top = ' + str(int(lhs_top_theory)))

        if lhs_top_theory < sciobj.data.shape[0] + 20 and lhs_bot_theory > -20:

            # Find the peaks in the shifted/subracted flat file near the theoretical peaks
            lhs_top = nh.get_actual_order_pos(tops, lhs_top_theory,
                                              data_dict['threshold'])

            logger.info('found a peak at ' + str(lhs_top))

            # Ensure the theoretical location and the actual location of 
            # the top of the order are close enough
            found_top = astro_math.actual_to_theory(lhs_top, lhs_top_theory,
                                                    data_dict['order_threshold'])

            if not found_top:
                logger.info('searching for top edge on a fainter flat')
                logger.info('  ---> this might affect the rectification')

                lhs_top = nh.get_actual_order_pos(tops, lhs_top_theory,
                                                  data_dict['threshold'] / 3)
                found_top = astro_math.actual_to_theory(lhs_top, lhs_top_theory,
                                                        data_dict['order_threshold'])

                if not found_top:
                    logger.info('Flat is too faint in this order, ' + str(order))
                    logger.info('Cannot find order location')
                    # order = order - 1
                    return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx

            # Find the peaks in the shifted/subracted flat file near the 
            # theoretical peaks
            lhs_bot = nh.get_actual_order_pos(bots, lhs_bot_theory,
                                              data_dict['threshold'])

            logger.info('found a peak at ' + str(lhs_bot))

            # Ensure the theoretical location and the actual location of 
            # the bot of the order are close enough
            found_bot = astro_math.actual_to_theory(lhs_bot, lhs_bot_theory,
                                                    data_dict['order_threshold'])

            if not found_bot:
                logger.info('searching for top edge on a fainter flat')
                logger.info('  ---> this might affect the rectification')

                lhs_bot = nh.get_actual_order_pos(bots, lhs_bot_theory,
                                                  data_dict['threshold'] / 3)
                found_bot = astro_math.actual_to_theory(lhs_bot, lhs_bot_theory,
                                                        data_dict['order_threshold'])

                if not found_bot:
                    logger.info('Flat is too faint in this order ' + str(order))
                    logger.info('Cannot find order location')
                    order -= 1
                    return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx

            logger.info('Measured -- left bot = ' + str(int(lhs_bot)))
            logger.info('            left top = ' + str(int(lhs_top)))

            if lhs_top < lhs_bot:
                logger.info('top of order cannot be below bottom')
                logger.info('skipping order: ' + str(order))
                # order = order-1
                return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx
        else:

            if lhs_bot_theory > 1100 or lhs_top_theory > 1200:
                lhs_top = lhs_top_theory
                lhs_bot = lhs_bot_theory
                # order = order-100
                # else: order = order-1
                # lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm,
                # self.neworder, highest_top, sciobj.dx
            return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx

        # # spectroid on top and bottom of order location #####
        # call flatobj.centroiding task, using the location of the peak and zero
        # background subtraction            
        # ct is the centroid fit to the top of the order
        # bft is bad fit top = number of times spectroid had to self-correct
        ct, bft = spectroid(tops, traceWidth=data_dict['spw'], backgroundWidth=0, startingLocation=lhs_top,
                            traceMean=True, traceLast=False, traceDelta=data_dict['trace_delta'])

        logger.info('had to self correct on top = ' + str(bft) + ' times ')
        if bft > 300.:
            traced_top = False
        else:
            traced_top = True

        try:  # where to cut out the array to isolate order
            highest_top = max(ct[0], ct[1010])
        except:
            traced_top = False

        cb, bfb = spectroid(bots, traceWidth=data_dict['spw'], backgroundWidth=0, startingLocation=lhs_bot,
                            traceMean=True, traceLast=False,
                            traceDelta=data_dict['trace_delta'])

        logger.info('had to self correct on bottom = ' + str(bfb) + ' times ')
        if bfb > 300.:
            traced_bot = False
        else:
            traced_bot = True

        try:
            lhs_bot = cb[1]
        except:
            traced_bot = False

        if not traced_bot or not traced_top:
            logger.info('order: ' + str(order) + ' not on detector  ')
            # don't keep tracing until order=1 once off of chip
            # if lhs_bot_theory > 1100 or lhs_top_theory > 1200:
            # order = order-100
            # else: order = order-1
            return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx

        cm = np.array([])

        if traced_top and traced_bot:
            # find the mean of the centroid along the top of the order and
            # along the bottom of the order
            cm = (ct + cb) / 2.
        elif traced_top:
            logger.info('using only top curve ')
            cm = ct - ((lhs_top - lhs_bot) / 2) + 1.
        elif traced_bot:
            logger.info('using only bottom curve ')
            cm = cb + ((lhs_top - lhs_bot) / 2) + 1.

        return lhs_bot, lhs_top, lhs_bot_theory, lhs_top_theory, cb, ct, cm, order, highest_top, dx



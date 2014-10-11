# -*- coding: utf-8 -*-
"""
array manipulation utilities

@author: jholt 20140315
"""
import fits
import numpy as np
import scipy as sp
import astro_math

reload(astro_math)
try:
    from scipy.signal import argrelextrema
except:
    print 'need to update scipy to get argrelextrema'
from spectroid import spectroid
from fudge_constants import Nirspec_fudge_constants


class BaseArray(object):
    """
    Base class for all science and flat array objects
    :param data - numpy array
    """

    def __init__(self, data):
        self.data = data

    def __mult__(self, mult):
        """ multiply array by mult = constant or mult = array """
        return self.data * mult

    def __div__(self, denominator):
        """ divide array by denominator = constant or array"""
        return self.data / denominator

    def __sub__(self, sub):
        """ subtract sub from array"""
        self.data, foo = fits.Handle_fits.ensure_array(sub)
        return self.data - sub

    def mask_off_order(self, on_order):
        """
       mask out off-order to avoid big numbers there messing up the scaling
       :param: on_order: array
           array of same size as self.data with off-order masked
       """
        self.data_masked = self.data * on_order
        # return data_masked

    def interp_shift(self, curve_array, orientation='vertical', pivot='middle'):
        """ shift array using scipy.interpolate 
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
            self.rectified = self.data
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
        if orientation == 'vertical':
            rectified = rectified.transpose()
        self.rectified = rectified

    def cosmic(self):
        """ call LA cosmic routine by Malte Tewes"""
        import cosmics  # LA cosmic routine

        reload(cosmics)
        c = cosmics.cosmicsimage(self.data, sigclip=self.sig_clip, sigfrac=self.sig_frac,
                                 objlim=self.obj_lim, verbose=False)
        c.run(self.max_iter)
        self.data = c.cleanarray

    def find_peak(self, order_shifted=False):
        """ this should be in SciArray 
        Parameters:
        ---------------
        order_shifted: Boolean
            if the order has been shifted or rectified use whole axis
        
        """

        if order_shifted:
            crosscut = self.data.sum(axis=1)
        else:
            crosscut = self.data[:30]

        if crosscut.any():
            peak = np.argmax(crosscut)
        else:
            peak = 0.0
        self.peak = peak
        self.crosscut = crosscut

    def cut_out(self, padding=30., lower=10., upper=30., orientation='landscape'):
        """ cut out orders -        
        adds padding to top and bottom of array
        ensures that bottom-lower is on array
        Parameters:
        -------------
        :param
            padding: float
                amount to add to upper location and subract from lower location before
                cutting
        :param
            lower: float
                lower location to cut
        :param
            upper: float
                upper location to cut
        :param
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

        self.order_slice = order_slice

    def mask_order(self, upper_edge, lower_edge):
        """  use the traces of the orders to make on and off order arrays  
         Used primarily to mask out the array data before normalizing or dividing   
         Parameters
         ----------
         :param
             lower_edge: array that traces the lower edge of the order
         :param
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

        below_top = y < upper_edge
        above_bot = y > lower_edge
        on_order = below_top & above_bot

        self.on_order = on_order
        self.off_order = off_order


class SciArray(BaseArray):
    """
    SciArray has all the methods of BaseArray and additional routines that do
    sky line fitting, sky vs continuum ranges, extraction
    """

    def __init__(self, data):
        """

        :param data: 
        """
        self.data = data

    def find_skyline_trace(self, sky_sigma=4.0, padding=30.):
        """
        this both finds peaks and measures centroids of lines -- should be separate methods
        :param
            sky_sigma
                number of sigma above the mean to look for sky lines
        :param
            padding
                the amount of padding added to the order, the higher the number, the farther in the array
                the routine looks for sky lines
        """
        centroid_sky = []

        # transpose the array because spectroid can only read horizontal peaks for now
        data_transposed = self.data.transpose()
        # The order cutout has padding on each side. In order to find the sky lines we should 
        # only look at the central section of the cut out array
        data_central = data_transposed[:, padding + 5:data_transposed.shape[1] - 5 - padding]

        cc = np.sum(data_central[:, 0:5], axis=1)

        loc_peaks = argrelextrema(cc, np.greater)
        loc_maxes = np.where(cc[loc_peaks[0]] > sky_sigma * cc.mean())
        maxes = np.array(loc_peaks[0][loc_maxes[0]])

        delete_list = []
        # remove 'close' real sky line values
        for i in range(1, len(maxes)):
            if abs(maxes[i] - maxes[i - 1]) < 5:
                delete_list.append(i)
        maxes = np.delete(maxes, delete_list, None)

        peaks = cc[maxes]
        sorted_order = np.argsort(peaks)
        maxes = maxes[sorted_order]
        maxes = maxes[::-1]

        sky_dict = {}
        centroid_sky_sum = np.array([])
        fit_number = 0

        for max_sky_loc in maxes:
            if 10 < max_sky_loc < 1010:
                # pl.plot(1,max_sky_loc,'r*')
                # pl.plot(centroid_sky,'c')
                # average up the good ones
                centroid_sky, bad_fit_num = spectroid(data_central, dloc=max_sky_loc, spw=3, bkw=0, trace_mean=True,
                                                      trace_last=False, trace_delta=0.8)

                if bad_fit_num == 9999:
                    continue  # skip this skyline

                if bad_fit_num < 10:
                    sky_dict[fit_number] = centroid_sky
                    fit_number += 1
                    # pl.plot(centroid_sky,'k')
                    if centroid_sky_sum.any():
                        centroid_sky_sum = centroid_sky_sum + centroid_sky - centroid_sky[0]
                    else:
                        centroid_sky_sum = centroid_sky - centroid_sky[0]
                    if fit_number > 2:
                        break

        if centroid_sky_sum.any():
            sky_line = centroid_sky_sum / fit_number
            self.sky_line_success = True
        else:
            sky_line = Nirspec_fudge_constants.badval
            self.sky_line_success = False

        self.sky_line = sky_line

    def setup_extraction_ranges(self, ext_height, sky_distance, sky_height,
                                peak, order_num, logger):
        """

        :param ext_height: The height (pixels) of the sky location
        :param sky_distance: The distance (pixels) away from the center of continuum to look for sky
        :param sky_height: The height (pixels) of the sky
        :param peak: The location of the peak (hopefully continuum)
        :param order_num: current order number, used only for logging purposes
        :param logger: logging instance, created using nirspec_util.setup_logger()
        :return: sets self.ext_range, self.sky_range_bot, self.sky_range_top
        """

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
            logger.error(' bottom sky range ' + str(sky_range_bot) + ' outside of order')
            logger.error('    trying smaller sky range and distance')
            sky_distance = min(peak - sky_range_bot, self.data.shape[0] - peak)
            sky_height = 2

            if ext_range[0] - sky_distance + 1 < self.data.shape[0]:
                sky_range_bot = range(ext_range[0] - sky_distance - sky_height + 1, ext_range[0] - sky_distance + 1)
            else:
                logger.error('not using a bottom sky subtraction, peak is too close to edge')
                sky_range_bot = []
                #sky_height /= 2.

        if (peak + sky_range_top[-1] + 1) > self.data.shape[0]:
            logger.error('top sky range ' + str(sky_range_top) + ' outside of order')
            logger.error('    trying smaller sky range and distance')
            sky_distance = 1
            sky_height = 4
            if peak + sky_distance + sky_height < self.data.shape[0]:
                sky_range_top = range(ext_range[-1] + sky_distance, ext_range[-1] + sky_distance + sky_height)
            else:
                logger.error('not using a top sky subtraction, peak is too close to edge')
                sky_range_top = []
                # sky_height = sky_height / 2.

        if self.data.shape[0] - peak < 2 or peak < 3:
            logger.error('WARNING cannot find a peak in order ' + str(order_num) + ' skipping extraction ')
            return 'ERROR', 'ERROR', 'ERROR'

        self.ext_range = ext_range
        self.sky_range_bot = sky_range_bot
        self.sky_range_top = sky_range_top

    def get_sky(self, peak, sky_range_bot, sky_range_top):
        """
        
        :param peak: 
        :param sky_range_bot: 
        :param sky_range_top: 
        :return:
        """
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
            self.sky = sky

        else:
            self.sky = Nirspec_fudge_constants.badval

    def sum_extract(self, ext_range, sky_range_bot, sky_range_top):
        """ extract peak 
        
        Parameters
        ----------
        :type self: object
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
            self.get_sky(peak, sky_range_bot, sky_range_top)
            self.continuum = extracted - self.sky  # subtract sky from peak
            self.sky = self.sky - self.sky.mean()
            self.extract_status = 1
        except:
            self.continuum = []
            self.sky = 'bad'
            self.extract_status = 0


class FlatArray(BaseArray):
    """
    FlatArray has all the methods of BaseArray and additional routines that do
    normalization of an order with on/off order masking, shifting and subtracting
    flat to create array of just tops and just bottoms of orders, and traces 
    along each top and bottom of order to find the order curvature and location
    """

    def __init__(self, data):
        self.data = data

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
            # flat blows up past 1000, mask out this part before finding order locations
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

        self.normalized = normalized
        self.flat_mean = flat_mean

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




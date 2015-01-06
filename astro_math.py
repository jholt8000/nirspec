# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 09:52:38 2014

@author: jholt
"""
import itertools

import numpy as np
from numpy import fft
import scipy as sp
from scipy import optimize

from astropy.io import fits

try:
    from scipy.signal import argrelextrema
except:
    print 'need to update scipy to get argrelextrema'

# collection of array math utilities that do not need self.data as in am classes

def fudge_padding(fcn, padding):
    # if order is very curved, add a bit more padding to
    # ensure we do not interpolate into the continuum.
    # ## This should be elsewhere ###
    if abs(fcn[0] - fcn[-1]) > 20.:
        padding += 10.

    if abs(fcn[0] - fcn[-1]) > 40.:
        padding += 10.
    return padding

def shift_order(top_spectroid, avg_spectroid, bot_spectroid, padding):

    if bot_spectroid.any():

        if float(bot_spectroid[0]) > float(padding):
            # shift the order edge trace to be in the correct place in
            # order cut out before rectification
            # centroid top , centroid middle , centroid bottom
            top_spectroid = top_spectroid - bot_spectroid[0] + padding
            avg_spectroid = avg_spectroid - bot_spectroid[0] + padding
            bot_spectroid = bot_spectroid - bot_spectroid[0] + padding

            order_shifted = True
        else:
            order_shifted = False

        return top_spectroid, avg_spectroid, bot_spectroid, order_shifted

def shift_order_back(top_spectroid, avg_spectroid, bot_spectroid, padding, order_shifted, lhs_bot):
    print 'order_shifted=',order_shifted
    if order_shifted:
        # remove the padding and start at lhs_bot to show plots in correct place
        avg_spectroid = avg_spectroid - padding + lhs_bot
        bot_spectroid = bot_spectroid - padding + lhs_bot
        top_spectroid = top_spectroid - padding + lhs_bot

    return top_spectroid, avg_spectroid, bot_spectroid

def conv_ang_to_mu(dx):
    mu_dx = dx / 10000.
    return mu_dx

def median_comb(cslist):
    """ stack up all arrays or fits files in cslist and take the median
    :param cslist:
    Parameters:
    --------------
    cslist : data list or array containing string FITS file names
    Returns:
    -------------
    median combined numpy array     
    """
    data_array = []
    for name in cslist:
        data = fits.Handle_fits.ensure_array(name)
        data_array.append(data)

    # images is a comma separated list of files to be median combined
    all_images = np.array(data_array)
    med_comb = np.median(all_images, axis=0)

    return med_comb

def arg_max_corr(a, b):
    """ Find the maximum of the cross-correlation - includes upsampling
    """
    if len(a.shape) > 1:
        raise ValueError('Needs a 1-dimensional array.')
    length = len(a)
    if not length % 2 == 0:
        raise ValueError('Needs an even length array.')

    if not a.shape == b.shape:
        raise ValueError('The 2 arrays need to be the same shape')

    # Start by finding the coarse discretised arg_max
    coarse_max = np.argmax(np.correlate(a, b, mode='full')) - length + 1

    omega = np.zeros(length)
    omega[0:length / 2] = (2 * np.pi * np.arange(length / 2)) / length
    omega[length / 2 + 1:] = (2 * np.pi *
                              (np.arange(length / 2 + 1, length) - length)) / length

    fft_a = fft.fft(a)

    def correlate_point(tau):
        rotate_vec = np.exp(1j * tau * omega)
        rotate_vec[length / 2] = np.cos(np.pi * tau)

        return np.sum((fft.ifft(fft_a * rotate_vec)).real * b)

    start_arg, end_arg = (float(coarse_max) - 1, float(coarse_max) + 1)

    max_arg = optimize.fminbound(lambda tau: -correlate_point(tau),
                                 start_arg, end_arg)
    # print 'coarse_max=',coarse_max,' max_arg=',max_arg

    return max_arg

def fit_poly(cm, xes='default', deg=4):
    """
    fit a polynomial of degree=degree to array cm (usually output from spectroid)
    
    """
    p0 = np.polyfit(np.arange(len(cm) - 10), cm[:-10], deg=deg)  # end always drops off
    if xes == 'default':
        xes = np.arange(len(cm))
    if deg == 4:
        cmfit = p0[0] * xes ** 4 + p0[1] * xes ** 3 + p0[2] * xes ** 2 + p0[3] * xes + p0[4]
    elif deg == 3:
        cmfit = p0[0] * xes ** 3 + p0[1] * xes ** 2 + p0[2] * xes + p0[3]
    elif deg == 2:
        cmfit = p0[0] * xes ** 2 + p0[1] * xes + p0[2]
    elif deg == 1:
        cmfit = p0[0] * xes + p0[1]
    else:
        return p0
    return cmfit, p0


def actual_to_theory(loc1, loc2, threshold='40'):
    if abs(loc1 - loc2) < threshold and loc1 > 1:
        return True
    else:
        return False



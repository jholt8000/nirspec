# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 09:52:38 2014

@author: jholt
"""
import itertools

import numpy as np
from numpy import fft
from scipy import optimize

from astropy.io import fits

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
        if isinstance(name, str):
            data = fits.getdata(name)
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

    import pylab
    #pylab.figure(2)
    #pylab.clf()
    #pylab.plot(np.correlate(a,b,mode='full'))

    pylab.figure(3)
    pylab.clf()
    pylab.plot(b)
    pylab.plot(a)
    pylab.show()


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

def smooth(x,window_len=11,window='hanning'):
   """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

   input:
       x: the input signal
       window_len: the dimension of the smoothing window; should be an odd integer
       window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
         flat window will produce a moving average smoothing.

   output:
       the smoothed signal

   example:

   t=linspace(-2,2,0.1)
   x=sin(t)+randn(len(t))*0.1
   y=smooth(x)

   see also:

   np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
   scipy.signal.lfilter

   TODO: the window parameter could be the window itself if an array instead of a string
   NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
   """

   if x.ndim != 1:
       raise ValueError, "smooth only accepts 1 dimension arrays."

   if x.size < window_len:
       raise ValueError, "Input vector needs to be bigger than window size."


   if window_len<3:
       return x


   if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
       raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


   s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
   #print(len(s))
   if window == 'flat': #moving average
       w=np.ones(window_len,'d')
   else:
       w=eval('np.'+window+'(window_len)')

   y=np.convolve(w/w.sum(),s,mode='valid')
   return y

def fit_poly(cm, xes='default', deg=4):
    """
    fit a polynomial of degree=degree to array cm (usually output from spectroid)
    
    """
    #cmfit=[]
    p0 = np.polyfit(np.arange(len(cm) - 10), cm[:-10], deg=deg)  # end always drops off

    if xes == 'default':
        xes = np.arange(len(cm))

    cmfit = np.polyval(p0,xes)

    return cmfit, p0


def actual_to_theory(loc1, loc2, threshold='40'):
    if abs(loc1 - loc2) < threshold and loc1 > 1:
        return True
    else:
        return False


def get_actual_order_pos(edges, theory, sigma):

        # older versions of scipy do not have this
        # I only require it when needed
        from scipy.signal import argrelextrema

        # take a vertical cut of edges
        magcrosscut = np.median(edges[:, 40:50], axis=1)

        import pylab as pl

        #pl.figure(7)
        #pl.clf()
        #pl.plot(magcrosscut)

        # find the highest peaks in crosscut, search +/- 15 pixels to narrow down list
        extrema = argrelextrema(magcrosscut, np.greater, order=35)[0]

        # find crosscut values at those extrema
        magcrosscutatextrema = magcrosscut[extrema]
        #pl.plot(extrema, magcrosscutatextrema,'r*')
        # narrow down extrema list to only ones over sigma

        peaks = np.where(magcrosscutatextrema > sigma)

        actualpeaks = extrema[peaks[0]]
        #pl.plot(actualpeaks, magcrosscut[actualpeaks], 'b*')
        #pl.show()

        if actualpeaks.any():
            actual = min((abs(theory - i), i) for i in actualpeaks)[1]
        else:
            return -3000

        return actual

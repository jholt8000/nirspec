# -*- coding: utf-8 -*-
"""
Created on Wed Jul 03 14:35:13 2013

@author: jholt
"""
import numpy as np
import spectroid

reload(spectroid)
# from numpy import fft
# from scipy import optimize
import astro_math
from fudge_constants import Nirspec_fudge_constants as nfc

try:
    from scipy.signal import argrelextrema
except:
    print 'need to update scipy to get argrelextrema'


class LineId(object):
    def __init__(self, low_disp, dx, sky):

        self.low_disp = low_disp
        self.dx = dx
        self.sky = sky

    def read_OH(self):
        """
        read sky line list and create lists of x and y locations of theoretical
        oh line strengths and locations
        sets self.ohx, ohy
        """
        if self.low_disp:
            ohdatfile = nfc.ohdatfile_low_disp
        else:
            ohdatfile = nfc.ohdatfile_high_disp

        f = open(ohdatfile)
        x = f.readlines()
        linesx = []
        linesy = []

        for line in x:
            line1 = line.split(" ")
            if float(line1[1]) > 0:
                linesx.append(float(line1[0]))
                linesy.append(float(line1[1]))

        self.ohx = linesx
        self.ohy = linesy

    def gauss_sky(self, ohx, ohy, sig):
        """create a synthetic sky spectrum using a linelist with postions and intensities
        make gaussians with sigma=sig at each line peak
        Parameters:
        -------------------
        ohx: list or numpy array 
            all x-values (Angstroms) of known sky line locations
        ohy: list or numpy array
            all y-values of known sky line locations
        sig: float
            sigma of gaussians  
        
        Returns:
        -------------------------
        synthetic sky line numpy array 
        """

        x = np.array(ohx)
        y = np.array(ohy)
        all_g = y[0]
        for i in np.arange(0, x.size):
            if y[i] > 0.01:
                g = y[i] * np.exp(-(self.dx - x[i]) ** 2 / (2. * sig ** 2))
                all_g = all_g + g

        self.fake_sky = all_g

    def find_xcorr_shift(self, ohgs):
        """ find shift between synthetic sky spectrum and real sky spectrum"""
        self.sky = self.sky - self.sky.mean()

        ohg = np.array([self.dx, ohgs])  # ohg is a synthetic spectrum of gaussians
        # at skylines in list

        xcorrshift = astro_math.arg_max_corr(ohg[1], self.sky)

        delta_x = (ohg[0][-1] - ohg[0][0]) / float(ohg[0].size)
        self.lambda_shift = xcorrshift * delta_x

    def identify(self, ohx, ohy):
        """
        Matches up max lines in list of real sky data and catalogue of known sky locations and intensities

        :param ohx: list of known line locations from table
        :param ohy: list of known line intensities from table
        :return: self.matchesdx,
            self. matchesohx
            self.matchesohy
            self.bigohx
            self.bigohy
            self.identify_status
            self.matchesidx

        """
        self.identify_status = 0
        debug = False
        dx = np.array(self.dx)

        # if dx.min() < 20500:
        # debug=True
        dy = np.array(self.sky)
        self.bigohx = np.array([])
        self.bigohy = np.array([])
        self.matchesohx = np.array([])
        self.matchesohy = np.array([])
        self.matchesdx = np.array([])
        self.matchesidx = np.array([])

        # only look at the part of sky line list that is around the theory locations
        locohx = np.intersect1d(np.where(ohx < dx[-1] + 20)[0],
                                np.where(ohx > dx[0] - 20)[0])

        ohxsized = np.array(ohx[locohx[0]:locohx[-1]])
        ohysized = np.array(ohy[locohx[0]:locohx[-1]])

        # ignore small lines in sky line list
        bigohy = ohysized[np.where(ohysized > nfc.sky_line_min)]
        bigohx = ohxsized[np.where(ohysized > nfc.sky_line_min)]

        deletelist = []

        # remove 'overlapping' sky lines 
        if bigohy.any():
            for i in range(1, len(bigohy)):
                if abs(bigohx[i] - bigohx[i - 1]) < nfc.sky_overlap_threshold:
                    deletelist.append(i)
            bigohy = np.delete(bigohy, deletelist, None)
            bigohx = np.delete(bigohx, deletelist, None)
        else:
            return np.array([]), np.array([]), np.array([]), bigohx, bigohy, 'nomatch'
        # look for relative maxes in dy (real sky line peak values) 
        if argrelextrema(dy, np.greater)[0].any():
            relx = dx[argrelextrema(dy, np.greater)[0]]
            rely = dy[argrelextrema(dy, np.greater)[0]]
            idx1 = argrelextrema(dy, np.greater)[0]

            # bixdx is the locations (in x) of any sky peak maximums greater than threshold sig
            bigdx = relx[np.where(rely > nfc.sky_threshold * rely.mean())]
            bigidx = idx1[np.where(rely > nfc.sky_threshold * rely.mean())]

        else:
            # couldn't find any relative maxes in sky 
            return np.array([]), np.array([]), np.array([]), bigohx, bigohy, 'nomatch'

        deletelist = []

        # remove 'overlapping' real sky line values
        for i in range(1, len(bigdx)):
            if abs(bigdx[i] - bigdx[i - 1]) < nfc.sky_overlap_threshold:
                deletelist.append(i)

        bigdx = np.delete(bigdx, deletelist, None)
        bigidx = np.delete(bigidx, deletelist, None)

        matchesohx = []
        matchesohy = []
        matchesdx = []
        matchesidx = []

        if bigohx.any() and bigdx.any():

            # search for shifted doublet
            bigdx2 = bigdx
            bigohx2 = bigohx
            bigohy2 = bigohy
            bigidx2 = bigidx
            happened = 0

            for i in range(0, len(bigdx) - 1):
                if bigdx[i + 1] - bigdx[i] < 2:
                    if debug:
                        print bigdx[i], ' and ', bigdx[i + 1], 'possible doublet'
                    # locx is the part of bigohx  within +/- 4 angstrom of the bigdx possible doublet
                    locx = np.intersect1d(np.where(bigohx2 > (bigdx[i] - 4))[0],
                                          np.where(bigohx2 < (bigdx[i + 1] + 4))[0])
                    if debug:
                        print 'locx=', locx

                    if len(locx) > 1:  # there has to be more than two lines within the range for doublet

                        if len(locx) > 2:
                            # found more than 2 possible sky lines to match with doublet
                            # 'happened' is how many doubles already removed from bigohy

                            # yslice is the part of bigohy that corresponds to bigohx (with a 'happened' fix for
                            # removed doublet fails)
                            yslice = np.array(bigohy2[locx[0] - 2 * happened:locx[-1] - 2 * happened + 1])

                            locxfix = np.zeros(2, dtype=np.int)

                            if len(yslice) > 0:
                                locxfix[0] = np.argmax(yslice)  # location of the peak in yslice
                            else:
                                return np.array([]), np.array([]), np.array([]), bigohx, bigohy, 'nomatch'

                            yslice = np.delete(yslice, locxfix[0])  # remove the max from yslice

                            locxfix[1] = np.argmax(
                                yslice)  # find the location of the next max; second biggest in original slice

                            if locxfix[1] <= locxfix[0]:
                                locxfix[1] += 1  # if lowest peak then highest peak

                            locx = locx[locxfix]
                            locx.sort()
                            if debug:
                                print 'locx=', locx

                        ohslice = np.array(bigohx2[locx[0] - 2 * happened:locx[1] - 2 * happened + 1])
                        if debug:
                            print 'ohslice=', ohslice, ' are in the same location as', bigdx[i], bigdx[i + 1]
                        if len(ohslice) > 1:
                            for j in range(0, 1):
                                if debug:
                                    print 'j=', j
                                if (ohslice[j + 1] - ohslice[j]) < 2 and abs(ohslice[j] - bigdx2[i - 2 * happened]) < 6\
                                        and abs(ohslice[j + 1] - bigdx2[i + 1 - 2 * happened]) < 6:
                                    if debug:
                                        print ohslice[j], ohslice[j + 1], 'is same doublet as ', \
                                            bigdx2[i - 2 * happened], bigdx2[i + 1 - 2 * happened]

                                    matchesohx.append(ohslice[j])
                                    matchesohx.append(ohslice[j + 1])
                                    matchesohy.append(bigohy2[locx[0] - 2 * happened + j])
                                    matchesohy.append(bigohy2[locx[0] - 2 * happened + j + 1])
                                    matchesdx.append(bigdx2[i - 2 * happened])
                                    matchesdx.append(bigdx2[i - 2 * happened + 1])
                                    matchesidx.append(bigidx[i - 2 * happened])
                                    matchesidx.append(bigidx[i - 2 * happened + 1])

                                    if debug:
                                        print 'removing bigdxs', bigdx2[i - 2 * happened], bigdx2[i - 2 * happened + 1]
                                        print 'removing bigoxs', bigohx2[locx[0] - 2 * happened + j], \
                                            bigohx2[locx[0] - 2 * happened + j + 1]
                                        print 'before dx2=', bigdx2
                                        print 'before oh2=', bigohx2

                                    bigdx2 = np.delete(bigdx2, i - 2 * happened)
                                    bigdx2 = np.delete(bigdx2, i - 2 * happened)  # this removes the "i+1"
                                    bigohx2 = np.delete(bigohx2, locx[0] - 2 * happened + j)
                                    bigohx2 = np.delete(bigohx2, locx[0] - 2 * happened + j)  # this removes the "j+1"
                                    bigohy2 = np.delete(bigohy2, locx[0] - 2 * happened + j)
                                    bigohy2 = np.delete(bigohy2, locx[0] - 2 * happened + j)  # this removes the "j+1"
                                    bigidx2 = np.delete(bigidx2, i - 2 * happened)
                                    bigidx2 = np.delete(bigidx2, i - 2 * happened)

                                    happened += 1

            bigdx = bigdx2
            bigidx = bigidx2
            bigohx = bigohx2
            bigohy = bigohy2

            if debug:
                print 'bigohx=', bigohx
                print 'bigohy=', bigohy
                print 'bigdx=', bigdx
                print 'bigidx=', bigidx

            for j in range(0, len(bigohx)):
                minimum = min((abs(bigohx[j] - i), i) for i in bigdx)

                if (minimum[0]) < 4.0:
                    matchesohx.append(bigohx[j])
                    matchesohy.append(bigohy[j])
                    matchesdx.append(minimum[1])

                    for idx in range(len(bigidx)):
                        if bigdx[idx] == minimum[1]:
                            matchesidx.append(bigidx[idx])

        if len(matchesdx) > 2:
            if debug:
                print 'matchesdx:', matchesdx
                print 'matchesohx:', matchesohx
                print 'matchesohy:', matchesohy
                print 'matchesidx:', matchesidx
            # check for duplicates
            matchesdx2 = matchesdx[:]
            matchesohx2 = matchesohx[:]
            matchesohy2 = matchesohy[:]
            matchesidx2 = matchesidx[:]
            happened = 0
            for j in range(0, len(matchesdx) - 1):

                if abs(matchesdx[j + 1] - matchesdx[j]) < 0.01:
                    if debug:
                        print 'duplicate=', matchesdx[j + 1], matchesdx[j]
                    # find which oh does it actually belongs to 

                    if min(matchesdx[j + 1] - matchesohx[j + 1], matchesdx[j + 1] - matchesohx[j]) == 0:
                        matchesdx2.pop(j + 1 - happened)
                        matchesidx2.pop(j + 1 - happened)
                        matchesohx2.pop(j + 1 - happened)
                        matchesohy2.pop(j + 1 - happened)
                    else:
                        matchesdx2.pop(j - happened)
                        matchesidx2.pop(j - happened)
                        matchesohx2.pop(j - happened)
                        matchesohy2.pop(j - happened)

                    happened += 1

            matchesdx = np.array(matchesdx2)
            matchesohx = np.array(matchesohx2)
            matchesohy = np.array(matchesohy2)
            matchesidx = np.array(matchesidx2)

            matchesdx.sort()
            matchesidx.sort()

            oh_sort_indices = matchesohx.argsort()
            matchesohy = matchesohy[oh_sort_indices]
            # matchesohx.sort()
            matchesohx = matchesohx[oh_sort_indices]

            self.matchesdx = matchesdx
            self.matchesohx = matchesohx
            self.matchesohy = matchesohy
            self.bigohx = bigohx
            self.bigohy = bigohy
            self.identify_status = 1
            self.matchesidx = matchesidx

def sanity_check(orig_pix_x, order_number_array, matched_sky_line):
    """
    tries to fit a line to each ID/OH value, throws out bad fits
    :param orig_pix_x:
    :param order_number_array:
    :param matched_sky_line:
    :return:
    """

    i=0
    while i < len(order_number_array):
        f1, residuals1, rank1, singular_values1, rcond1 = np.polyfit(orig_pix_x[i], matched_sky_line[i], 1, full=True)
        f2, residuals2, rank2, singular_values2, rcond2 = np.polyfit(orig_pix_x[i], matched_sky_line[i], 2, full=True)

        if float(residuals2) > 500:
            print 'order number ',order_number_array[i][0],' is a bad fit'
            orig_pix_x.pop(i)
            order_number_array.pop(i)
            matched_sky_line.pop(i)
        i+=1

    return orig_pix_x, order_number_array, matched_sky_linex



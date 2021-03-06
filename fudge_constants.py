# -*- coding: utf-8 -*-
"""
Created on Tue Jul 01 10:43:26 2014

@author: jholt
"""


class NirspecFudgeConstants(object):
    lambda_linear_fit_threshold = 0.05
    disp_lower_limit = 0.95
    disp_upper_limit = 1.05
    badval = 9999
    max_shift_from_theory = 50.  # pixels
    ohdatfile_high_disp = '/net/filer2/vol/home/jholt/nirspec/nirspec_2015/ir_ohlines.dat'
    #ohdatfile_high_disp = 'ir_ohlines.dat'
    # used in nirspec_wavelength_utils identify
    sky_threshold = 3.
    sky_overlap_threshold = 0.6
    sky_line_min = 10.
    # if theory if total_chip_padding  
    total_chip_padding = 20.
    badfit_limit = 200.

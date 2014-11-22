# -*- coding: utf-8 -*-
"""
Spectroid finds a trace of the centroid of a 2 dimensional data array in the 
x-direction. It can subtract a background gradient value and self correct if 
centroid locations have too high of a slope. It returns the centroid array and
the counter of how many self corrections occurred. The second value is a 
reasonable measure of the goodness of the calculation for most smooth curves. 

Created on Thu Mar 28 10:30:41 2013

@author: jholt

2014/03/02 improve background logic
2014/07/25 use scipy center-of-mass routine instead of handmade algorithm
"""
import numpy as np
from scipy import ndimage

def spectroid(data_array, spw=10, bkw=30, dloc=926, trace_mean=True,
              trace_last=False, trace_delta=1.9):
    """ Find centroids for each column in data_array -- only works horizontally
    
     data_array: 2d numpy data array, usually read in from pyfits
     spw=10: is distance away, up and down, from center to search for centroid
     bkw=30: is background distance from center to subtract before finding centroid
     dloc=y: starting location for centroid at x=0
     trace_mean=True: if trace_mean is true and centroid found is 
           >trace_delta*abs(previous centroid found) then use mean of previous 
           3 centroids  
           warning: only use trace_mean if curve to be fit does not have large jumps or 
           discontinuities
     trace_last=False: if trace_last is true and centroid found is 
          >trace_delta*abs(previous centroid found) then use previous centroid  
           warning: only use trace_last if curve to be fit does not have large jumps or 
           discontinuities
           trace_last will take precedence over trace_mean if both are set to True
     trace_delta: threshold for unacceptable jump in centroid location
     
    returns: the centroid trace along the x-direction, the number of times the 
    centroiding went outside of centroid[i-1]*trace_delta
    """

    # initialize "center of mass" array to be populated with centroids
    # cm=np.zeros(data_array.shape[1])
    cm = np.zeros(data_array.shape[1])

    cm2 = []
    #R -= R.sum(0) / len(R)

    # badfit is a measure of how many times routine had to correct for bad centroid
    badfit = 0

    # first centroid will be starting location input parameter
    # cm[0]=dloc
    cm[0] = dloc
    cm2.append(dloc)
    xmin = 1
    xmax = data_array.shape[1]

    for i in range(xmin, xmax):
        # i is each "x" or column
        # lowest height (or row) to search for centroid

        ymin = int(cm[i - 1] - spw)
        # maximum height (or row) to search for centroid
        ymax = int(cm[i - 1] + spw)

        # if the max is bigger than the top, use the top of array
        if abs(ymax) > data_array.shape[0]:
            ymax = int(data_array.shape[0])

        # if the min is smaller than the
        #  bottom, use the bottom of array
        if ymin < 1: #dont let it trace the bottom of the detector
            ymin = 1
        if ymax <= 0:
            ymax = int(cm[i] + spw)+1

        if bkw > 0:  # Subtract background
            # bkmin is the center found for the previous column minus the bkw (background width)
            # location of lower background sample
            bkmin = cm[i - 1] - bkw
            # location of top background sample
            bkmax = cm[i - 1] + bkw

            # do not go outside of the array
            if bkmax > data_array.shape[0]:
                bkmax = data_array.shape[0] - 1
            if bkmin < 0:
                bkmin = 0

            # find the values to subtract from the array 
            bkmax_val = data_array[bkmax, i]
            bkmin_val = data_array[bkmin, i]
            bk_mean = (bkmax_val + bkmin_val) / 2.

        else:
            bk_mean = 0.0  # do not subtract background value

        # centroid eqn for column 1: 
        # (Sum(data_array[j,1]) * j / Sum(data_array[j,1]))

        # for each column (i)  find the sum (centroid_value_j * y_j) 
        # where centroid_value_j is the value at y=j,x=i         

        cm[i] = ndimage.measurements.center_of_mass(data_array[int(ymin):int(ymax)+1, i] - bk_mean)[0] + ymin

        # return data_array[:,i]
        if cm[i] is np.inf or cm[i] is -np.inf : # This happens when we went off array
            return np.array(0.0), 9999

        # centroid jumped more than trace_delta
        if np.abs(cm[i] - cm[i - 1]) > trace_delta :
            badfit += 1
            if trace_mean:
                if i > 4:
                    # jump is past beginning, use past three centroids

                    cm[i] = cm[i-3:i-1].mean()
                elif i > 1:
                    # average as many cms as we have gone through
                    cm[i] = cm[i-2:i-1].mean()
                else:
                    # use the first one found
                    cm[i] = cm[i-1]

            elif trace_last:
                cm[i] = cm[i - 1]

    return cm, badfit



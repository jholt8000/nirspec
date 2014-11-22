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


    print 'type of data array=',type(data_array)
    # initialize "center of mass" array to be populated with centroids
    # cm=np.zeros(data_array.shape[1])
    cm = []
    # badfit is a measure of how many times routine had to correct for bad centroid
    badfit = 0

    # first centroid will be starting location input parameter
    # cm[0]=dloc
    cm.append(dloc)
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

        # if the min is smaller than the bottom, use the bottom of array
        if ymin < 0:
            ymin = 0

        if bkw > 0:  # Subtract background
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
            bk_mean = 0.0  # do not subract background value

        # centroid eqn for column 1: 
        # (Sum(data_array[j,1]) * j / Sum(data_array[j,1]))

        # for each column (i)  find the sum (centroid_value_j * y_j) 
        # where centroid_value_j is the value at y=j,x=i         

        cm_n = sum([(data_array[j, i] - bk_mean) * j for j in range(int(ymin), int(ymax) + 1)])
        cm_d = sum([(data_array[j, i] - bk_mean) for j in range(int(ymin), int(ymax) + 1)])

        # test = data_array[int(ymin):int(ymax),i]
        # cm[i]= ndimage.measurements.center_of_mass(test-bk_mean)[0]+ymin
        # if i < 20: print cm[i], ymin, ymax

        # return data_array[:,i]
        if cm_d == 0:
            # return from function before dividing by zero
            return np.array(0.0), 9999

        cm.append(float(cm_n) / cm_d)

        # centroid jumped more than trace_delta
        if np.abs(cm[i] - cm[i - 1]) > trace_delta:

            badfit += 1

            if trace_mean:

                if len(cm) > 4:
                    # jump is past beginning, use past three centroids
                    kmax = 3
                elif len(cm) > 1:
                    # average as many cms as we have gone through
                    kmax = len(cm)
                else:
                    # use the first one found
                    kmax = 2

                cm[i] = (sum([cm[i - k] for k in range(1, kmax)])) / (float(kmax) - 1)

            elif trace_last:
                cm[i] = cm[i - 1]

    cm = np.array(cm)

    return cm, badfit



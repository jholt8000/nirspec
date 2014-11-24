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


def spectroid(dataArray, traceWidth=10, backgroundWidth=30, startingLocation=926, traceMean=True,
              traceLast=False, traceDelta=1.9):
    """ Find centroids for each column in dataArray -- only works horizontally
    
     :parameter dataArray: 2d numpy data array, usually read in from pyfits
     :parameter traceWidth=10: is distance away, up and down, from center to search for centroid
     :parameter backgroundWidth=30: is background distance from center to subtract before finding centroid
     :parameter  startingLocation=y: starting location for centroid at x=0
     :parameter traceMean=True: if traceMean is true and centroid found is
           >traceDelta*abs(previous centroid found) then use mean of previous 
           3 centroids  
           warning: only use traceMean if curve to be fit does not have large jumps or 
           discontinuities
     :parameter traceLast=False: if traceLast is true and centroid found is
          >traceDelta*abs(previous centroid found) then use previous centroid  
           warning: only use traceLast if curve to be fit does not have large jumps or 
           discontinuities
           traceLast will take precedence over traceMean if both are set to True
     :parameter traceDelta: threshold for unacceptable jump in centroid location
     
    :returns: the centroid trace along the x-direction, the number of times the
    centroiding went outside of centroid[i-1]*traceDelta
    """

    # initialize "center of mass" array to be populated with centroids
    # cm=np.zeros(dataArray.shape[1])
    cm = np.zeros(dataArray.shape[1])

    cm2 = []
    # R -= R.sum(0) / len(R)

    # badfit is a measure of how many times routine had to correct for bad centroid
    badfit = 0

    # first centroid will be starting location input parameter
    # cm[0]=startingLocation
    cm[0] = startingLocation
    cm2.append(startingLocation)
    xmin = 1
    xmax = dataArray.shape[1]

    for i in range(xmin, xmax):
        # i is each "x" or column
        # lowest height (or row) to search for centroid

        ymin = int(cm[i - 1] - traceWidth)
        # maximum height (or row) to search for centroid
        ymax = int(cm[i - 1] + traceWidth)

        # if the max is bigger than the top, use the top of array
        if abs(ymax) > dataArray.shape[0]:
            ymax = int(dataArray.shape[0])

        # if the min is smaller than the
        # bottom, use the bottom of array
        if ymin < 1:  # don't let it trace the bottom of the detector
            ymin = 1
        if ymax <= 0:
            ymax = int(cm[i] + traceWidth) + 1

        if backgroundWidth > 0:  # Subtract background
            # backgroundMin is the center found for the previous column minus the backgroundWidth (background width)
            # location of lower background sample
            backgroundMin = cm[i - 1] - backgroundWidth
            # location of top background sample
            backgroundMax = cm[i - 1] + backgroundWidth

            # do not go outside of the array
            if backgroundMax > dataArray.shape[0]:
                backgroundMax = dataArray.shape[0] - 1
            if backgroundMin < 0:
                backgroundMin = 0

            # find the values to subtract from the array 
            backgroundMax_val = dataArray[backgroundMax, i]
            backgroundMin_val = dataArray[backgroundMin, i]
            bk_mean = (backgroundMax_val + backgroundMin_val) / 2.

        else:
            bk_mean = 0.0  # do not subtract background value

        # centroid eqn for column 1: 
        # (Sum(dataArray[j,1]) * j / Sum(dataArray[j,1]))

        # for each column (i)  find the sum (centroid_value_j * y_j) 
        # where centroid_value_j is the value at y=j,x=i         

        cm[i] = ndimage.measurements.center_of_mass(dataArray[int(ymin):int(ymax) + 1, i] - bk_mean)[0] + ymin

        # return dataArray[:,i]
        if cm[i] is np.inf or cm[i] is -np.inf:  # This happens when we went off array
            return np.array(0.0), 9999

        # centroid jumped more than traceDelta
        if np.abs(cm[i] - cm[i - 1]) > traceDelta:
            badfit += 1
            if traceMean:
                if i > 4:
                    # jump is past beginning, use past three centroids
                    cm[i] = cm[i - 3:i - 1].mean()
                elif i > 1:
                    # average as many cms as we have gone through
                    cm[i] = cm[i - 2:i - 1].mean()
                else:
                    # use the first one found
                    cm[i] = cm[i - 1]

            elif traceLast:
                cm[i] = cm[i - 1]

    return cm, badfit



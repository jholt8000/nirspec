# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:41:32 2014

@author: jholt
"""
import sys
import nirspec_util
reload(nirspec_util)
import astro_math
reload(astro_math)
import array_manipulate
reload(array_manipulate)
import nirspec_wavelength_utils
reload(nirspec_wavelength_utils)
import robust
reload(robust)

import trace_order
reload(trace_order)
from fudge_constants import Nirspec_fudge_constants

try: 
    import numpy as np
except: 
    print "ERROR: you do not have pyfits, pylab, numpy, and/or scipy installed!"
    print "       exiting...."
    sys.exit()

# This should be more stand-alone than it is, the input args should be spelled
# out better
        
def reduce_order(reduction, sciobj, flatobj):
    '''
    redue each order found
    Parameters:
    --------------
    reduction: object
        The reduction object from the wrapper nirspecOO, has attributes: 
            reduction.do_extract
            reduction.order_num (should be order_num)
            reduction.logger
            reduction.padding
            
    sciobj: object
        The science array object, contains sciobj.data array
    flatobj: object
        The flat array object, contains flatoj.data
    '''
    ### make class? attributes could be padding and dx as those change throughout ###
    
    order_data=[]
    found_wavelength=False

    # make order-specific extraction bool
    do_extract = reduction.do_extract
    
    # keep this order number (reduction.order_num will change)
    reduction.order_num = reduction.order_num

    reduction.logger.info('-------Starting on Order '+str(reduction.order_num)+'---------\n \n')       
    
    # find order position on detector #### 
    
    # Use's grating eqn empirical fitting to determine the theoretical 
    # starting location of the top and bottom of order as well as 
    # starting wavelengths for that order. 
    #  lhs_bot_theory,lhs_top_theory, dx, lhs_bot, lhs_top, bottom_spectroid,
    #top_spectroid, avg_spectroid, highest_top
    traceobj = trace_order.Trace_order_utils(reduction, flatobj, sciobj)
    traceobj.trace_order()

    # Ensure that order was traced successfully and is on detector    
    if not traceobj.trace_success:    
        reduction.logger.info('order '+str(reduction.order_num)+' not on array')
        # return to main reduction and skip this order
        return order_data, traceobj.lhs_top

    padding = traceobj.padding
    
    ## cut out order and spatially rectify #### 
    if traceobj.trace_success:
 
        reduction.logger.info('cutting out order'+str(reduction.order_num)+' from science and flat')

        ### cut out the order from the science and the flat ###
        ### include padding on the top and bottom to ensure order is on the cutout
        ### array and to avoid cutting into the science when order is straightened
        
        print 'padding=',padding
        print 'lhsbot,lhstop',traceobj.lhs_bot, traceobj.highest_top
        # sets sciobj.order_slice that contains only the cutout around the order
        sciobj.cut_out(padding = padding, lower = traceobj.lhs_bot, 
                       upper = traceobj.highest_top)
         
        # sets flatobj.order_slice that contains only the cutout around the order          
        flatobj.cut_out(padding = padding, lower = traceobj.lhs_bot, 
                        upper = traceobj.highest_top)
                                              
        #make instances of array manip class using just the order slice
        sciorder = array_manipulate.SciArray(sciobj.order_slice)
        flatorder = array_manipulate.FlatArray(flatobj.order_slice)

        # copy new order-specific dx array to be attribute of sciorder
        sciorder.dx = traceobj.dx
        
        reduction.logger.info('masking out off-order locations')
        
        ### normalize flatorder data ###
        
        # specify the locations of actual order and off order
        # sets flatorder.on_order and flatorder.off_order
        flatorder.mask_order(traceobj.top_spectroid, traceobj.bottom_spectroid)
           
        reduction.logger.info('normalizing flat order '+str(reduction.order_num))
        
        # sets flatorder.normalized and flatorder.flat_mean
        import pylab as pl
        pl.clf()
        pl.imshow(flatorder.on_order)
        pl.imshow(flatorder.off_order)
        pl.show()
        flatorder.normalize(flatorder.on_order, flatorder.off_order,
                            mask = True, instr = "NIRSPEC")
        
        reduction.logger.info('flatfielding science')  
        
        # sets sciorder.masked
        sciorder.mask_off_order(flatorder.on_order)
                              
        # Where is this flat corrected science array used? should norm_data be data?
        sciorder.norm_data = sciorder.data / flatorder.normalized

        # rectify the order slice, sets sciorder.rectified        
        sciorder.interp_shift(traceobj.avg_spectroid, orientation = 'vertical', pivot = 'middle')

        #remove the padding and start at lhs_bot to show plots in correct place
        traceobj.shift_order_back()
            
        # overwrite data array as rectified array        
        sciorder.data = sciorder.rectified 
        
        order_rectified=True
 
    else: # could not find fit along flat order edge              
        reduction.logger.info('WARNING: did not find a good order fit, not rectifying spatially ')
        reduction.logger.info('         and therefore not extracting continuum ')
        do_extract = False
        order_rectified = False
        return order_data, traceobj.lhs_top

        
    ## determine sky and continuum locations  #### 
    
    #take a the crosscut of the order and determine place off-peak for sky lines
    # sets sciorder.crosscut and sciorder.peak
    sciorder.find_peak(order_rectified)
        
    # sets sciorder.ext_range, sciorder.sky_range_bot, sciorder.sky_range_top      
    sciorder.setup_extraction_ranges(reduction.ext_height, reduction.sky_distance, 
                                         reduction.sky_height, sciorder.peak, reduction.order_num, 
                                         reduction.logger)
           
    ### Horizontal order rectification ###
    
    # sets sciorder.sky_line to use in rectification using sky lines       
    sciorder.find_skyline_trace(sky_sigma = reduction.sky_sigma, 
                                           padding=padding)

    try: 
        sky_line_fit, foo = astro_math.fit_poly(sciorder.sky_line, 
                                     xes=np.arange(sciorder.data.shape[0]),
                                     deg=1)
        sky_line_fit_success = True
                            
    except:
        sky_line_fit_success = False
        
        
    if sky_line_fit_success: # or len(sky_line_fit[20:]) > 20:
         ## rectify spectral direction using sky line centroid fits ####        
         # sets sciorder.rectified
        sciorder.interp_shift(sky_line_fit, orientation = 'horizonal',
                              pivot = 'peak')  
        
        sciorder.data = sciorder.rectified
        
    else: 
        text = 'WARNING could not fit sky lines, not rectifying spectral dimension '
        reduction.logger.info(text)

    # Extract spectrum and subtract sky ## 
    #Gain = 4 e- / ADU  and readnoise = 625 e- = 156 ADU
    if do_extract:      
        
        # sets sciorder.cont, sciorder.skys, sciorder.extract_status        
        sciorder.sum_extract(sciorder.ext_range, sciorder.sky_range_bot, 
                             sciorder.sky_range_top)
                             
        if sciorder.extract_status ==  0:
            reduction.logger.error('could not extract order '+str(reduction.order_num))

            return order_data, traceobj.lhs_top

            
        # identify sky lines with catalogue sky line locations

        #### Wavelength skyline identification ######## 
        lineobj = nirspec_wavelength_utils.Line_id(reduction.low_disp,
                                                      sciorder.dx, 
                                                      sciorder.skys)
        #Test commit
        ## find and apply wavelength shift ###
        
        # Read in the sky line list. sets lineobj.ohx, lineobj.ohy
        # skyline 
        # list is determined using low_disp and fudge_constants ohlinelist
        lineobj.read_OH()
        
        # make a synthetic sky spectrum using line list information with width
        # the size of the data and with sigma = 0.2 (found empirically)
        # sets lineobj.fake_sky
        lineobj.gaussOH(lineobj.ohx, lineobj.ohy, 0.2)
                                                 
        # Cross correlate data: lineobj.dx and lineobj.fake_sky
        # with the synthetic sky, makes lineobj.lambda_shift 
        lineobj.find_xcorr_shift(lineobj.fake_sky)
    
        reduction.logger.info('shift between sky list and real sky lines = '+str(int(lineobj.lambda_shift))) 
        
        if abs(lineobj.lambda_shift) < Nirspec_fudge_constants.max_shift_from_theory:
            sciorder.dx = sciorder.dx + lineobj.lambda_shift     
            lineobj.dx = sciorder.dx # i hate this, identify should have dx as a parameter -JH
            reduction.logger.info('applied the shift = '+str(lineobj.lambda_shift)+' to aid in sky line identification')      
            
        # match sky lines
        # sets lineobj.matchesdx, lineobj.matchesohx, lineobj.matchesohy, 
        # lineobj.bigohx, lineobj.bigohy, lineobj.identify_status 
        lineobj.identify(lineobj.ohx, lineobj.ohy)
            
        if lineobj.identify_status < 1:
            reduction.logger.info('problem with sky line identification')

            if abs(lineobj.lambda_shift) < Nirspec_fudge_constants.max_shift_from_theory:
                sciorder.dx = sciorder.dx-lineobj.lambda_shift      
                lineobj.dx = sciorder.dx
                reduction.logger.info('removed the shift since sky is unreliable : '+str(-lineobj.lambda_shift))   
                reduction.logger.error('Could not find sky lines: only doing zeroith order wavelength calibration ')
                found_wavelength=False
        
            
        else:
            # find the solution between current zeroith order solution and real lambda
            #foo, (disp,offset) = astro_math.fit_poly(lineobj.matchesohx, xes=lineobj.matchesdx, deg=1, cutoff=False)
            p0=np.polyfit(lineobj.matchesohx,lineobj.matchesdx, deg=1) 

            disp=p0[0]
            offset = p0[1]
            
            reduction.logger.info('linear solution bw matches disp='+str(disp)+' offset ='+str(offset))
            
            # if the matches between theory and sky line list are too far off from a linear fit, dont use matches
            if abs(disp) < Nirspec_fudge_constants.disp_upper_limit and abs(disp) > Nirspec_fudge_constants.disp_lower_limit: 
                found_wavelength=True
                #astro_math.order_wavelength_solution(lineobj.matchesdx, lineobj.matchesohx, sciorder.dx, order)

            else:
                reduction.logger.error('bad fit: only doing zeroith order wavelength calibration ')
                found_wavelength=False


        #need to store all the pre-wavelength fixed data to make out files

        order_data=([reduction.order_num, sciorder, lineobj, flatobj, traceobj, found_wavelength])


                            
    return order_data, traceobj.lhs_top

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 01 15:11:52 2014

@author: jholt
"""
from fudge_constants import Nirspec_fudge_constants as nfc
#reload(fudge_constants)
import fudge_constants
reload(fudge_constants)
from spectroid import spectroid
import numpy as np
import astro_math

class Trace_order_utils(object):
    
    def __init__(self, reductionobj, flatobj, sciobj):
        self.reductionobj = reductionobj
        self.flatobj = flatobj
        self.sciobj = sciobj
        self.trace_success = False
    
    def trace_order(self):

        self.lhs_bot_theory = -1
        self.lhs_top_theory = -1
        self.dx = []

        self.lhs_bot = -1
        self.lhs_top = -1
 
        self.cb = np.array(0)
        self.ct = np.array(0)
        self.cm = np.array(0)
        self.highest_top=-1 
        
        self.trace_success = False
        
        # this finds theoretical order position, finds the peaks in tops,bots, 
        # runs spectroid on those peaks, dx is theoretical wavelength array
        
        self.lhs_top_theory, self.lhs_bot_theory, self.dx = self.reductionobj.nh.get_theory_order_pos(self.reductionobj.order)
   
        self.reductionobj.logger.info('Theory -- left bot = '+str(int(self.lhs_bot_theory))+\
                    ' left top = '+str(int(self.lhs_top_theory)))
        
        if self.lhs_top_theory < self.sciobj.data.shape[0]+fudge_constants.Nirspec_fudge_constants.total_chip_padding and self.lhs_bot_theory > -fudge_constants.Nirspec_fudge_constants.total_chip_padding: 
             
            self.lhs_top = self.determine_lhs_edge_pos(self.flatobj.tops, self.lhs_top_theory, self.reductionobj.data_dict)                               
               
            self.lhs_bot = self.determine_lhs_edge_pos(self.flatobj.bots, self.lhs_bot_theory, self.reductionobj.data_dict)

            self.reductionobj.logger.info('Measured -- left bot = '+str(int(self.lhs_bot)))
            self.reductionobj.logger.info('            left top = '+str(int(self.lhs_top)))
            
            if self.lhs_top > self.lhs_bot:
         
                ## spectroid on top and bottom of order location #####
                #call centroiding task, using the location of the peak and zero 
                # background subtraction            
                # ct is the centroid fit to the top of the order
                # bft is bad fit top = number of times spectroid had to self-correct
                self.ct, bft = spectroid(self.flatobj.tops, dloc = self.lhs_top, spw = self.reductionobj.data_dict['spw'], 
                                    bkw = 0, trace_mean = True, trace_last = False, 
                                    trace_delta = self.reductionobj.data_dict['trace_delta']) 
                                    
                self.reductionobj.logger.info('had to self correct on top = '+str(bft)+' times ')
                if bft > fudge_constants.Nirspec_fudge_constants.badfit_limit: self.traced_top = False
                else: self.traced_top = True
               
                try:         #where to cut out the array to isolate order
                    self.highest_top = max(self.ct[0], self.ct[1010]) 
                except: self.traced_top = False    
                    
                self.cb, bfb = spectroid(self.flatobj.bots, dloc = self.lhs_bot, spw = self.reductionobj.data_dict['spw'], 
                                            bkw = 0, trace_mean = True, trace_last = False, 
                                            trace_delta = self.reductionobj.data_dict['trace_delta'])
                                            
                self.reductionobj.logger.info('had to self correct on bottom = '+str(bfb)+' times ')
                if bfb > fudge_constants.Nirspec_fudge_constants.badfit_limit: self.traced_bot = False
                else: self.traced_bot = True
                
                try: self.lhs_bot = self.cb[1]
                except: self.traced_bot = False               

                if self.traced_top and self.traced_bot:
                        #find the mean of the centroid along the top of the order and 
                        # along the bottom of the order
                        self.cm = (self.ct + self.cb) / 2.        

                elif self.traced_top:  
                    self.reductionobj.logger.info('using only top curve ')
                    self.cm = self.ct-((self.lhs_top-self.lhs_bot)/2)+1.
                elif self.traced_bot:
                    self.reductionobj.logger.info('using only bottom curve ')
                    self.cm = self.cb+((self.lhs_top-self.lhs_bot)/2)+1.
                    
                else:
                    self.reductionobj.logger.info('could not trace order '+str(self.reductionobj.order))
                    self.trace_success = False
                    self.lhs_top = self.lhs_top_theory
                    
     
            else:
                self.reductionobj.logger.info('top of order cannot be below bottom')
                self.reductionobj.logger.info ('skipping order: '+str(self.reductionobj.order))
                self.trace_success = False
                self.lhs_top = self.lhs_top_theory

        else:
            self.reductionobj.logger.info('order: '+str(self.reductionobj.order)+' not on detector  ')
            self.trace_success = False
            self.lhs_top = self.lhs_top_theory
        if self.trace_success:
            self.fudge_padding()
            self.smooth_spectroid()
        
        if self.cm.any(): self.trace_success = True
   
    def fudge_padding(self):
           
            #if order is very curved, add a bit more reductionobj.padding to 
            # ensure we do not interpolate into the continuum.  
            ### This should be elsewhere ###
            if abs(self.cm[0]-self.cm[-1]) > 20.:
                self.reductionobj.padding += 10. 
                
            if abs(self.cm[0]-self.cm[-1]) > 40.:
                self.reductionobj.padding += 10.
                       
    def smooth_spectroid(self):
            self.reductionobj.logger.info('fitting a polynomial function from spectroid ')
    
            # smooth out the spectroid fit         
            #fitting a 3rd order polynomial using least squares
            self.cm, p0 =  astro_math.fit_poly(self.cm, xes='default', 
                                                   deg=3)
    def shift_order(self): 
            if self.cb[0] > self.reductionobj.padding: 
                # shift the order edge trace to be in the correct place in 
                # order cut out before rectification
                # centroid top -> ct, centroid middle -> cm, centroid bottom->cb
                 self.ct = self.ct - self.cb[0] + self.reductionobj.padding
                 self.cm = self.cm - self.cb[0] + self.reductionobj.padding
                 self.cb = self.cb - self.cb[0] + self.reductionobj.padding
                 
                 self.order_shifted=True
             
            else: self.order_shifte=False

    def shift_order_back(self):
        if self.order_shifted:
            #remove the padding and start at lhs_bot to show plots in correct place
            self.cm=self.cm - self.reductionobj.padding + self.lhs_bot
            self.cb=self.cb - self.reductionobj.padding + self.lhs_bot
            self.ct=self.ct - self.reductionobj.padding + self.lhs_bot
            
    def determine_lhs_edge_pos(self, edges, theory_lhs, data_dict):
        """  find location of either top or bottom of the order using theoretical position
        as a starting point """        
        # Find the peaks in the shifted/subracted flat file near the theoretical peaks
        lhs_edge_pos = self.reductionobj.nh.get_actual_order_pos(edges, theory_lhs, 
                                          data_dict['threshold'])
 
        self.reductionobj.logger.info('found a peak at '+str(lhs_edge_pos))
        
        # Ensure the theoretical location and the actual location of 
        # the top of the order are close enough
        found_edge = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs, 
                                       data_dict['order_threshold'])

        if not found_edge:
            """ lower threshold and try again """
            self.reductionobj.logger.info('searching for top edge on a fainter flat')
            self.reductionobj.logger.info('  ---> this might affect the rectification')
            
            lhs_edge_pos = self.reductionobj.nh.get_actual_order_pos(edges, theory_lhs, 
                                              data_dict['threshold']/3)
            found_top = astro_math.actual_to_theory(lhs_edge_pos, theory_lhs, 
                                            data_dict['order_threshold'])
    
            if not found_top:
                self.reductionobj.logger.info('Flat is too faint in this order ')
                self.reductionobj.logger.info('Cannot find order location')
           
        # return variable because it could be for top or bottom location
        return lhs_edge_pos
         
 
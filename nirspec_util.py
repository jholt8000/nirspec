# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 13:23:47 2014

@author: jholt
"""
import logging, os
import numpy as np

import pylab as pl
#import pyfits as pf

class NirspecBookkeeping(object):
    '''
    a collection of methods that handle the bookkeeping details of a nirspec specific reduction
    Parameters:
    --------------
    fitsname: str
        FITS file name of science
    outpath: str
        (optional) Location of output files, defaults to current directory
    verbose: str
        (optional) print messages, defaults to True
    '''
    def __init__(self, fitsname, outpath = "", verbose = True):
        ''' Nirspec_bookkeeping needs only the FITS file name'''
        self.fitsname = fitsname
        self.outpath = outpath
        self.verbose = verbose

    @staticmethod
    def make_nirspec_final_FITS_and_plots(allreduceobj, order_num, sciorder, lineobj, traceobj, flatobj, sciobj):
        '''
        Creates final reduction outputs
        
       Parameters:
        --------------
        allreduceobj: object
            FITS file name of science
        ordernum: int
            todo
        sciorder: object
            todo
        lineobj: object
            todo
        traceobj: object
            todo
        flatobj: object
            todo
        sciobj: object
            todo

        '''
        
        pl.figure(1)
        pl.clf()
        pl.title("Order="+str(order_num))
        pl.plot(sciorder.dx, sciorder.skys,'b',label='original theory sky')
        if len(lineobj.bigohx) > 0:                    
          for i in np.arange(0,len(lineobj.matchesohx)):
            pl.text(lineobj.matchesohx[i],sciorder.skys.max()+100,str(lineobj.matchesohx[i]),rotation='vertical',fontsize='x-small')
            #pl.text(lineobj.ohx[i],sciorder.skys.max()+100,str(lineobj.ohx[i]),rotation='vertical',fontsize='x-small')
          pl.plot(lineobj.bigohx,lineobj.bigohy,'r+')
          pl.plot(lineobj.matchesohx,np.array(lineobj.matchesohy),'k+')
    
        pl.legend(loc=4)
        pl.xlabel("$\mu$")
        if allreduceobj.write_plots:
          pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'sky_order_'+str(order_num)+'.png',bbox_inches=0)
        #if write_fits:
        #  pf.writeto(outpath+sciname+'extracted_order_num'+str(order_num)+'.fits',np.array((sciorder.dx,ex2)),reductionobj.sciheader,output_verify='warn',clobber=True)

        
        pl.figure(2)
        pl.clf()
        pl.title("Order="+str(order_num))
        pl.plot(sciorder.dx_2dfit,sciorder.skys,'b',label='sky after 2dfit')
        if len(lineobj.bigohx) > 0:                    
          for i in np.arange(0,len(lineobj.matchesohx)):
            pl.text(lineobj.matchesohx[i],sciorder.skys.max()+100,str(lineobj.matchesohx[i]),rotation='vertical',fontsize='x-small')
            #pl.text(lineobj.ohx[i],sciorder.skys.max()+100,str(lineobj.ohx[i]),rotation='vertical',fontsize='x-small')
          pl.plot(lineobj.bigohx,lineobj.bigohy,'r+')
          pl.plot(lineobj.matchesohx,np.array(lineobj.matchesohy),'k+')
    
        pl.legend(loc=4)
        pl.xlabel("$\mu$")
        if allreduceobj.write_plots:
          pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'sky_order_'+str(order_num)+'.png',bbox_inches=0)
        #if write_fits:
        #  pf.writeto(outpath+sciname+'extracted_order_num'+str(order_num)+'.fits',np.array((sciorder.dx_2dfit,ex2)),reductionobj.sciheader,output_verify='warn',clobber=True)
    
        pl.figure(3)
        pl.clf()
        pl.title("Order="+str(order_num))
        pl.plot(sciorder.crosscut,'g')
        #pl.plot(fit,'k.-',label='gaussian fit')
        pl.plot(sciorder.peak,sciorder.crosscut[sciorder.peak],'k*',markersize=5,label='peak')
        pl.ylim((sciorder.crosscut[1].min()-20000,1.1*sciorder.crosscut[sciorder.peak]))
        pl.plot((sciorder.peak+sciorder.ext_range[0],sciorder.peak+sciorder.ext_range[-1]),(sciorder.crosscut[sciorder.peak],sciorder.crosscut[sciorder.peak]),'r',label='extraction window')
        if sciorder.sky_range_top: pl.plot((sciorder.peak+sciorder.sky_range_top[0],sciorder.peak+sciorder.sky_range_top[-1]),(sciorder.crosscut[sciorder.peak],sciorder.crosscut[sciorder.peak]),'b',label='sky')
        if sciorder.sky_range_bot: pl.plot((sciorder.peak+sciorder.sky_range_bot[0],sciorder.peak+sciorder.sky_range_bot[-1]),(sciorder.crosscut[sciorder.peak],sciorder.crosscut[sciorder.peak]),'b')
        pl.legend(loc=4,prop={'size':8})
        if allreduceobj.write_plots:
            pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'crosscut_order'+str(order_num)+'.png',bbox_inches=0)
        if allreduceobj.write_fits:
            #pf.writeto(allreduceobj.outpath+allreduceobj.sciname+'crosscut_order'+str(order_num)+'.fits',np.array([sciorder.crosscut]),allreduceobj.sciheader,output_verify='warn',clobber=True)                           
            c1=pf.Column(name='cross cut', format='E', array=sciorder.crosscut, bscale=4.4, ascii=True)
            hdu = pf.TableHDU.from_columns([c1])
            hdu.writefits(allreduceobj.outpath+allreduceobj.sciname+'crosscut_order'+str(order_num)+'.fits')
       
        if sciorder.cont.any():
            pl.figure(4,figsize=(15,8))
            pl.clf()
            pl.title("Order="+str(order_num))
            ax1=pl.subplot(211)
            pl.title("Order="+str(order_num))
            pl.plot(sciorder.dx_2dfit,sciorder.cont,'r',label='avg of central rows')
            pl.xlim([sciorder.dx_2dfit[0],sciorder.dx_2dfit[-1]])
            pl.subplot(212,sharex=ax1)
            
            pl.xlabel("$\mu$")
        else:
            pl.figure(4,figsize=(15,1.5))
            pl.title("Order="+str(order_num))
            pl.xlabel("$\mu$")
        pl.imshow(sciorder.rectified,origin='lower',extent=[sciorder.dx_2dfit[0],sciorder.dx_2dfit[-1],0,sciorder.rectified.shape[0]],aspect='auto')
        pl.xlabel("$\mu$")
        if allreduceobj.write_plots:
            pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'rectified_order'+str(order_num)+'.png',bbox_inches=0)
        if allreduceobj.write_fits:
            pf.writeto(allreduceobj.outpath+allreduceobj.sciname+'rectified_order'+str(order_num)+'.fits',sciorder.rectified,allreduceobj.sciheader, output_verify='warn',clobber=True)
                
        if True:
            if not traceobj.traced_bot: traceobj.cb=0
            if not traceobj.traced_top: traceobj.ct=0
        
            pl.figure(18)
            pl.clf()
            pl.imshow(flatobj.data,origin='lower')
            pl.plot(10,traceobj.lhs_bot,'c*')
            pl.plot(10,traceobj.lhs_top,'g*')
            pl.plot(10,traceobj.lhs_top_theory,'r+')
            pl.plot(10,traceobj.lhs_bot_theory,'r+')    
            pl.plot(traceobj.cm,'k',lw=2)            
            pl.plot(traceobj.ct,'g',lw=2)
            pl.plot(traceobj.cb,'g',lw=2)
            if allreduceobj.write_plots:
                pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'allflat'+str(order_num)+'.png',bbox_inches=0)       
                            
            pl.figure(19)
            pl.clf()
            pl.imshow(sciobj.data,origin='lower')
            pl.plot(10,traceobj.lhs_bot,'c*')
            pl.plot(10,traceobj.lhs_top,'g*')
            pl.plot(10,traceobj.lhs_top_theory,'r+')
            pl.plot(10,traceobj.lhs_bot_theory,'r+')    
            #pl.plot(cm,'k')            
            pl.plot(traceobj.ct,'g')
            pl.plot(traceobj.cb,'g')
            if allreduceobj.write_plots:
                pl.savefig(allreduceobj.outpath+allreduceobj.sciname+'allsci'+str(order_num)+'.png',bbox_inches=0)       
             
            if allreduceobj.show_plot:               
                pl.show() 


    def close_logger(self):
        '''
        close the logging utility
        '''
        log = logging.getLogger(__name__)
        x = list(log.handlers)
        for i in x:
            log.removeHandler(i)
            i.flush()
            i.close()
            
    def setup_logger(self):
        '''
        creates logging utility, returns logger object
        '''
        log=logging.getLogger(__name__)
        x=list(log.handlers)   
        for i in x:
            log.removeHandler(i)
            i.flush()
            i.close()
    
        if '/' in self.fitsname: fitsnameforlog=self.fitsname.split('/')[-1]
        else: fitsnameforlog=self.fitsname
        if os.path.isfile(self.outpath+'/drp_'+fitsnameforlog.rstrip('.fits')+'.log'): 
            if os.path.isfile(self.outpath+'/drp_'+fitsnameforlog.rstrip('.fits')+'_old.log'):
                os.remove(self.outpath+'/drp_'+fitsnameforlog.rstrip('.fits')+'_old.log')
            os.rename(self.outpath+'/drp_'+fitsnameforlog.rstrip('.fits')+'.log',self.outpath+'drp_'+fitsnameforlog.rstrip('.fits')+'_old.log')
               
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter=logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')    
        
        if not logger.handlers:
            fh = logging.FileHandler(filename=self.outpath+'drp_'+fitsnameforlog.rstrip('.fits')+'.log') 
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)      
          
            if self.verbose: 
                ch = logging.StreamHandler()
                ch.setLevel(logging.DEBUG)
                ch.setFormatter(formatter)
                logger.addHandler(ch)    
            else:
                logger.propagate = False
        
        return logger
    

class NirspecHeader(object):
    '''
    A collection of methods that deal with NIRSPEC header bookeeping and
    creation of tables for order locations based on header keywords.
    '''
    
    def __init__(self, header):
        self.header=header
        self.filter_name = self.header['filname']

    def get_data_dict(self):
        '''
        create a data dictionary of constants to use in NIRSPEC reduction
        '''
        # if tall orders and low filter, do not pad for interpolated shift
        # padding around order before rectification, ensures interp shift doesn't cut into data
        if '0.288x24' in self.header['slitname']:
            padding_modify = 0.0
        else:
            padding_modify = 1.0
            
        if 'NIRSPEC-1' in self.filter_name: 
            data_dict={'filt':1,'startord':80,'padding':0 * padding_modify,'threshold':1000,'spw':5.,'trace_delta':1.5}
        elif 'NIRSPEC-2' in self.filter_name: 
            data_dict={'filt':2,'startord':70,'padding':0 * padding_modify,'threshold':3000,'spw':5.,'trace_delta':1.5}
        elif 'NIRSPEC-3' in self.filter_name:
            data_dict={'filt':3,'startord':67,'padding':0 * padding_modify,'threshold':10000,'spw':3.,'trace_delta':1.1}
        elif 'NIRSPEC-4' in self.filter_name:
            data_dict={'filt':4,'startord':61,'padding':10 + 5 * padding_modify,'threshold':10000,'spw':3.,'trace_delta':1.1}
        elif 'NIRSPEC-5' in self.filter_name:
            data_dict={'filt':5,'startord':53,'padding':10 + 5 * padding_modify,'threshold':10000,'spw':3.,'trace_delta':1.1}
        elif 'NIRSPEC-6' in self.filter_name:
            data_dict={'filt':6,'startord':49,'padding':25,'spw':3.,'threshold':10000,'trace_delta':1.1}
        elif 'NIRSPEC-7' in self.filter_name:
            data_dict={'filt':7,'startord':41,'padding':30,'threshold':30000,'spw':3.,'trace_delta':1.1}   
        else:
            data_dict={'filt':self.filter_name,'startord':9999,'padding':9999,'threshold':9999,'spw':9999,'trace_delta':9999}
        
        if '24' in self.header['slitname']:
            data_dict['order_threshold'] = 50
            data_dict['order_threshold_faint_flat'] = 60
        else:
            data_dict['order_threshold'] = 40
            data_dict['order_threshold_faint_flat'] = 40
        
        print 'data_dict1=',data_dict
        return data_dict        

    def get_theory_order_pos(self, order, A_to_mu=True):
        '''
        use grating equation with coefficients empirically found for each 
        filter, grating angle, and echelle angle to determine starting 
        locations to look for each order on array
        
        '''
        if 'NIRSPEC-7' in self.filter_name:
          c1=0.24792775; c2=-35906.947; y0=15955.4515;r1=0.23482994; r2=-33591.707; z0=14891.3158
        elif 'NIRSPEC-6' in self.filter_name:
           c1=0.24986411; c2=-35961.453; y0=15944.8337;r1=0.23686484; r2=-33685.61; z0=14901.32
        elif 'NIRSPEC-5' in self.filter_name:
           c1=0.37030507; c2=-36304.523; y0=16215.28;r1=0.35065575; r2=-33928.81; z0=15115.5859
        elif 'NIRSPEC-4' in self.filter_name:
          c1=0.37318182; c2=-36031.358; y0=16009.5052;r1=0.35333187; r2=-33646.524; z0=14908.2313
        elif 'NIRSPEC-3' in self.filter_name:
           c1=0.37495832; c2=-36067.086; y0=15987.1927;r1=0.35514352; r2=-36175.709; z0=16283.995 # subtract 20 from z0,y0
        elif 'NIRSPEC-2' in self.filter_name:
           c1=0.49449345; c2=-35965.39; y0=15987.1423;r1=0.46779492; r2=-31416.355; z0=13601.3478 # subtract 20 from z0,y0
        elif 'NIRSPEC-1' in self.filter_name:
           c1=0.49777509; c2=-38653.878; y0=17488.344;r1=0.4713783; r2=-38876.842; z0=17880.5877 # subtract 20 from z0,y0
 
        x1=np.arange(1024.)
        
        if A_to_mu: multiplier = 1/10000.
        else: multiplier = 1.
        
        WL50=(852800.*np.sin(np.radians(self.header['echlpos']))- \
            24.13*(512-x1)*np.cos(np.radians(self.header['echlpos'])))/order
        
        WL50_2=(852800.*np.sin(np.radians(self.header['echlpos']))- \
            24.13*(512-50)*np.cos(np.radians(self.header['echlpos'])))/order    
        
        lhs_mid=c1*WL50_2+c2*np.sin(np.radians(self.header['disppos']))+y0
        #rhs_mid=r1*WL50_2+r2*np.sin(np.radians(self.header['disppos']))+z0
        
        lhs_top=lhs_mid+0.5*(99.488-(1.0517*self.header['disppos']))
        lhs_bot=lhs_mid-0.5*(99.488-(1.0517*self.header['disppos']))
        
        #rhs_top=rhs_mid+0.5*(99.488-(1.0517*self.header['disppos']))
        #rhs_bot=rhs_mid-0.5*(99.488-(1.0517*self.header['disppos']))
    
        # if slit is "taller", account for that
        if '24' in self.header['slitname']: 
            lhs_top=lhs_top+20.
            lhs_bot=lhs_bot-20                
    
        elif '42x' in self.header['slitname']: 
            lhs_top=lhs_top+30.
            lhs_bot=lhs_bot-30.       
     
        # The starting wavelength solution drifts with order 
        c = 0.
        
        if order > 55:#filter=3     
            WL50=WL50-137. #
            WL50=0.96465*WL50 + 412.67 # This might not work with new WL50 multiplier
            
        elif order < 38:
            c=70.
        elif order < 40:
            c=70.
        elif order < 45:
            c=60.
        elif order < 51:
            c=50.
        elif order < 55:
            c=50.
        else:
            c=20.
            
        WL50 = WL50+c
        #WL50 = multiplier * WL50
        
        return lhs_top, lhs_bot, WL50
                        
    def sanity_check(self, flatheader):
        '''
        ensure the flat and science have matching grating and echelle angles
        filter names and slitname
        requires a logger to have been set up with 
        '''
       # Sanity checking, this will not be needed if main is called from a calib manager 
        if self.header['disppos'] != flatheader['disppos'] or \
          self.header['echlpos'] != flatheader['echlpos'] or \
            self.header['filname'] != flatheader['filname'] or \
            self.header['slitname'] != flatheader['slitname']:
                    str1 = '''ERROR: the flat field does not match science. 
                    science : 
                        disppos = str(self.header['disppos'])
                        echlpos = str(self.header['echlpos'])
                        filtname = str(self.header['filtname'])
                        slitname = str(self.header['slitname'])
                    flat :
                        disppos =  str(self.flatheader['disppos'])
                        echlpos = str(self.flatheader['echlpos'])
                        filtname = str(self.flatheader['filtname'])
                        slitname = str(self.flatheader['slitname'])
                    '''
        else:
                    str1 = ''
        return str1

    @staticmethod
    def get_actual_order_pos(edges, theory, threshold):
        
        # older versions of scipy do not have this
        # I only require it when needed
        from scipy.signal import argrelextrema

        # take a vertical cut of edges 
        magcrosscut=np.sum(edges[:,0:10],axis=1)
        
        # find the highest peaks in crosscut, search +/- 15 pixels to narrow down list
        extrema=argrelextrema(magcrosscut,np.greater,order=15)[0]
        
        # find crosscut values at those extrema
        magcrosscutatextrema=magcrosscut[extrema]
        
        # narrow down extrema list to only ones over threshold
        peaks=np.where(magcrosscutatextrema > (threshold))
             
        actualpeaks=extrema[peaks[0]]
        print 'magcrosscut[extrema]=',magcrosscut[extrema]
        print 'peaks = ',peaks
        print 'actual peaks=',actualpeaks
        pl.figure(13)
        pl.clf()
        
        pl.plot(magcrosscut)
        pl.plot(extrema,magcrosscutatextrema,'bx')
        pl.plot(actualpeaks,magcrosscut[actualpeaks],'r*')
        pl.show()
        
        if actualpeaks.any(): 
            actual=min((abs(theory - i),i) for i in actualpeaks)[1]        
        else: return -3000
     
        return actual    
 
 
 


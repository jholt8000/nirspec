#!/home/koaadmin/Ureka/variants/common/bin/python

##!/usr/local/edp-python-2.7.2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:55:18 2013

@author: jholt
"""
import nirspecOO,os,subprocess
reload(nirspecOO)
import nirspec_utils
reload(nirspec_utils)
import array_manipulate

import pyfits as pf
from datetime import datetime
import optparse

def main(datadir,outdir,utdate='default'):
    '''
    wrapper script
    goes through datadir and finds any fits files, for each lev0 sorts into science and flats using KOAIMTYP 
    matches a flat to each science using echelle pos, grating position, filter, and disperser header keywords
    runs nirspec reduction script
    -> flat is only used to determine order edge locations, data are rectified, wavelength cal, and extracted
    '''
    if outdir[-1] != '/': outdir=outdir+'/'
    if datadir[-1] != '/': datadir=indir+'/'
    
    startTime=datetime.now()            
    #os.system("logger -p local3.debug nirespec_drp: Started")
    if not os.path.exists(datadir):
        print 'no dir '+datadir
        #os.system("logger -p local3.debug nirspec_drp: Directory not found"+datadir)
        return
    if utdate=='default': # first check the data directory name
      try: 
        utdate=str(datadir.split("/")[-4])
      except:
        #use today's date
        dt_ut=datetime.utcnow()
        utdate=str(dt_ut.year)+str(dt_ut.month).zfill(2)+str(dt_ut.day).zfill(2)

    year=int(utdate[0:4])
    month=int(utdate[4:6])
    day=int(utdate[6:8])
    
    if (year < 1990 or year > 2100) or (month < 1 or month >12) or (day < 1 or day > 31):
         #os.system("logger -p local3.debug nirspec_drp: Incorrect date format")
         print 'incorrect date '+utdate
         
    lev0list,err=subprocess.Popen(["find "+datadir+" -name \*fits\* | sort"],stdout=subprocess.PIPE,shell=True).communicate()
    #lev0list=os.system("find "+datadir+" -name \*fits\* | sort")
    lev0list=lev0list.split('\n')
    if not lev0list: 
        print 'no lev0 fits'
        #os.system("logger -p local3.debug nirspec_drp: No lev0 fits files found in "+datadir)
        return
        
    #make a list of koaimtype objects
    #make a list of koaimtype flatlamps
    science_files={}
    flat_files={}
    dark_files={}
    #print 'lev0list=',lev0list
    for fitsfile in lev0list:
        if len(fitsfile) < 3: 
            continue
    
        if fitsfile.endswith('gz'): os.system('gunzip '+fitsfile)
        fitsfile=fitsfile.rstrip('.gz')
        header=pf.getheader(fitsfile)
        try:
            koaimtyp = header['IMAGETYP']
            disppos = header['disppos']  
            thetaE = header['echlpos']
            filname = header['filname']
            slitname = header['slitname']  
            dispers=header['dispers']
            elaptime = header['ELAPTIME']
        except:
            print 'not lev0 NIRSPEC data'
            return
        
        if dispers == 'low':
            print fitsfile+'low dispersion'
            continue 
        #science_files[fitsfile]=[]
        # make lists of all object (science) files and flat files in the directory
        if koaimtyp=='object':
            science_files[fitsfile]=[]
            science_files[fitsfile].append((disppos,thetaE,filname,slitname,elaptime))
        elif koaimtyp=='flatlamp':
            flat_files[fitsfile]=[]
            flat_files[fitsfile].append((disppos,thetaE,filname,slitname))
        elif koaimtyp == 'dark':
            dark_files[fitsfile]=[]
            dark_files[fitsfile].append((elaptime))            
        else:
            print fitsfile+' koaimtype=',koaimtyp
    
    print 'science_files=',science_files
    # associate a matching dark and flat with each science frame
    for name,val in science_files.iteritems():
        print science_files[name]
        
        science_files[name].append([])
        
        for dname, dval in dark_files.iteritems():
            if val:                
                try: 
                    if val[0][4]==dval[0]: # elaptime
                        science_files[name][1].append(dname)
                except:
                    pass
                
        #placeholder in tuple for file name if no matching darks found

        if len(science_files[name]) < 2: 
            science_files[name][1].append('nodarks')

        science_files[name].append([])

        for fname,fval in flat_files.iteritems():
            if val: 
                try: 
                    if val[0][:4]==fval[0]: 
                        science_files[name][2].append(fname)
                except:
                    pass

        if len(science_files[name]) < 3: 
            science_files[name][2].append('noflats')
        
        print 'flats found=',science_files[name][2]
        print 'darks found=',science_files[name][1]
        

    return science_files, flat_files, dark_files
        
        
    for sfile in science_files.keys():
        print 'sfile'     
        if len(science_files[name][2]) > 10: flats = science_files[name][2][:10]
        else: flats = science_files[name][2]
        if len(science_files[name][1]) > 10: darks = science_files[name][1][:10]
        else: darks = science_files[name][1]
        
        master_flat = array_manipulate.Astro_math.median_comb(flats)
        master_dark = array_manipulate.Astro_math.median_comb(darks)   
        
        #master_flat=nirspec_utils.median_comb(flats)
        ##master_dark=nirspec_utils.median_comb(darks)
        
        
        if len(science_files[sfile]) > 1: #found matching dark
          print 'starting reduction for science file ',sfile
          filestart=datetime.now()
          print 'call=','nirspec.main('+sfile+', ',science_files[sfile][1]+')'
          if '1' in filname: 
              sky_distance=4
              sky_height=10
          else: 
              sky_distance=5
              sky_height=5
              
          c= nirspecOO.main(sfile, master_flat, dark=master_dark,
                                 do_extract=True, verbose=True, show_plot=False, 
                                 cosmic_clean=True, write_fits=False, write_plot=True,
                                 outpath=outdir,sky_distance=sky_distance, 
                                 sky_height=sky_height)
          c.reduce_nirspec()
          
          print datetime.now()-filestart
    print datetime.now()-startTime
    
if __name__=="__main__":
    usage="""
    %prog drp_nirspec.py input_dir output_dir utdate(optional)
    """
    p=optparse.OptionParser(usage=usage)
    (options,args)=p.parse_args()
    if len(args) >= 2:
        indir=args[0]
        outdir=args[1]
        if len(args) == 3: utdate=args[2]
        else: utdate='default'
    else:
        p.error('Wrong arguments')
        
 
    main(indir,outdir,utdate)
#!/home/koaadmin/Ureka/variants/common/bin/python

# #!/usr/local/edp-python-2.7.2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:55:18 2013

@author: jholt
"""
import os
import subprocess

import nirspecOO
reload(nirspecOO)
print (nirspecOO)
import astro_math
reload(astro_math)

from astropy.io import fits
from datetime import datetime
import optparse

def main(datadir, outdir, utdate='default'):
    """
    wrapper script
    goes through datadir and finds any fits files, for each lev0 sorts into science and flats using KOAIMTYP 
    matches a flat to each science using echelle pos, grating position, filter, and disperser header keywords
    runs nirspec reduction script
    -> flat is only used to determine order edge locations, data are rectified, wavelength cal, and extracted
    """
    if outdir[-1] != '/':
        outdir += '/'
    if datadir[-1] != '/':
        datadir += '/'

    if not os.path.exists(outdir):
        try: os.mkdir(outdir)
        except: print 'cannot make outdir ',outdir
        return

    startTime = datetime.now()
    if not os.path.exists(datadir):
        print 'no dir ' + datadir
        # os.system("logger -p local3.debug nirspec_drp: Directory not found"+datadir)
        return
    if utdate == 'default':  # first check the data directory name
        try:
            utdate = str(datadir.split("/")[-4])
        except:
            # use today's date
            dt_ut = datetime.utcnow()
            utdate = str(dt_ut.year) + str(dt_ut.month).zfill(2) + str(dt_ut.day).zfill(2)

    year = int(utdate[0:4])
    month = int(utdate[4:6])
    day = int(utdate[6:8])

    if (year < 1990 or year > 2100) or (month < 1 or month > 12) or (day < 1 or day > 31):
        # os.system("logger -p local3.debug nirspec_drp: Incorrect date format")
        print 'incorrect date ' + utdate

    lev0list, err = subprocess.Popen(["find " + datadir + " -name \*fits\* | sort"], stdout=subprocess.PIPE,
                                     shell=True).communicate()

    lev0list = lev0list.split('\n')
    if not lev0list:
        print 'no lev0 fits'
        # os.system("logger -p local3.debug nirspec_drp: No lev0 fits files found in "+datadir)
        return

    # make a list of koaimtyp objects
    # make a list of koaimtyp flatlamps
    science_files = {}
    flat_files = {}
    dark_files = {}
    # print 'lev0list=',lev0list
    for fitsfile in lev0list:
        if len(fitsfile) < 3:
            continue

        if fitsfile.endswith('gz'):
            os.system('gunzip ' + fitsfile)

        fitsfile = fitsfile.rstrip('.gz')
        header = fits.getheader(fitsfile)
        try:
            koaimtyp = header['IMAGETYP']
            disppos = header['disppos']
            thetaE = header['echlpos']
            filname = header['filname']
            slitname = header['slitname']
            dispers = header['dispers']
            elaptime = header['ELAPTIME']
        except:
            print 'not lev0 NIRSPEC data'
            return

        if dispers == 'low':
            print fitsfile + ' low dispersion'
            continue
            # science_files[fitsfile]=[]
        # make lists of all object (science) files and flat files in the directory
        if koaimtyp == 'object':
            science_files[fitsfile] = []
            science_files[fitsfile].append((disppos, thetaE, filname, slitname, elaptime))
        elif koaimtyp == 'flatlamp':
            flat_files[fitsfile] = []
            flat_files[fitsfile].append((disppos, thetaE, filname, slitname))
        elif koaimtyp == 'dark':
            dark_files[fitsfile] = []
            dark_files[fitsfile].append(elaptime)
        else:
            print fitsfile + ' koaimtype=', koaimtyp

    # associate a matching dark and flat with each science frame
    for name, val in science_files.iteritems():
        science_files[name].append([])

        for dname, dval in dark_files.iteritems():
            if val:
                try:
                    if val[0][4] == dval[0]:  # elaptime
                        science_files[name][1].append(dname)
                except:
                    pass

        # placeholder in tuple for file name if no matching darks found

        if len(science_files[name]) < 2:
            science_files[name][1].append('nodarks')

        science_files[name].append([])

        for fname, fval in flat_files.iteritems():
            if val:
                try:
                    if val[0][:4] == fval[0]:
                        science_files[name][2].append(fname)
                except:
                    pass

        if len(science_files[name]) < 3:
            science_files[name][2].append('noflats')

    # return science_files, flat_files, dark_files

    for sfile in science_files.keys():

        if len(science_files[sfile][2]) > 10:
            flats = science_files[sfile][2][:10]
        else:
            flats = science_files[sfile][2]
        if len(science_files[sfile][1]) > 10:
            darks = science_files[sfile][1][:10]
        else:
            darks = science_files[sfile][1]

        if len(flats) > 1:
            master_flat = astro_math.median_comb(flats)
        elif len(flats) == 1:
            master_flat = fits.getdata(flats[0])
        else:
            print "cannot reduce"+sfile+" with no matching flat!"
            continue
        if len(darks) > 1:
            master_dark = astro_math.median_comb(darks)
        elif len(darks) == 1:
            master_dark = fits.getdata(darks[0])
        else:
            master_dark = []

        if len(science_files[sfile]) > 1:
            # found matching
            print 'starting reduction for science file ', sfile, ' using flats ',flats
            filestart = datetime.now()
            if '1' in filname:
                sky_distance = 4
                sky_height = 10
            else:
                sky_distance = 5
                sky_height = 5
            #try:
            if True:
                c = nirspecOO.Main(sfile, flat_array=master_flat, dark_array=master_dark, do_extract=True, verbose=True,
                               show_plot=False, cosmic_clean=True, write_fits=False, write_plots=True, outpath=outdir,
                               sky_distance=sky_distance, sky_height=sky_height,  max_iter=3, sig_clip=5.0, sig_frac=0.3,
                               obj_lim=5.0, ext_height=3, sky_sigma=2.25, traceWidth=10,
                               backgroundWidth=30,  traceMean=True, traceLast=False, traceDelta=1.9)

                c.reduce_nirspec()
            #except:
            #    print 'problem reducing file :',sfile

            print datetime.now() - filestart
    print datetime.now() - startTime


if __name__ == "__main__":
    usage = """
    %prog input_dir output_dir utdate(optional)
    """
    p = optparse.OptionParser(usage=usage)
    (options, args) = p.parse_args()
    if len(args) >= 2:
        datadir = args[0]
        outdir = args[1]
        if len(args) == 3:
            utdate = args[2]
        else:
            utdate = 'default'
    else:
        p.error('Wrong arguments')

    main(datadir, outdir, utdate)
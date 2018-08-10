#!/usr/bin/env python

import sys,os
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table, Column


def wcs_corr(input_fits,input_prm,path=os.getcwd(),output_file=''):

    '''
    input_fits: str
        The fully reduced datacube Pampelmuse has been run on

    input_prm: str
        The prm file produced by Pampelmuse

    path: str, optional
        I/O path, Default: current directory

    output_file: str, optional
        outputfile name, Default: inputname_corr.fits

    '''

    cube=fits.open(path+'/'+input_fits)
    pmr=fits.open(path+'/'+input_prm)

    A=np.nanmedian(pmr[4].data[0][1])
    B=np.nanmedian(pmr[4].data[1][1])
    C=np.nanmedian(pmr[4].data[2][1])
    D=np.nanmedian(pmr[4].data[3][1])
    x0=np.nanmedian(pmr[4].data[4][1])
    y0=np.nanmedian(pmr[4].data[5][1])
    
    print(A,B,C,D,x0,y0)
    ### as described in the Pampelmuse manual:

    dx=  A + C + x0 + 0.5
    dy=  D - B + y0 + 0.5

    print(dx,dy)

    prihdr = cube[1].header

    prihdr['CRPIX1']=cube[1].header['CRPIX1']+dx
    prihdr['CRPIX2']=cube[1].header['CRPIX2']+dy

    cube.writeto(path+'/'+input_fits[:-5]+'_corr.fits',overwrite=True)




def pampelmuse_cat(ra, dec, mag, filter, idx=None, path=os.getcwd(),sat = 0.,mag_sat = 10.):
    
    '''
    
    May be used to create the inout catalog for pampelmuse
    
    ra: float
        RA coordinates of the stars
    
    dec: float
        Dec coordinates of the stars
    
    mag: float
        magnitudes of the stars
    
    filter: str
        filter used for the magnitudes
    
    idx: float, optional
        index of the stars, default counting up from 1
    
    sat: float, optional
        value assigned to saturated sources in the catalog, default: 0.
    
    mag_sat: float, optional
        magnitude that replaces saturated sources in the catalog, default: 10.
    
    path: str, optional
        path of the output file, default: current directory
    
    '''
    
    if idx == None:
        id=np.arange(len(ra))+1
    else
        id = idx
    
    sat_source = np.where(mag == sat)
        
    tab=Table([id,ra,dec,mag,names=('id','ra','dec',str(filter).lower()),dtype=('i4', 'f8', 'f8','f8'))    
    tab[str(filter).lower()][sat] = new_mag_sat
    
    tab.write(path+'/'+str(filter).upper(),format='ascii.basic',delimiter=',',overwrite=True)
    
    
    
    
    
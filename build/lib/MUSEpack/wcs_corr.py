#!/usr/bin/env python

import sys,os
import numpy as np
from astropy.io import fits


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

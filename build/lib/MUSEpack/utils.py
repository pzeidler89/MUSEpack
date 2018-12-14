#!/usr/bin/env python
'''
utils.py

Copyright 2018-2018 Peter Zeidler

This file contains multiple utility functions for MUSEpack.
This module is part of MUSEpack and intended to do the RV determination of MUSE spectra

MUSEpack is a free software package: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

MUSEpack is distributed in the hope that it will be useful for working with
MUSE data and spectra, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with MUSEpack. If not, see <http://www.gnu.org/licenses/>.

'''

import sys,os,shutil,warnings
import pyspeckit
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table, Column
from astropy import units as u
from astropy import log
from astropy.stats import median_absolute_deviation as MAD
import pandas as pd
import time
import logging
from pysynphot import observation,spectrum
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import find_peaks


def initial_guesses(lines, linestrength=100., sigma_gauss=2.4, sigma_lorentz=2.4, absorption=True, llimits = [2.,2.]):

    '''
    lines: ARRAY, central wavelengths of the spectral lines
    linestrength: FLOAT, initial guess for the line strength
                         default: 100
    sigma_gauss: FLOAT, initial guess for the 1 sigma of the gaussian part
                        of the Voigt profile
                        default: 2.4
    sigma_lorentz: FLOAT, initial guess for the 1 sigma of the lorentz part
                          of the Voigt profile
                          default: 2.4
    absorption: BOOL, defines if line is in absorption or in emission
                    default: True (absorption)
    return: guesses,limits,limited
    '''

    if absorption == True: linestrength = linestrength * (-1)

    guesses=np.zeros(4*len(lines))
    limits=[]
    limited=[]

    for idx,line in enumerate(lines):
    
        guesses[4*idx+0]=linestrength
        guesses[4*idx+1]=line
        guesses[4*idx+2]=sigma_gauss
        guesses[4*idx+3]=sigma_lorentz
    
        limits.append((0,0)) #Amplitude
        limited.append((False,True)) #Amplitude

        limits.append((line-llimits[0],line+llimits[1]))
        limited.append((True,True)) #Center

        limits.append((0,0)) #width (gausssigma)
        limited.append((True,False)) #width

        limits.append((0,0)) #width (lorentzsigma)
        limited.append((True,False)) #width
    
    return guesses,limits,limited


def baseline_extract(specfit,poly_order,dispersion=1.25):
    coeff=np.zeros(poly_order+1)
    power=np.arange(poly_order+1)
    for order in range(poly_order+1): coeff[order]=specfit.header['BLCOEF'+str('{:02d}'.format(order))]
    baseline=np.poly1d(coeff/dispersion**power[::-1])
    return baseline

def voigt_FWHM(sigma_g, gamma_l):
    #calculating the FWHM of the voigt, gauss, and lorentz profile
    
    c_0 = 2.0056
    c_1 = 1.0593
    fwhm_g = 2.*np.sqrt(2.*np.log(2.))*sigma_g
    fwhm_l = 2.*gamma_l
    phi = fwhm_l/fwhm_g
    fwhm_v = fwhm_g*(1.-c_0*c_1+np.sqrt(phi**2+2.*c_1*phi+(c_0*c_1)**2))
    for i,fwhm_g in enumerate(fwhm_g):
        if fwhm_g == 0: fwhm_v[i] = fwhm_l[i]

    return fwhm_g,fwhm_l,fwhm_v


def spec_res_downgrade(l_in,spec_in,l_out):
    dispersion=2.4
    templ_spec = spectrum.ArraySourceSpectrum(wave=l_in, flux=spec_in)
    white_filter = spectrum.ArraySpectralElement(l_out, np.ones(len(l_out)), waveunits='angstrom')
    convolved_spec = observation.Observation(templ_spec, white_filter, binset=l_out, force='taper').binflux
    return convolved_spec

def line_clipping(self,x,line_significants,sigma = 3):
    
    mask = np.zeros_like(x,dtype = int)
    ind_sig = np.where(self.cat.loc[:,'significance'].values.astype(np.float64) < line_significants)
    mask[ind_sig] = 1
    
    if len(x) >= 3:
        
        if len(x) == 3:
            median = np.median(x)
            mad = MAD(x)
        
        if len(x) > 3:
            x_red = np.delete(x,[np.argmin(x),np.argmax(x)])
            median = np.median(x_red)
            mad = MAD(x_red)
    
        ind = np.where((x < median - sigma*mad) | (x > median + sigma*mad))[0]
        mask[ind] = 1
    
    x_masked = np.ma.masked_array(x, mask=[mask])
    
    return x_masked
    
def clipped_MAD(x):
    if len(x) <= 3: mad = MAD(x)
    if len(x) > 3:
        x = np.delete(x,[np.argmin(x),np.argmax(x)])
        mad = MAD(x)
    
    return mad
    
    
def continuum_deviation(l_in,f_in,baseline,contorder):
    ## Calculates the relative defiation between the baseline and a rudimentary continuum
    
    peaks = find_peaks(f_in*(-1),width=2)[0]
    mask = np.zeros_like(f_in)
    
    for p in peaks: mask[p-3:p+4] = 1
    
    remove = ~np.ma.array(l_in,mask = mask).mask
    
    l_masked = l_in[remove]
    f_masked = f_in[remove]
    baseline_masked = baseline[remove]
    
    masked_fit = np.polyfit(l_masked,f_masked,contorder)
    masked_fitfunc = np.poly1d(masked_fit)
    
    cont_dev = MAD(f_masked/np.median(f_masked)-baseline_masked/np.median(baseline_masked))
    
    return cont_dev

    
    
    
    

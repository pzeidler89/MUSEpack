#!/usr/bin/env python
'''
radial_velocities.py

Copyright 2018-2018 Peter Zeidler

This file contains ultiple utility functions for main class RV_spectrum class.
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

import sys,os,shutil,warnings,pyspeckit
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table, Column
from astropy import units as u
from astropy import log
from astropy.stats import median_absolute_deviation as MAD
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from multiprocessing import Pool,cpu_count
from functools import partial
from scipy.special import wofz
from ppxf import ppxf
from ppxf import ppxf_util
import dask
from multiprocessing import cpu_count
import time
import logging

from pysynphot import observation,spectrum
from scipy.ndimage.filters import gaussian_filter


def initial_guesses(lines,linestrength=100.,sigma_gauss=2.4,sigma_lorentz=2.4,absorption=True):

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

        limits.append((line-2.,line+2.))
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
    line_width=1.25
    templ_spec = spectrum.ArraySourceSpectrum(wave=l_in, flux=spec_in)
    white_filter = spectrum.ArraySpectralElement(l_out, np.ones(len(l_out)), waveunits='angstrom')
    obs = observation.Observation(templ_spec, white_filter, binset=l_out, force='taper')

    obs_dispersion=line_width/(2.*np.sqrt(2.*np.log(2)))
    convolved_spec=gaussian_filter(obs.binflux,obs_dispersion)
    
    return convolved_spec
    
    
def line_fitting_dask(self,linecat,line_idx,niter,resid_level,max_contorder,adjust_preference):

    self.logger.info('Started line '+line_idx)
    
    std_resid = resid_level + 1.
    lstart = self.cat.loc[line_idx,'l_start']
    lend = self.cat.loc[line_idx,'l_end']
    contorder = self.cat.loc[line_idx,'cont_order']
    
    while resid_level < std_resid:
        
        
        lines_select = linecat[np.where((linecat >= lstart) & (linecat < lend))]
        spec_select_idx = np.where((self.spec_lambda >= lstart) &\
                                   (self.spec_lambda < lend))

        linefit_guess, linefit_limits, linefit_limited = initial_guesses(lines_select)
        exponent=int(np.log10(np.nanmedian(self.spec_f[spec_select_idx]))-4.) # Obtain larger flux numbers to make it numerically stable

        factor=10**(-exponent)

        template_f  = np.zeros_like(self.spec_lambda,dtype=float)
        fit_f       = np.zeros_like(self.spec_lambda,dtype=float)
        temp_f      = np.array(self.spec_f[spec_select_idx]*factor,dtype=float)
        temp_err    = np.array(self.spec_err[spec_select_idx]*factor,dtype=float)
        temp_lambda = np.array(self.spec_lambda[spec_select_idx],dtype=float)

        ### Initiate the spectrum for pyspeckit

        
        sp = pyspeckit.Spectrum(data=temp_f,error=temp_err,xarr=temp_lambda,\
                                unit=r'flux [$10^{'+str(exponent)+'}$ erg$^{-1}$ s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',header={})
        sp.xarr.set_unit(u.AA)
        sp.xarr.xtype='wavelength'

        ### Initiate the spectral fit

        exclusion_level=0.01
        iterations=0
        
        while iterations < niter:

            with log.log_to_list() as log_list:
                sp.baseline(order=contorder,plot_baseline=False,subtract=False, annotate=False,highlight_fitregion=False, excludefit=True,\
                            save=True,exclude=None,fit_plotted_area=False,exclusionlevel=exclusion_level)
                sp.specfit.multifit(fittype='voigt', renormalize='auto', annotate=False, show_components=False,verbose=False, color=None, guesses=list(linefit_guess),\
                                    parinfo=None, reset_fitspec = True,plot=False, limits=linefit_limits, limited=linefit_limited)
            
                if len(log_list) > 0 and log_list[0].message[:8] == 'gnorm=0.':
                    exclusion_level += 0.01
                    self.logger.info(line_idx+': Adjusting exclusion level to: '+str(exclusion_level))
                    iterations = 0
                    sp = pyspeckit.Spectrum(data=temp_f,error=temp_err,xarr=temp_lambda, unit=r'flux [$10^{'+str(exponent)+'}$ erg$^{-1}$ s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]',header={})
                    sp.xarr.set_unit(u.AA)
                    sp.xarr.xtype='wavelength'

                if len(log_list) == 0 or log_list[0].message[:8] != 'gnorm=0.':
                    if iterations > 0 and sp.specfit.optimal_chi2(reduced=True, threshold='auto') == chi2: break ##leavin the loop if it prematurely convereged
                    chi2=sp.specfit.optimal_chi2(reduced=True, threshold='auto')
                    self.logger.debug(line_idx+':Iteration: '+str(iterations)+' Chi2: '+str(chi2))
                    iterations += 1
        
        
        self.logger.info(line_idx+': Converged after '+str(iterations)+' iterations; Chi2: '+str(chi2))
        

        par_extract_idx = []
        for lab_line_idx,lab_lines in enumerate(np.array([self.cat.loc[line_idx,'l_lab']])):
            par_extract_idx=np.concatenate((par_extract_idx,np.where(lines_select == lab_lines)[0]))

        temp_sg = np.zeros_like(np.array([self.cat.loc[line_idx,'l_lab']]))
        temp_sl = np.zeros_like(np.array([self.cat.loc[line_idx,'l_lab']]))
        temp_a  = np.zeros_like(np.array([self.cat.loc[line_idx,'l_lab']]))
        temp_l  = np.zeros_like(np.array([self.cat.loc[line_idx,'l_lab']]))

        temp_idx = 0
        highres_lambda= np.linspace(self.spec_lambda[0], self.spec_lambda[-1], 100000)
        
        for i in range(len(lines_select)):
            xcen = lines_select[i]
            gamma = sp.specfit.parinfo[int(4*i+3)]['value']
            sigma = sp.specfit.parinfo[int(4*i+2)]['value']
            amp = sp.specfit.parinfo[int(4*i+0)]['value']
                        
            # z = ((self.spec_lambda-xcen) + 1j*gamma) / (sigma * np.sqrt(2))
            # lineflux = amp * np.real(wofz(z))
            # re = np.real(wofz(z))
            
            z = ((highres_lambda-xcen) + 1j*gamma) / (sigma * np.sqrt(2))
            highres_f = np.real(wofz(z))
            lowres_f = spec_res_downgrade(highres_lambda,highres_f ,self.spec_lambda)
            lineflux = amp * lowres_f
            template_f += lineflux/factor
            
            z = ((highres_lambda-sp.specfit.parinfo[int(4*i+1)]['value']) + 1j*gamma) / (sigma * np.sqrt(2))
            highres_f = np.real(wofz(z))
            lowres_f = spec_res_downgrade(highres_lambda,highres_f,self.spec_lambda)
            lineflux1 = amp * lowres_f
            fit_f += lineflux1/factor
                
                
            if (i == par_extract_idx).any():

                temp_a[temp_idx] = sp.specfit.parinfo[int(4*i+0)]['value']
                temp_l[temp_idx] = sp.specfit.parinfo[int(4*i+1)]['value']
                temp_sg[temp_idx] = sp.specfit.parinfo[int(4*i+2)]['value']
                temp_sl[temp_idx] = sp.specfit.parinfo[int(4*i+3)]['value']
                temp_idx += 1


        baseline_temp=baseline_extract(sp,contorder)
        baseline_sub_spec=sp.data/baseline_temp(temp_lambda-temp_lambda[0])
        continuum = baseline_temp(temp_lambda-temp_lambda[0])/factor
        resid = (self.spec_f[spec_select_idx] - fit_f[spec_select_idx])/continuum
        
            
        std_resid = np.std(resid)
        std_resid =0.01
        if std_resid > resid_level:
            max_contorder_reached = False

            if  adjust_preference == 'contorder':
                if contorder < max_contorder:
                    contorder += 1
                    self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' => Adjusting continuum order to: '+str(contorder))

                else:
                    max_contorder_reached = True
                    self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' => maximum cont order reached: Adjusting wavelength range')

            if  adjust_preference == 'min_lambda' or adjust_preference == 'minmax_lambda' or max_contorder_reached:
                lstart += 5.
                self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' => Adjusting lower wavelength range to: '+str(lstart))

            if  adjust_preference == 'max_lambda' or adjust_preference == 'minmax_lambda' or max_contorder_reached:
                lend += 5.
                self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' => Adjusting upper wavelength range to: '+str(lend))
    
    self.logger.info('Finished line '+line_idx)
    
    return line_idx,temp_l,temp_a,temp_sl,temp_sg,spec_select_idx,template_f,continuum,lstart,lend,contorder,fit_f


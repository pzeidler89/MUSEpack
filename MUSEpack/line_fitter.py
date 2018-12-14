'''
line_fitter.py

Copyright 2018-2018 Peter Zeidler

This module is part of MUSEpack and intended to do the line fitting of MUSE spectra

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
from astropy import units as u
from astropy import log
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from scipy.special import wofz
import time
import logging

from MUSEpack.utils import *


def line_fitter(self,linecat,line_idx,niter,input_resid_level,max_contorder,max_ladjust,adjust_preference,input_continuum_deviation,llimits):

    self.logger.info('Started line '+line_idx)
    
    resid_level = 5.
    std_resid = resid_level + 1.
    
    continuum_dev = input_continuum_deviation + 1.
    
    lstart = self.cat.loc[line_idx,'l_start']
    lend = self.cat.loc[line_idx,'l_end']
    contorder = self.cat.loc[line_idx,'cont_order']

    n_ladjust = 0
    max_cont_reached = False
    max_n_ladjust_reached = False
    significant = True

    
    while std_resid > resid_level or continuum_dev > input_continuum_deviation:
        
        lines_select = linecat[np.where((linecat >= lstart) & (linecat < lend))]
        spec_select_idx = np.where((self.spec_lambda >= lstart) &\
                                   (self.spec_lambda < lend))

        linefit_guess, linefit_limits, linefit_limited = initial_guesses(lines_select,llimits = llimits)
        exponent=int(np.log10(np.nanmedian(self.spec_f[spec_select_idx]))-4.) # Obtain larger flux numbers to make it numerically stable

        factor=10**(-exponent)

        template_f  = np.zeros_like(self.spec_lambda,dtype=float)
        fit_f       = np.zeros_like(self.spec_lambda,dtype=float)
        temp_f      = np.array(self.spec_f[spec_select_idx]*factor,dtype=float)
        temp_err    = np.array(self.spec_err[spec_select_idx]*factor,dtype=float)
        temp_lambda = np.array(self.spec_lambda[spec_select_idx],dtype=float)

        if input_resid_level == None:
            resid_level = np.std(temp_err/np.max(temp_err))
        else:
            resid_level = input_resid_level
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
                    # if iterations > 0 and (sp.specfit.optimal_chi2(reduced=True, threshold='auto') == chi2 or np.abs(re_chi2_change) < 0.05): break ##leavin the loop if it prematurely convereged
                    if iterations > 0:
                        chi2_change = (chi2-sp.specfit.optimal_chi2(reduced=True, threshold='auto'))/sp.specfit.optimal_chi2(reduced=True, threshold='auto')
                        if sp.specfit.optimal_chi2(reduced=True, threshold='auto') == chi2 or np.abs(chi2_change) < 0.05: break ##leavin the loop if it prematurely convereged
                        self.logger.debug(line_idx+':Iteration: '+str(iterations)+' delta Chi2: '+str('{:2.3f}'.format(chi2_change)))
                        
                    chi2=sp.specfit.optimal_chi2(reduced=True, threshold='auto')
                    self.logger.debug(line_idx+':Iteration: '+str(iterations)+' Chi2: '+str('{:3.3f}'.format(chi2)))
                    iterations += 1
        
        if iterations < niter: self.logger.info(line_idx+': Converged after '+str(iterations)+' iterations; Chi2: '+str('{:3.3f}'.format(chi2))+' delta Chi2 [%]: '+str('{:2.3f}'.format(chi2_change)))
        if iterations == niter: self.logger.info(line_idx+': Maximum number of iterations ('+str(iterations)+') reached; Chi2: '+str(chi2))
        

        par_extract_idx = []
        for lab_line_idx,lab_lines in enumerate(np.array([self.cat.loc[line_idx,'l_lab']])):
            par_extract_idx=np.concatenate((par_extract_idx,np.where(lines_select == lab_lines)[0]))

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

                temp_a = sp.specfit.parinfo[int(4*i+0)]['value']/factor
                temp_l = sp.specfit.parinfo[int(4*i+1)]['value']
                temp_sg = sp.specfit.parinfo[int(4*i+2)]['value']
                temp_sl = sp.specfit.parinfo[int(4*i+3)]['value']

        baseline_temp=baseline_extract(sp,contorder)
        baseline_sub_spec=sp.data/baseline_temp(temp_lambda-temp_lambda[0])
        continuum = baseline_temp(temp_lambda-temp_lambda[0])/factor
        resid = (self.spec_f[spec_select_idx] - fit_f[spec_select_idx])/continuum
            
        std_resid = np.std(resid)
        
        continuum_dev = continuum_deviation(temp_lambda,temp_f,continuum,contorder)

        if std_resid > resid_level or continuum_dev > input_continuum_deviation:

            if std_resid > resid_level: self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' targeted: '+str('{:2.4f}'.format(resid_level)))
            if continuum_dev > input_continuum_deviation: self.logger.info(line_idx+': Continuum deviation: '+str('{:2.4f}'.format(continuum_dev))+' targeted: '+str('{:2.4f}'.format(input_continuum_deviation)))
            
            if contorder == max_contorder: max_cont_reached = True
            if n_ladjust == max_ladjust: max_n_ladjust_reached = True
            
            if max_cont_reached and max_n_ladjust_reached:
                self.logger.warning(line_idx+": Maximum adjustments reached: Check the output")
                break
                
            if adjust_preference == 'contorder':
                
                if not max_cont_reached:
                    contorder += 1
                    self.logger.info(line_idx+': Adjusting continuum order to: '+str(contorder))
                    
                if max_cont_reached and not max_n_ladjust_reached:
                    lstart -= 5.
                    lend += 5.
                    self.logger.info(line_idx+': maximum continuum order reached => Adjusting wavelength range to: ['+str(lstart)+','+str(lend)+']')
                    n_ladjust += 1
                    contorder = self.cat.loc[line_idx,'cont_order']
                     
            if adjust_preference != 'contorder':
            
                contorder = self.cat.loc[line_idx,'cont_order']
            
                if not max_n_ladjust_reached:
                    
                    if adjust_preference == 'min_lambda' or adjust_preference == 'minmax_lambda':
                        lstart -= 5.
                        self.logger.info(line_idx+': Adjusting lower wavelength range to: '+str(lstart))

                    if adjust_preference == 'max_lambda' or adjust_preference == 'minmax_lambda':
                        lend += 5.
                        self.logger.info(line_idx+': Adjusting upper wavelength range to: '+str(lend))
                    
                    n_ladjust += 1
                if max_n_ladjust_reached and not max_cont_reached:
                    
                    contorder += 1
                    self.logger.info(line_idx+': maximum wavelength adjustment reached => Adjusting continuum order to: '+str(contorder))
                                     
                    lstart = self.cat.loc[line_idx,'l_start']
                    lend = self.cat.loc[line_idx,'l_end']

    significance = np.abs(temp_a)/np.median(temp_err/factor)
    
    self.logger.info(line_idx+': STD of residuals: '+str('{:2.4f}'.format(std_resid))+' targeted: '+str('{:2.4f}'.format(resid_level)))
    self.logger.info(line_idx+': Continuum deviation: '+str('{:2.4f}'.format(continuum_dev))+' targeted: '+str('{:2.4f}'.format(input_continuum_deviation)))
    self.logger.info(line_idx+' Line significance: '+str('{:3.2f}'.format(significance)))
    self.logger.info('Finished line '+line_idx)
    
    return line_idx,temp_l,temp_a,temp_sl,temp_sg,spec_select_idx,template_f,continuum,lstart,lend,contorder,fit_f,significance
###            0      1        2    3       4           5           6           7         8     9       10      11      12
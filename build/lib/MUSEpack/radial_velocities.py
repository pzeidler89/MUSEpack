#!/usr/bin/env python

"""
radial_velocities.py

Copyright 2018-2018 Peter Zeidler

This file contains the main class for the RV fitting.
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

Required non-standard packages: ppxf, pyspeckit, pysynphot

"""

__version__ = '0.1.0'

__revision__ = '20191226'



import sys,os,shutil,warnings,pyspeckit,pickle
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column
from astropy import units as u
from astropy import constants as const
from astropy import log
from astropy.stats import median_absolute_deviation as MAD
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from multiprocessing import cpu_count
from functools import partial
from scipy.special import wofz
from ppxf import ppxf
from ppxf import ppxf_util
import dask
import time
import logging

from pysynphot import observation,spectrum
from scipy.ndimage.filters import gaussian_filter

''' internal modules'''
from MUSEpack.ppxf_MC import ppxf_MC
from MUSEpack.utils import *
from MUSEpack.line_fitter import *



logging.basicConfig(format = '%(asctime)s [%(levelname)8s]: %(message)s [%(funcName)s]',datefmt="%Y-%m-%d %H:%M:%S")

class RV_spectrum:
    """
    Args:
        spec_id (str): unique identifier of the spectrum
        spec_f (array): flux array of the spectrum in erg/s/cm^2/Angstrom
        spec_err (array): flux uncertainty array of the spectrum in erg/s/cm^2/Angstrom
        spec_lambda (array): wavelength_array of the spectrum in Angstrom
        loglevel (str): outputlevel of the logger: INFO, DEBUG optional, default: INFO
            when DEBUG is activated, all functions run on a single core to obtain
            proper output in the correct order.
    
    """
    
    def __init__(self,spec_id,spec_f,spec_err,spec_lambda,logger=None,loglevel = "INFO"):
        self.spec_id = spec_id          # ID of input spectrum
        self.spec_f = spec_f            #input spectrum
        self.spec_err = spec_err        #1s uncertainty of input spectrum
        self.spec_SNR = spec_f/spec_err #SNR of the spectrum
        self.spec_lambda = spec_lambda  #lambda array for input spectrum
        
        self.fit_residuals = np.zeros_like(spec_f)*np.nan  #plot residuals from the specplot
        self.template_f = np.zeros_like(spec_f)      # artificially created spectrum after fit with catalog wavelengths.
        self.fit_f = np.zeros_like(spec_f)      # fitted spectrum.
        self.continuum = np.zeros_like(spec_f)     # the continuum of the artificially created spectrum.
        
        self.cat = None # pandas dataframe for all the line fit parameters
        self.rv = None      #radial velocity of the star
        self.erv = None     #1sigma radial velocity uncertainty
        
        self.rv_peak = None      #radial velocity of the star
        self.erv_peak = None     #1sigma radial velocity uncertainty
        
        self.loglevel = loglevel
        
        self.logger = logging.getLogger(self.spec_id)
        
        if self.loglevel == "INFO": self.logger.setLevel(logging.INFO)
        if self.loglevel == "DEBUG": self.logger.setLevel(logging.DEBUG)
        
        fh = logging.FileHandler(str(self.spec_id)+'.log',mode='w')
        
        fh.setLevel(self.loglevel)
        formatter = logging.Formatter('%(asctime)s [%(levelname)8s]: %(message)s [%(funcName)s]',datefmt="%Y-%m-%d %H:%M:%S")
        fh.setFormatter(formatter)        
        
        self.logger.addHandler(fh)
        self.logger.info('Initiate instance of RV_spectrum for: '+str(self.spec_id))
        self.logger.debug('DEBUG mode => all modules run on one core')
        
        
    def clean(self):
        
        """
        this cleans the output in case one wants to repeat the spectral fit
        """
        
        self.logger.info('The fit output is cleaned:')
        
        self.fit_residuals = np.zeros_like(self.spec_f)
        self.template_f = np.zeros_like(self.spec_f)
        self.continuum = np.zeros_like(self.spec_f)
    
    def catalog(self, initcat = None, save = False, load = None, printcat = False):
        """
        Initializes the catalog
        
        Args:
            line_file (str): reading in the line catalog of the ascii format name,lambda,start,end
        """
        
        if load:
            self.logger.info('Load catalog from file: '+str(load))
            if initcat:
                initcat = None
                self.logger.warning('Autoset initcat = None when loading a catalog ')
                
            self.cat.read_csv(path_or_buf = load)
            
        if save:
            self.logger.info('Save catalog to file: '+str(load))
            self.cat.to_csv(path_or_buf = self.spec_id+'.cat')
            
        if printcat: print(self.cat)
        
        if initcat:
            self.logger.info('Initiating the catalog')
        
            temp = ascii.read(initcat)
        
            self.logger.info('Adding lines: '+str(temp['name'].data))

            d = {'l_lab': temp['lambda'],'l_start': temp['start'],'l_end': temp['end'],\
            'l_fit':np.empty_like(temp['lambda']),'a_fit':np.empty_like(temp['lambda']),\
            'sg_fit':np.empty_like(temp['lambda']),'sl_fit':np.empty_like(temp['lambda']),\
            'cont_order':temp['contorder'],'RV':np.empty_like(temp['lambda']),\
            'eRV':np.empty_like(temp['lambda']),'used':np.empty_like(temp['name']),\
            'significance':np.empty_like(temp['lambda'])}
        
            self.cat = pd.DataFrame(d,index = temp['name'].data)
        
    def plot(self):
        
        self.logger.info('Initiate plotting')
        
        n_lines = len(self.cat.index)
        
        nrows = int(len(self.cat.index)/3.) + 1
        if len(self.cat.index)%3 > 0: nrows += 1
        
        plthight = 2 * nrows + 1
        
        spec_plot = plt.figure(self.spec_id,figsize = (8,plthight))  
        plotgrid = gridspec.GridSpec(nrows,3)
        ax_spec = plt.subplot(plotgrid[0,:])
        
        
        non_inf_cont  = ~np.isinf(self.template_f/self.continuum)
        
        exponent=int(np.log10(np.nanmedian(self.spec_f[non_inf_cont]))-1.)  ## to better plot larger numbers
        factor=float(10**(-exponent))
        
        ax_spec.plot(self.spec_lambda[non_inf_cont],self.spec_f[non_inf_cont]*factor,c='black',label='data',linewidth=2)
        ax_spec.plot(self.spec_lambda[non_inf_cont],self.template_f[non_inf_cont]*factor,c='red',label='template',linewidth=2)
        ax_spec.plot(self.spec_lambda[non_inf_cont],self.fit_residuals[non_inf_cont]*factor,c='green',label='residual',linewidth=2)
        ax_spec.legend(fontsize=12)
        ax_spec.set_ylabel(r'flux [$10^{'+str(exponent)+r'}$ erg/s/cm$^2$/${\rm \AA}$]', fontsize=12)
        ax_spec.set_xlabel(r'wavelength $[{\rm \AA}]$', fontsize=12)
        ax_spec.tick_params(axis='both', which='major', labelsize=12, width = 2)
        ax_spec.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        
        for n in np.arange(len(self.cat.index)):
        
            nrow = int(n/3) + 1
            ncol = n%3
            el = self.cat.index[n]
            ax_norm = plt.subplot(plotgrid[nrow,ncol])
        
            ax_norm.plot(self.spec_lambda,self.spec_f/self.continuum,c='black',label ='data',linewidth=2)
            ax_norm.plot(self.spec_lambda,self.template_f/self.continuum,c='red',label ='template',linewidth=2)
            ax_norm.plot(self.spec_lambda,self.fit_f/self.continuum,'--',c='lightblue',label ='fit',linewidth=2)
            ax_norm.plot(self.spec_lambda,self.fit_residuals/self.continuum,c='green',label='residual',linewidth=2)
            ax_norm.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax_norm.yaxis.set_major_locator(plt.MaxNLocator(3))
            ax_norm.xaxis.set_major_locator(plt.MaxNLocator(2))
            
                
            ax_norm.set_xlabel(r'wavelength $[{\rm \AA}]$', fontsize=12)
            ax_norm.tick_params(axis='both', which='major', labelsize=12, width = 2)
            ax_norm.tick_params(axis='y', which='major', labelsize=12, width = 2,direction='in', left=True, right=True)
            if ncol == 0: ax_norm.set_ylabel(r'norm. flux', fontsize=12)
            if ncol != 0: ax_norm.yaxis.set_ticklabels([])
                
                
            ax_norm.set_xlim(self.cat.loc[el,'l_start'],self.cat.loc[el,'l_end'])
            ax_norm.annotate(el,(0.1,0.5),xycoords = 'axes fraction',fontsize = 12,fontweight = 'bold')
            
        
        plt.tight_layout(w_pad=-0.25,h_pad=-1.8)
        plt.savefig(self.spec_id+'_plot.png',dpi=600,overwrite=True)
        plt.close()

    def line_fitting(self,input_cat,line_idxs,niter = 5, contorder = [1], n_CPU = -1, resid_level = None, max_contorder = 2,max_ladjust = 4, adjust_preference = 'contorder',input_continuum_deviation = 0.05, llimits = [2.,2.], max_exclusion_level = 0.3):
        
        """
        Initializes the line fitting
        
        Args:
            input_cat (:array:'float'): input spectral line catalog wavelengths only in Angstrom
            line_idxs (:array:'str'): str array of the linen names for which RV fits 
                should be performed, must be identical to "name" in init_cat
            niter (int, optional, default=5): maximum number of iterations for the line fitting
            contorder (:array:'int, optional, default=[1]): Array of the continuum orders per line fit.
        """
        
        self.logger.info('Starting the line fitting')
        if self.loglevel == "DEBUG": n_CPU = 1
        
        
        start_time  = time.time()
        
        warnings.filterwarnings('ignore',category=RuntimeWarning,message='divide by zero encountered in true_divide')
        warnings.filterwarnings('ignore',category=RuntimeWarning,message='invalid value encountered in true_divide')
        warnings.filterwarnings('ignore',category=RuntimeWarning,message='invalid value encountered in greater_equal')
        
        if n_CPU == -1: n_CPU = cpu_count()
        self.logger.info('Max number of cores: '+str(n_CPU))
        
        if len(contorder) == 1: contorder = np.full_like(line_idxs,fill_value=contorder,dtype=int)
        for ix, el in enumerate(line_idxs): self.cat.loc[el,'cont_order'] = contorder[ix]
        
        result = [dask.delayed(line_fitter)(self,input_cat,ldx,niter,resid_level,max_contorder,max_ladjust,\
                  adjust_preference,input_continuum_deviation,llimits,max_exclusion_level) for ldx in np.array(line_idxs)]
        
        if self.loglevel == "DEBUG": results = dask.compute(*result,num_worker=8,scheduler='single-threaded')
        if not self.loglevel == "DEBUG": results = dask.compute(*result,num_worker=len(line_idxs),scheduler='processes')
        
        
        for result in results:
            
            if not result[13]:
                self.template_f[result[5]] = result[6][result[5]]
                self.fit_f[result[5]] = result[11][result[5]]
                self.cat.loc[result[0],'l_fit'] = result[1]
                self.cat.loc[result[0],'a_fit'] = result[2]
                self.cat.loc[result[0],'sl_fit'] = result[3]
                self.cat.loc[result[0],'sg_fit'] = result[4]

                self.continuum[result[5]] = result[7]
                self.template_f[result[5]] += result[7]
                self.fit_f[result[5]] += result[7]
                self.fit_residuals[result[5]] = self.spec_f[result[5]] - self.fit_f[result[5]]
        
                self.cat.loc[result[0],'l_start'] = result[8]
                self.cat.loc[result[0],'l_end'] = result[9]
                self.cat.loc[result[0],'cont_order'] = result[10]
            
                self.cat.loc[result[0],'significance'] = result[12]
        
                if self.cat.loc[result[0], 'l_lab']\
                - self.cat.loc[result[0], 'l_fit'] > 0.8 * llimits[0]\
                or self.cat.loc[result[0],'l_fit']\
                - self.cat.loc[result[0], 'l_lab'] > 0.8 * llimits[1]:
                    self.logger.warning(result[0]\
                    + ' exceeds 0.8 of lambda limits; lfit: '\
                    + str(self.cat.loc[result[0], 'l_fit'])\
                    + ' llab: '+str(self.cat.loc[result[0], 'l_lab'])\
                    + ' llimits: '+str(llimits[0])+' '+str(llimits[1]))
            
            if result[13]:
                self.logger.warning(result[0]+' failed and will be marked in the catalog')
                self.cat.loc[result[0],'used'] = 'f'
            
        
        elapsed_time = time.time() - start_time
        self.logger.info('Finished line fitting')
        self.logger.info('Elapsed time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    def rv_fit_peak(self,line_sigma = 3,line_significants = 5):
        
        self.logger.info('Calculating RVs based on peaks only')
        v = np.array((self.cat.loc[:,'l_fit'].values/self.cat.loc[:,'l_lab'] - 1.)*const.c.to('km/s').value)
        
        remaining_lines = line_clipping(self,v,line_significants,sigma = line_sigma)
        
        if remaining_lines.mask.all():
            self.logger.error('NO USABLE LINE FOUND WITH SET PARAMETER !!')
        else:
            self.rv_peak = np.median(v[~remaining_lines.mask])
            self.erv_peak = MAD(v[~remaining_lines.mask])
            
            self.logger.info('Finished Calculating RVs based on peaks only: RV=('+str('{:4.2f}'.format(self.rv_peak))+\
                             '+-'+str('{:4.3f}'.format(self.erv_peak))+')km/s based on '+str(len(v[~remaining_lines.mask]))+'/'+str(len(v))+' lines')

    def rv_fit(self, guesses, niter = 10000, line_sigma = 3, n_CPU = -1, line_significants = 5):
        
        """
        Initializes the RV fit
        
        Args:
            guesses (:array:'float') :array[vel,sig], initial guesses for the RV fit
            niter (int, optional, default = 10000): number of iterations for resampling the spectrum
            line_exclude (:list:'str', optional, default=False): line names that should be manually excluded
            line_sigma (float, optional, default=3): sigma for the RV clipping for the individual lines        
            n_CPU (int, optiona, default=-1): number of CPUs used. -1 for all available
        
        """
        
        self.logger.info('Starting the RV fitting')
        self.logger.info('Settings: niter='+str(niter)+'; sigma RV for excluding lines: '+str(line_sigma))
        
        start_time  = time.time()
        
        if self.loglevel == "DEBUG": n_CPU = 1
        if n_CPU == -1: n_CPU = cpu_count()
        self.logger.info('Max number of cores: '+str(n_CPU))
        
        ### transfer the spectrum to log-space
        log_spec_f, logspec_lambda, velscale_spec = ppxf_util.log_rebin([self.spec_lambda[0],self.spec_lambda[-1]],self.spec_f)
        log_template_f, log_template_lambda, velscale_template = ppxf_util.log_rebin([self.spec_lambda[0],self.spec_lambda[-1]],self.template_f)
        log_spec_err, log_spec_err_lambda, velscale_spec_err = ppxf_util.log_rebin([self.spec_lambda[0],self.spec_lambda[-1]],self.spec_err)
        
        log_spec_err[~np.isfinite(log_spec_err)] = 1. # This happens if AO was used and there is a gap in the spec
        
        exponent=int(np.log10(np.nanmedian(log_spec_f))-4.) # Obtain larger flux numbers to make it numerically stable
        factor=float(10**(-exponent))
 
 
        log_spec_f = log_spec_f*factor
        log_template_f = log_template_f*factor
        log_spec_err = log_spec_err*factor
        
        l_fit = self.cat.loc[:,'l_fit'].values.astype(np.float64)
        l_lab = self.cat.loc[:,'l_lab'].values.astype(np.float64)
        line_name = self.cat.index
        
        fwhm_g, fwhm_l, fwhm_v = voigt_FWHM(self.cat.loc[:,'sg_fit'].values.astype(np.float64),self.cat.loc[:,'sl_fit'].values.astype(np.float64))
        
        l_fit = np.delete(l_fit,np.where(self.cat.loc[:,'used'] == 'f'))
        l_lab = np.delete(l_lab,np.where(self.cat.loc[:,'used'] == 'f'))
        line_name = np.delete(line_name,np.where(self.cat.loc[:,'used'] == 'f'))
        fwhm_g = np.delete(fwhm_g,np.where(self.cat.loc[:,'used'] == 'f'))
        fwhm_l = np.delete(fwhm_l,np.where(self.cat.loc[:,'used'] == 'f'))
        fwhm_v = np.delete(fwhm_v,np.where(self.cat.loc[:,'used'] == 'f'))
        
        v = np.zeros_like(l_lab)
        ev = np.zeros_like(l_lab)
        
        for i,line in enumerate(l_lab):
            
            self.logger.info('Started line '+line_name[i])
            
            mask = np.zeros(len(log_template_f), dtype=bool)
            ind_min = np.argmin(np.abs(log_template_lambda-np.log(line-0.5*fwhm_v[i])))-1
            ind_max = np.argmin(np.abs(log_template_lambda-np.log(line+0.5*fwhm_v[i])))+1
            ind_center=np.argmin(np.abs(log_template_lambda-np.log(line)))
            mask[ind_min:ind_max+1] = True

            pp_outliers_init=ppxf.ppxf(log_template_f,log_spec_f,log_spec_err,velscale_spec,guesses,mask = mask,degree=-1,clean=False,quiet=True,plot=False)
            
            v[i], ev[i] = ppxf_MC(log_template_f,log_spec_f,log_spec_err,velscale_spec, guesses,nrand= 0.5*niter,goodpixels=pp_outliers_init.goodpixels, degree=-1, moments=2, n_CPU = n_CPU)
            
            self.cat.loc[line_name[i],'RV'] = v[i]
            self.cat.loc[line_name[i],'eRV'] = ev[i]
            
            self.logger.info('Finished line '+line_name[i]+': RV=('+str('{:4.2f}'.format(v[i]))+'+-'+str('{:4.3f}'.format(ev[i]))+')km/s')
    
        remaining_lines = line_clipping(self,v,line_significants,sigma = line_sigma)

        if remaining_lines.mask.all():
            self.logger.error('NO USABLE LINE FOUND WITH SET PARAMETER !!')
        else:
            self.cat.loc[line_name[~remaining_lines.mask],'used'] = 'x'

            l_fit = l_fit[~remaining_lines.mask]
            fwhm_v = fwhm_v[~remaining_lines.mask]

            mask=np.zeros(len(self.spec_lambda),dtype=bool)
        
            for i,line in enumerate(l_fit):
                ind_min = np.argmin(np.abs(logspec_lambda-np.log(line-0.5*fwhm_v[i])))-1
                ind_max = np.argmin(np.abs(logspec_lambda-np.log(line+0.5*fwhm_v[i])))+1
                ind_center = np.argmin(np.abs(logspec_lambda-np.log(line)))
                mask[ind_min:ind_max+1] = True
        
            pp_final_plot = plt.figure(str(self.spec_id)+'_ppxf_fit_final',figsize=(10,3))
            if len(l_lab) > 1: pp_final_init = ppxf.ppxf(log_template_f,log_spec_f,log_spec_err,velscale_spec,guesses,degree=-1,clean=False,mask=mask,quiet=True)
            if len(l_lab) == 1: pp_final_init = pp_outliers_init
            pp_final_init.plot()
            plt.tight_layout()
            plt.savefig(self.spec_id+'_ppxf_fit_final.png',dpi=600,overwrite=True)
            plt.close()
        
            self.logger.info('Started final RV fit')
        
            self.rv, self.erv = ppxf_MC(log_template_f,log_spec_f,log_spec_err,velscale_spec,guesses,\
                                        nrand=niter,goodpixels=pp_final_init.goodpixels,degree=-1,moments=2,spec_id = self.spec_id, n_CPU = n_CPU)

            elapsed_time = time.time() - start_time
            self.logger.info('Used lines for RV fit: '+', '.join(line_name[~remaining_lines.mask]))
            self.logger.info('Finished RV fitting: RV=('+str('{:4.2f}'.format(self.rv))+'+-'+str('{:4.3f}'.format(self.erv))+')km/s based on '+str(len(l_fit))+'/'+str(len(line_name))+' lines')
            self.logger.info('Elapsed time: '+time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

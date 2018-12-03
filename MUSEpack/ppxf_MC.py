#!/usr/bin/env python
"""
ppxf_MC.py

Copyright 2018-2018 Peter Zeidler

This module is part of MUSEpack and adds and subtracts the spectral uncertainties in a radom fashion
and then for each random draw the routine recalculates RV.
The final result is the mean of the gaussian dristibution and its uncertanty is
the 1sigma width.

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

"""
import numpy as np
from astropy.stats import sigma_clip
from lmfit import Model
from ppxf import ppxf
from multiprocessing import Pool, cpu_count
from functools import partial
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation as MAD
import logging

def gaussian_fit(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))

def ppxf_bootstrap(log_template_f, log_spec_f, log_spec_err, velscale, degree, goodpixels, guesses,moments, vsyst,iter_spec,noise,k):

    # Use negligible BIAS=0 to estimate Bootstrap errors
    # See Section 3.4 of Cappellari & Emsellem (2004)
    
    ### add or subtract the uncertainties
    
    pm = np.random.randint(1,3,size=len(log_spec_err))
    pm[np.where(pm == 2)] = -1
    pm_log_spec_err = pm*log_spec_err
    
    #random draw of uncertainties and add them
    np.random.seed()
    scramble = goodpixels[(np.random.rand(len(goodpixels))*len(goodpixels)).astype(int)]
    new_log_spec_err = pm_log_spec_err[scramble]
    iter_spec[goodpixels] = log_spec_f[goodpixels] + new_log_spec_err

    pp1 = ppxf.ppxf(log_template_f, iter_spec, noise, velscale, guesses, goodpixels=goodpixels,
                    plot=False, degree=degree, moments=moments, vsyst=vsyst, quiet=1,bias=0,velscale_ratio=1) 

    if moments <= 2: uncert=np.array([pp1.sol[0],pp1.sol[1]])
    if moments > 2: uncert=np.array([pp1.sol[0],pp1.sol[1],pp1.sol[2],pp1.sol[3]])
    
    return uncert

def ppxf_MC(log_template_f, log_spec_f, log_spec_err, velscale, guesses, nrand=100,
                    degree=4, goodpixels=None, moments=4, vsyst=0, sigma=5, spec_id = None,n_CPU = -1):

    '''
    log_template_f: array, logarithmically binned template spectrum
    log_spec_f:     array, logarithmically binned spectrum
    log_spec_err:   array, logarithmically binned uncertainties
    velscale:       float, velocity scale

    THE ABOVE DEFINED VALUES CAN BE CALCULAYED WITH log_rebin OF THE ppxf_util package

    degree         float, optional, default =4, degree of the additive continuum
    goodpixels     array, optional, bool, default: None; same dimension as log_template_f, array of spectral bins used for fit
    moments=4      int, optional, default = 4, moments to calculate velcoities
    vsyst=0        flat, optional, default: 0, velocity of the system

    THE ABOVE PARAMETERS ARE ALSO DEFINED AND FURTHER EXPLAINED IN Cappellari & Emsellem (2004)

    guesses:        array[velocity,1sigma], starting value for RV fit
    nrand:          float, default: 100, number of iterations for MC
    sigma:          float, optional, default: 5, sigma for sigma clipping of velocity results for plotting the gaussian
    spec_id:        str, optional, default: -1, identifier of the spectrum to name the plot, if plot is named it is also safed
    n_CPU:          int, optional, default: None, number of cores to use for calculation                
                    
    '''
    noise = np.ones_like(log_spec_f) # Constant error spectrum
    if goodpixels.any() == None: goodpixels = np.arange(len(log_spec_f))
    iter_spec = log_spec_f.copy()

    if n_CPU > 0: cpu = n_CPU
    if n_CPU ==-1 : cpu = cpu_count()
    
    pool=Pool(cpu)
    uncert_ppxf=pool.map(partial(ppxf_bootstrap,log_template_f, log_spec_f, log_spec_err, velscale,
                                 degree, goodpixels, guesses, moments, vsyst,iter_spec,noise),np.arange(nrand))

    pool.close()
    pool.join()
   
    clipped=sigma_clip(np.array(uncert_ppxf)[:,0], sigma=sigma, stdfunc=MAD)
    clipped=clipped.data[~clipped.mask]
    
    
    if not spec_id:
        ret = np.histogram(clipped,bins=int(20),density=True)
        bins = [ret[1][i]+0.5*(ret[1][i+1]-ret[1][i]) for i in range(len(ret[0]))]
        
        result=Model(gaussian_fit).fit(ret[0], x=bins, a=np.max(ret[0]), b=np.mean(clipped), c=np.std(clipped))
        popt=[result.best_values['a'],result.best_values['b'],result.best_values['c']]
    
    if spec_id:
        plt.figure(spec_id)
        ret = plt.hist(clipped, bins=int(20),density=True,label='n realizations: '+str(nrand))
        bins = [ret[1][i]+0.5*(ret[1][i+1]-ret[1][i]) for i in range(len(ret[0]))]
        
        binlength=bins[1]-bins[0]
        
        result=Model(gaussian_fit).fit(ret[0], x=bins, a=np.max(ret[0]), b=np.mean(clipped), c=np.std(clipped))
        popt=[result.best_values['a'],result.best_values['b'],result.best_values['c']]
        
        x=((np.arange(20000*len(bins))-10000*len(bins))/(100*len(bins)))+popt[1]
        
        plt.plot(x,gaussian_fit(x,*popt),lw=3)
        plt.axvline(x=popt[1],c='r',lw=3,
                    label=r'mean$=$'+str('{:4.2f}'.format(popt[1]))+r'$\,\frac{\rm km}{\rm s}$')
        plt.axvline(x=popt[1]+popt[2],c='r',linestyle='--',lw=3,
                    label=r'$1\sigma = $'+str('{:5.3f}'.format(popt[2]))+r'$\,\frac{\rm km}{\rm s}$')
        plt.axvline(x=popt[1]-popt[2],c='r',linestyle='--',lw=3)
        
        plt.xlim(popt[1]-3*popt[2],popt[1]+3*popt[2])
        plt.xlabel(r'velocity [$\frac{\rm km}{\rm s}$]')
        plt.ylabel(r'relative number')
        plt.legend()
        plt.tight_layout()
        plt.savefig(spec_id+'_v_dist.png',dpi=600)
        plt.close()
     
    return popt[1], popt[2]
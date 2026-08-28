 #!/usr/bin/env python


__version__ = '1.0.0'

__revision__ = '20260814'

import sys
import os
import gc
import shutil
import warnings
import pyspeckit
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column
from astropy import units as u
from astropy import constants as const
from astropy import log
from astropy.stats import median_absolute_deviation as MAD
import pandas as pd
import time
import logging
import stsynphot as stsyn
import synphot as syn
from specutils import Spectrum1D
from astropy.coordinates import SkyCoord
from pathlib import Path

from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.special import wofz
from scipy.interpolate import UnivariateSpline

from itertools import combinations as inter_comb

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import matplotlib.gridspec as gridspec

import dask
from dask.distributed import Client, progress, wait

def initial_guesses(self, lines, blends=None, linestrength=100.,\
                    llimits=[2. * (-1), 2.]):

    '''
    Creates the initial guesses for the line fitter

    Args:
        lines : :func:`numpy.array`
            central wavelengths of the spectral lines

    Kwargs:
        linestrength : :obj:`float` (default: 100)
            initial guess for the line strength

        blends : :obj:`str` (optional)
            A file containing the a list of blended lines in the
            format: ** List is coming soon**

    return:
        guesses : :obj:`list`
            lists containing the guesses

        limits : :obj:`list`
            lists containing the limits

        limited : :obj:`list`
            lists containing the limited
    '''

    if self.linetype == 'absorption':
        linestrength = linestrength * (-1)

    guesses = np.zeros(4 * len(lines))
    limits = []
    limited = []

    for idx, line in enumerate(lines):

        if self.rv_sys != 0.:
            line = lambda_rv_shift(self, line)

        guesses[4 * idx + 0] = linestrength
        guesses[4 * idx + 1] = line
        guesses[4 * idx + 2] = self.dispersion
        guesses[4 * idx + 3] = self.dispersion

        limits.append((0, 0)) #Amplitude
        if self.linetype == 'absorption':
            limited.append((False, True)) #Amplitude
        if self.linetype == 'emission':
            limited.append((True, False)) #Amplitude
        if self.linetype == 'both':
            limited.append((False, False)) #Amplitude

        limits.append((line + llimits[0], line + llimits[1]))
        limited.append((True, True)) #Center

        limits.append((0, 0)) #width (gausssigma)
        limited.append((True, False)) #width

        limits.append((0, 0)) #width (lorentzsigma)
        limited.append((True, False)) #width

    if blends:
        blendlist = ascii.read(blends)

        for blendidx in range(len(blendlist)):

            # lprime = blendlist[blendidx]['lprime']
            # lsec = blendlist[blendidx]['lsec']

            if self.rv_sys == 0.:
                lprime = blendlist[blendidx]['lprime']
                lsec = blendlist[blendidx]['lsec']
            else:
                lprime = lambda_rv_shift(self,
                blendlist[blendidx]['lprime'])
                lsec = lambda_rv_shift(self,
                blendlist[blendidx]['lsec'])

            primeidx = np.where(lprime == guesses)[0]
            secidx = np.where(lsec == guesses)[0]

            if len(primeidx) == 1:

                self.logger.info('Blend detected: '\
                + str(blendlist[blendidx]['prime'])\
                + ' and ' + str(blendlist[blendidx]['sec']) + ' with ratio '\
                + str(blendlist[blendidx]['ratio']))

                guesses[secidx - 1] = blendlist[blendidx]['ratio']\
                * guesses[primeidx - 1]

                if lprime < lsec and (lprime + llimits[1]) > lsec + llimits[0]:
                    limits[secidx[0]]\
                    = (limits[primeidx[0]][1] + 0.01, limits[secidx[0]][1])

                    if limits[primeidx[0]][1] + 0.01 >= lsec:
                        limits[secidx[0]] = (lsec - 0.5, limits[secidx[0]][1])

                if lprime > lsec and (lprime + llimits[0]) < lsec + llimits[1]:
                    limits[secidx[0]] = (limits[secidx[0]][0],\
                                         limits[primeidx[0]][0] - 0.01)
                    if limits[primeidx[0]][0] - 0.01 <= lsec:
                        limits[secidx[0]] = (limits[secidx[0]][0], lsec + 0.5)

    return guesses, limits, limited


def update_parinfo(self, guesses, llimits, line_idx, blends,
                   parinfo, autoadjust, fwhm_block):

    r"""
    Updates the parinfo file, created by pyspeckit.

    Args:
        guesses : :func:`numpy.array`
            The initial guesses for the the radial velocity fit guesses in
            the form [RV,sepctral_dispersion]

        llimits : :obj:`list`
            the limits for the wavelength fit as set in ``ppxf``

        line_idx : :obj:`str`
            Name of the primary line

        blends : :obj:`ascii`-file or :obj:`None`
           A file with primary lines that contain blends to provide a maximum
           amplitude ratio of the primary and the blend to prevent that the
           blend becomes the dominant line in the fit.

        parinfo: :obj:`dict`
            the parinfo file created by pyspeckit, which contains the fitted
            parameters for all input lines

        autoadjust : :obj:`bool`
            :obj:`True`: the wavelength limits ``llimit`` will be adjusted to
            the fit of the previous iteration. All other wavelength range are
            adjusted accordingly taking into account the proper velocity
            corrected shift :math:`\Delta \lambda/\lambda`. This is especially
            important to detect hyper-velocity stars.

            :obj:`False`: no adjustment to the limits done

        fwhm_block : :obj:`bool:obj:`
            :obj:`True`: The minimum fwhm of the voigt profiles of the fitted
            lines is the instrument's dispersion

            :obj:`False`: The minimum fwhm of the voigt profiles of the fitted
            lines is zero

    """

    if self.rv_sys == 0.:
        lprime = self.cat.loc[line_idx, 'l_lab']
        primeidx = np.where(lprime == guesses)[0]
    else:
        lprime = lambda_rv_shift(self, self.cat.loc[line_idx, 'l_lab'])
        primeidx = np.where(lprime == guesses)[0]

    if autoadjust:
        lshift_in = parinfo[int(primeidx[0])]['value']\
        - guesses[int(primeidx[0])]

        for adj_idx in range(int(len(guesses) / 4.)):

            lshift = _lambda_shift(lshift_in, guesses[int(primeidx[0])],\
            guesses[int(4 * adj_idx + 1)])

            parinfo[int(4 * adj_idx + 1)]['limits']\
            = (guesses[int(4 * adj_idx + 1)] + llimits[0]\
            + lshift, guesses[int(4 * adj_idx + 1)]\
            + llimits[1] + lshift)

            parinfo[int(4 * adj_idx + 1)]['value']\
            = guesses[int(4 * adj_idx + 1)] + lshift

    #### check that fwhm is always > dispersion:

    if fwhm_block:
        for paridx in range(int(len(guesses) / 4.)):

            fwhm_g, fwhm_l, fwhm_v\
            = voigt_FWHM(np.array([parinfo[int(4 * paridx + 2)]['value']],\
            dtype=float), np.array([parinfo[int(4 * paridx + 3)]['value']],\
            dtype=float))
            if fwhm_v[0] < self.dispersion:
                target_fwhm_g = self.dispersion - fwhm_l[0]

                s = target_fwhm_g / 2. / np.sqrt(2. * np.log(2.))
                parinfo[int(4 * paridx + 2)]['limits']\
                = (s, parinfo[int(4 * paridx + 2)]['limits'][1])
                parinfo[int(4 * paridx + 2)]['value'] = s

    if blends:
        blendlist = ascii.read(blends)

        for blendidx in range(len(blendlist)):
            blendration = blendlist[blendidx]['ratio']
            lsec = blendlist[blendidx]['lsec']
            secidx = np.where(lsec == guesses)[0]
            if len(secidx) == 1:
                prime_profile = voigt_funct(self.spec_lambda_highres,\
                                    parinfo[int(primeidx[0])]['value'],\
                                    parinfo[int(primeidx[0] - 1)]['value'],\
                                    parinfo[int(primeidx[0] + 1)]['value'],\
                                    parinfo[int(primeidx[0] + 2)]['value'])

                sec_profile = voigt_funct(self.spec_lambda_highres,\
                                    parinfo[int(secidx[0])]['value'],\
                                    parinfo[int(secidx[0] - 1)]['value'],\
                                    parinfo[int(secidx[0] + 1)]['value'],\
                                    parinfo[int(secidx[0] + 2)]['value'])

                ind_max_prime = np.argmax(np.abs(prime_profile))
                ind_max_sec = np.argmax(np.abs(sec_profile))

                if np.max(np.abs(sec_profile)) >= \
                blendration * np.max(np.abs(prime_profile)):

                # if parinfo[int(secidx[0] - 1)]['value'] >= \
                # blendration * parinfo[int(primeidx[0] - 1)]['value']:

                    parinfo[int(secidx[0] - 1)]['limited'] = (True, True)

                    parinfo[int(secidx[0] - 1)]['limits'] \
                    = (blendration * parinfo[int(primeidx[0] - 1)]['value'], \
                    parinfo[int(secidx[0] - 1)]['limits'][1])

                    parinfo[int(secidx[0] - 1)]['value']\
                    = 0.99 * blendration * prime_profile[ind_max_prime]
                    # = 0.99 * blendration\
                    # * parinfo[int(primeidx[0] - 1)]['value']

                if parinfo[int(primeidx[0] - 1)]['value'] == 0.:
                    parinfo[int(secidx[0] - 1)]['limits'] \
                    = (0.01 * (-1), \
                    parinfo[int(secidx[0] - 1)]['limits'][1])
                    parinfo[int(secidx[0] - 1)]['value']\
                    = blendration * (0.01) * (-1.)

    return parinfo


def plotcolor(n):
    colarray = np.array(['Lime', 'Blue', 'Magenta', 'Olive', 'Maroon',
    'indigo', 'orange', 'Cyan', 'Yellow', 'Silver'])
    colselect = colarray[n % len(colarray)]
    return colselect


def baseline_extract(self, specfit, poly_order):
    coeff = np.zeros(poly_order + 1)
    power = np.arange(poly_order + 1)
    for order in range(poly_order + 1):
        coeff[order] = specfit.header['BLCOEF' + str('{:02d}'.format(order))]
    baseline = np.poly1d(coeff / (self.specbinsize ** power[::-1]))
    return baseline


def voigt_FWHM(sigma_g, gamma_l):
    #calculating the FWHM of the voigt, gauss, and lorentz profile
    c_0 = 2.0056
    c_1 = 1.0593
    fwhm_g = 2. * np.sqrt(2. * np.log(2.)) * sigma_g
    fwhm_l = 2. * gamma_l
    phi = fwhm_l / fwhm_g
    fwhm_v = np.zeros_like(sigma_g, dtype=float)
    for i in range(len(fwhm_g)):
        if fwhm_g[i] == 0:
            fwhm_v[i] = fwhm_l[i]
        else:
            fwhm_v[i] = fwhm_g[i] * (1. - c_0 * c_1 + np.sqrt(phi[i] ** 2\
            + 2. * c_1 * phi[i] + (c_0 * c_1) ** 2))

    return fwhm_g, fwhm_l, fwhm_v


def voigt_funct(x_array, x_cen, amplitude, sigma, gamma):

    if sigma > 0 and gamma > 0:
        z = ((x_array - x_cen) + 1j * gamma) / (sigma * np.sqrt(2.))
        y_arr = amplitude * np.real(wofz(z))
    if sigma == 0:
        y_arr = _lorentzian(x_array, x_cen, amplitude, gamma)
    if gamma == 0:
        y_arr = _gaussian(x_array, x_cen, amplitude, sigma)

    return y_arr


def spec_res_downgrade(l_in, spec_in, l_out):
    templ_spec = syn.SourceSpectrum(syn.Empirical1D, points=l_in, lookup_table=spec_in, keep_neg=True)
    white_filter = syn.SpectralElement(syn.models.Empirical1D, points=templ_spec.waveset,
                                      lookup_table=np.ones(len(templ_spec.waveset)), keep_neg=True)
    obs_convolved = syn.Observation(templ_spec, white_filter, force='taper', binset=l_out)
    convolved_spec = obs_convolved.as_spectrum(binned=True)
    #
    # convolved_spec = Spectrum1D(spectral_axis=obs_convolved.waveset,
    #                             flux=syn.units.convert_flux(obs_convolved.waveset,
    #                                               fluxes=obs_convolved(obs_convolved.waveset)))
    return convolved_spec(convolved_spec.waveset).value


def line_clipping(self, x, line_significants, sigma=3):

    mask = np.zeros_like(x, dtype=int)
    ind_sig = np.where((self.cat.loc[:,\
    'significance'].values.astype(np.float64) < line_significants)\
                       | (self.cat.loc[:, 'used'].values == 'f'))
    mask[ind_sig] = 1

    x_red1 = np.delete(x, np.where(mask == 1))

    if len(x_red1) >= 3:
        if len(x_red1) == 3:
            median = np.median(x_red1)
            mad = MAD(x_red1)
            ind = np.where((x < median - sigma * mad) | (x > median\
            + sigma * mad))[0]
        if len(x) > 3:
            gr = np.array(list(inter_comb(x_red1, len(x_red1) - 2)))

            std_gr = [np.std(it) for it in gr]
            x_red = gr[np.argmin(std_gr)]

            median = np.median(x_red)
            mad = MAD(x_red)
            ind = np.where((x < median - sigma * mad) | (x > median\
            + sigma * mad))[0]

        mask[ind] = 1
    x_masked = np.ma.masked_array(x, mask=[mask])
    return x_masked


def clipped_MAD(x):
    if len(x) <= 3:
        mad = MAD(x)
    if len(x) > 3:
        x = np.delete(x, [np.argmin(x), np.argmax(x)])
        mad = MAD(x)
    return mad


def continuum_deviation(self, l_in, f_in, baseline, contorder):

    if self.linetype == 'absorption':
        peaks = find_peaks(f_in * (-1), width=2)[0]
    else:
        peaks = find_peaks(f_in, width=2)[0]
    mask = np.zeros_like(f_in)

    for p in peaks:
        mask[p - 3:p + 4] = 1

    remove = ~np.ma.array(l_in, mask=mask).mask

    l_masked = l_in[remove]
    f_masked = f_in[remove]
    baseline_masked = baseline[remove]

    masked_fit = np.polyfit(l_masked, f_masked, contorder)
    masked_fitfunc = np.poly1d(masked_fit)

    cont_dev = MAD(f_masked / np.median(f_masked)
    - baseline_masked / np.median(baseline_masked))

    return cont_dev


def ABtoVega(instrument, bandpass):

    bp = stsyn.band(str(instrument) + ',wfc1,' + str(bandpass) + ',mjd#57754')
    obs = syn.Observation(stsyn.Vega, bp, binset=bp.binset)
    photflam = bp.unit_response(stsyn.conf.area)
    photplam = bp.pivot()

    zp_vega = (-1) * obs.effstim(flux_unit='obmag', area=stsyn.conf.area).value
    zp_ab = (-2.5 * np.log10(photflam.value) - 21.1
             - 5 * np.log10(photplam.value) + 18.6921)
    zp_st = -2.5 * np.log10(photflam.value) - 21.1

    difference = zp_vega - zp_ab

    return difference


def lambda_rv_shift(self, lam):
    lambda_new = lam + self.rv_sys * lam / const.c.to('km/s').value
    if self.rv_sys == 0:
        lambda_new = lam
    return lambda_new


def _lambda_shift(dlin, lin, lout):
    '''
       This functions shifts the wavelength limits if the central wavelength
       is shifted. This assures that lambda never runs out of bounds

    '''

    dlout = dlin * lout / lin
    return dlout


def _gaussian(x, mu, A, sigma):
    prefact = 1. / np.sqrt(2. * np.pi * sigma ** 2)
    exponetial = np.exp((-1.) * ((x - mu) ** 2) / (2. * sigma ** 2))
    returnfunction = A * exponetial
    return returnfunction


def _lorentzian(x, mu, A, gamma):
    returnfunction = A * gamma / np.pi / ((x - mu) ** 2 + gamma ** 2)
    return returnfunction


def _any_bit_in_number(arr, num):
    for elem in arr:
        if elem & num:  # Check if any bit is shared
            return False
    return True


def _make_circle_mask(shape, wcs, ra_deg, dec_deg, radius_arcsec):
    """Boolean mask, True inside a sky-circle of given radius."""
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    x0, y0 = wcs.world_to_pixel_values(coord.ra, coord.dec)

    # Pixel scale in arcsec/pix from the CD matrix (or CDELT)
    pixscale = np.sqrt(np.abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600.0
    radius_pix = radius_arcsec / pixscale
    print(f"  source pixel center: ({x0:.1f}, {y0:.1f})")
    print(f"  pixel scale: {pixscale:.4f} arcsec/pix")
    print(f"  source radius: {radius_arcsec} arcsec = {radius_pix:.1f} pix")

    yy, xx = np.ogrid[:shape[0], :shape[1]]
    r2 = (xx - x0) ** 2 + (yy - y0) ** 2
    return r2 < radius_pix ** 2, (x0, y0, radius_pix)


def _make_rectangle_mask(shape, wcs, ra_deg, dec_deg,
                         pa_deg, length_arcsec, width_arcsec):
    """
    Boolean mask, True inside a rotated rectangle centered on (ra, dec).
    pa_deg is the position angle of the rectangle's long axis,
    measured east of north on the sky.
    """
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    x0, y0 = wcs.world_to_pixel_values(coord.ra, coord.dec)

    pixscale = np.sqrt(np.abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600.0
    L = length_arcsec / pixscale
    W = width_arcsec  / pixscale

    # Determine the on-image angle corresponding to PA on the sky.
    # We compute it numerically: take a small step north (delta-Dec > 0)
    # from the source, see where it lands in pixel coordinates, and that
    # gives us the "north" direction in the image. Then PA is measured
    # east of north.
    step_arcsec = 1.0
    north = SkyCoord(ra=ra_deg * u.deg,
                     dec=(dec_deg + step_arcsec / 3600.0) * u.deg,
                     frame="icrs")
    xn, yn = wcs.world_to_pixel_values(north.ra, north.dec)
    # angle of "north" in image coords (radians, measured from +x axis CCW)
    north_angle = np.arctan2(yn - y0, xn - x0)
    # east is 90 deg CCW from north on the sky; on the image this is
    # north_angle + pi/2 (because we're working in standard math convention)
    pa_rad = np.deg2rad(pa_deg)
    # rectangle long-axis direction in image coords
    angle = north_angle + pa_rad - np.pi / 2.0  # PA east of north -> image angle

    cos_a, sin_a = np.cos(angle), np.sin(angle)

    yy, xx = np.mgrid[:shape[0], :shape[1]]
    dx = xx - x0
    dy = yy - y0
    # Project into the rectangle's local frame
    u_along = dx * cos_a + dy * sin_a   # along the long axis
    v_perp  = -dx * sin_a + dy * cos_a  # perpendicular to it

    inside = (np.abs(u_along) < L / 2.0) & (np.abs(v_perp) < W / 2.0)

    print(f"  bleed center: ({x0:.1f}, {y0:.1f})")
    print(f"  bleed length: {length_arcsec} arcsec = {L:.1f} pix")
    print(f"  bleed width:  {width_arcsec} arcsec = {W:.1f} pix")
    print(f"  bleed image-frame angle: {np.rad2deg(angle):.1f} deg")

    return inside, (x0, y0, L, W, angle)

def _make_diagnostic_plot(img, mask, src_geom_dict, bleed_geom_dict, outfile):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: white-light image with mask regions overlaid
    finite = np.isfinite(img)
    if finite.any():
        vmin, vmax = np.nanpercentile(np.asarray(img[finite]), [5, 99])
    else:
        vmin, vmax = 0, 1
    axes[0].imshow(img, origin="lower", cmap="gray_r", vmin=vmin, vmax=vmax)
    axes[0].set_title("White-light image + mask outlines")

    for geom_dict_key in src_geom_dict.keys():

        x0, y0, R = src_geom_dict[geom_dict_key]

        axes[0].add_patch(Circle((x0, y0), R, fill=False,
                                 edgecolor="red", linewidth=1.5))
        if bleed_geom_dict[geom_dict_key] is not None:
            bx, by, L, W, ang = bleed_geom_dict[geom_dict_key]
            # Rectangle patch is anchored at a corner; rotate around the corner.
            rect = Rectangle(
                (float(bx - L / 2), float(by - W / 2)),
                float(L), float(W),
                angle=float(np.rad2deg(ang)),
                rotation_point=(float(bx), float(by)),
                fill=False, edgecolor="orange", linewidth=1.5,
            )
            axes[0].add_patch(rect)

    # Panel 2: the mask itself
    axes[1].imshow(mask, origin="lower", cmap="gray", vmin=0, vmax=1)
    axes[1].set_title("SKY_MASK (white = sky, black = excluded)")

    for ax in axes:
        ax.set_xlabel("x [pix]")
        ax.set_ylabel("y [pix]")
    plt.tight_layout()
    plt.savefig(outfile, dpi=120)
    print(f"\nWrote {outfile}")

    plt.close(fig)

def _skyfit(refsky, wav, sky, sigma=2., s=1., order=3, plot=True, plot_axis=None):

    ratio = sky/refsky

    median_ratio = median_filter(ratio, size=111, mode='constant', cval=0.0)

    residuals = ratio - median_ratio
    std_dev_residuals = np.std(residuals)

    mask = np.abs(residuals) < (sigma * std_dev_residuals)
    wav_c, ratio_c = wav[mask], ratio[mask]

    ratio_c = median_filter(ratio_c, size=501, mode='reflect')
    ratio_c = savgol_filter(ratio_c, window_length=21, polyorder=2)

    if plot:
        plot_axis.plot(wav_c,ratio_c,color='green', zorder=1, label='smoothed continuum')

    spline = UnivariateSpline(wav_c, ratio_c, s=s, k=order)
    yspline = spline(wav)

    return yspline, spline


def _apply_spline(pixtable, spline_funct_dict, outdir, key):
     pixtable_filename = Path(pixtable).name
     print("   ... Apply continuum correction to:", pixtable_filename)

     with fits.open(pixtable) as hdu:
         data = hdu['data'].data
         wav = hdu['lambda'].data
         spline_funct = spline_funct_dict[key]
         spline = spline_funct(wav)

         print("   ... dividing by spline ...")
         data /= spline
         hdu['data'].data = data

         outfile = os.path.join(outdir, pixtable_filename.replace('.fits','_CONTCORR.fits'))
         print("   ... writing:", outfile)
         hdu.writeto(outfile, overwrite=True)

     del pixtable, spline_funct_dict
     gc.collect()

def _intra_cube_continuum_correction(pixeltables, skycontinuum, pixtab_dir, plot=True, n_CPU=1):

    '''
    This module is intended to correct for continuum differences between individual dithers
    The correction is directly implemented on the PIXTABLES

    Args:
        pixtables : :obj:`list`
                The list of input pixeltables

        skycontinuum : :obj:`list`
                The list of skycontinua that match the pixeltables

        pixtab_dir : :obj:`str`
                The directory where the PIXELTABLES are located

    Kwargs:

        plot : :obj:`bool`
                If set to true plots will be created (default: :obj:`true`)

        n_CPU : :obj:`int`
                The number of cores used for the multicore processing (default: 1)

    '''

    sky_level = {}
    pixtable_dict = {}

    for i, (sky_file, pixtab_filename) in enumerate(zip(skycontinuum, pixeltables)):

        hdu_sky = fits.open(os.path.join(pixtab_dir, sky_file))
        sky_level[pixtab_filename] = hdu_sky[1].data
        pixtable_dict[pixtab_filename] = os.path.join(pixtab_dir, pixtab_filename)

    sky_level_cube = np.vstack([sky_level[pixeltables[0]]['flux'], sky_level[pixeltables[1]]['flux'], sky_level[pixeltables[2]]['flux']])

    refsky = np.mean(sky_level_cube, axis=0)

    if plot:
        print('   ... creating the plot')

        plt.figure('median_slice', figsize=(20, 20))

        plotgrid = gridspec.GridSpec(5, 1)
        ax0 = plt.subplot(plotgrid[0, 0])
        ax1 = plt.subplot(plotgrid[1, 0])
        ax2 = plt.subplot(plotgrid[2, 0])
        ax3 = plt.subplot(plotgrid[3, 0])
        ax4 = plt.subplot(plotgrid[4, 0])
        spline_ax = [ax1, ax2, ax3]

    spline_dict = {}
    spline_funct_dict = {}

    for ax, skyname in zip(spline_ax, np.sort(list(sky_level.keys()))):
        wav = sky_level[skyname]['lambda']
        sky = sky_level[skyname]['flux']

        print('   ... calculating the spline for: ', skyname)

        skyspline, spline_funct = _skyfit(refsky, wav, sky, sigma=1, s=0.1, order=5, plot=True, plot_axis=ax)
        spline_dict[skyname] = skyspline
        spline_funct_dict[skyname] = spline_funct

        if plot:
            ax.plot(wav, sky / refsky, c='k', zorder=0, label='continuum')
            ax.plot(wav, skyspline, c='red', zorder=2, label='spline')
            ax.set_ylim(0.5, 1.5)
            ax.set_ylabel('rel. flux')

            ax.legend()

            ax0.plot(wav, sky / np.median(sky), alpha=0.5)
            ax4.plot(wav, sky / skyspline / np.median(sky / skyspline), alpha=0.5)

            ax0.set_ylim(0.5, 1.5)
            ax4.set_ylim(0.5, 1.5)
            ax4.set_xlabel('wavelength [A]')
            ax0.set_ylabel('norm. flux')
            ax4.set_ylabel('cor. norm. flux')
            ax0.set_title(skyname)

    if plot:
        print('   ... saving the plot:', os.path.join(pixtab_dir, 'intra_cube_flux_cor.png'))
        plt.savefig(os.path.join(pixtab_dir, 'intra_cube_flux_cor.png'), dpi=300)
        plt.close()

    with Client(n_workers=n_CPU, threads_per_worker=1, memory_limit='24GB') as client:

        pixtable_future = client.scatter(pixtable_dict)

        tasks = [dask.delayed(_apply_spline)(pixtable_future[key_to_process], spline_funct_dict, pixtab_dir, key_to_process) for key_to_process in list(sky_level.keys())]
        futures = client.compute(tasks)
        progress(futures)
        wait(futures)
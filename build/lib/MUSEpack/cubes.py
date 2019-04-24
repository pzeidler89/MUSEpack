#!/usr/bin/env python

import sys
import os
import shutil
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from spectral_cube import SpectralCube
import montage_wrapper as montage


def wcs_corr(input_fits, input_prm, path=os.getcwd(), prm_path=os.getcwd(),
            output_file=None, out_frame=None, in_frame=None,
            wcsname='Pampelmuse', correctiontype='shift'):

    '''
    Args:
        inputfits : :obj:`str`
            The fully reduced datacube, whose WCS has to be corrected

        input_prm : :obj:`str`
            The prm file produced by Pampelmuse

    Kwargs:
        path : :obj:`str` (optional, Default: current directory)
            I/O path

        prm_path : :obj:`str` (optional, default: current directory)
            I/O path of prm file

        output_file : :obj:`str` (optional, default: update fitsheader of input
            fits)
            outputfile name

        output_frame : :obj:`str` (optional, default : input frame)
            coordinate frame of the output cube in case one want to
            change

        in_frame : :obj:`str` (optional, default: input frame)
            coordinate frame of the output cube in case it cannot be
            determined from the header information or it has to be manually
            changed

        wcsname : :obj:`str` (optional, default: ``Pampelmuse``)
            the name of the new, updated wcs

        correctiontype : :obj:`str` (optional, default: ``shift``)
            the type of distortion correction

            ``full``: the full 2D CD matrix,

            ``shift``: shift in XY only.

    '''


    cube = fits.open(path + '/' + input_fits + '.fits')
    pmr = fits.open(prm_path + '/' + input_prm + '.prm.fits')

    print('processing observation: ' + path + '/' + input_fits + '.fits')
    print('using prm file: ' + prm_path + '/' + input_prm + '.prm.fits')

    assert (len(cube) != 3 or len(cube) != 1),\
    'fits file has currently unsupported extensions: Please check'

    if len(cube) == 3:
        prihdr = cube[0].header
        sechdr = cube[1].header
        assert prihdr['INSTRUME'] != 'MUSE    ',\
        ' This is not a MUSE cube. Please check'

        if not in_frame:
            in_frame = prihdr['RADECSYS'].lower()
        print('MUSE cube detected')

    if len(cube) == 1:
        print('No MUSE cube all info in one extension')
        prihdr = cube[0].header
        sechdr = cube[0].header
        if 'RADECSYS' in prihdr:
            in_frame = prihdr['RADECSYS'].lower()

    assert in_frame is not None, 'No WCS frame provided'

    print(' Input WCS frame: ', in_frame)
    if out_frame is None:
        print('Output WCS frame: ', in_frame)
    else:
        print('Output WCS frame: ', out_frame)
    print('')

    A = np.nanmedian(pmr[4].data[0][1])
    B = np.nanmedian(pmr[4].data[1][1])
    C = np.nanmedian(pmr[4].data[2][1])
    D = np.nanmedian(pmr[4].data[3][1])
    x0 = np.nanmedian(pmr[4].data[4][1])
    y0 = np.nanmedian(pmr[4].data[5][1])

    CD = np.array([[A, C], [B, D]])
    r = np.array([[x0, 0.], [0., y0]])

    w = WCS(sechdr)
    w.wcs.name = 'MUSE'
    orig_hdr = w.to_header(relax=False, key='A')
    orig_hdr.remove('RADESYSA')

    #### coord sys change
    if out_frame != None:
        ref_ra = sechdr['CRVAL1']
        ref_dec = sechdr['CRVAL2']
        ref_coord = SkyCoord(ra=ref_ra * u.degree,\
        dec=ref_dec * u.degree, frame=in_frame)
        trans_ref_coord = ref_coord.transform_to(out_frame)

        sechdr['CRVAL1'] = trans_ref_coord.ra.value
        sechdr['CRVAL2'] = trans_ref_coord.dec.value
        sechdr.set('RADESYS', out_frame.upper())

    ref_xy = np.array([sechdr['CRPIX1'], sechdr['CRPIX2']])
    ref_xy_new = (np.dot(r, np.ones(ref_xy.T.shape))\
    + np.dot(CD, ref_xy.T)).swapaxes(-1, -0)

    ref_shift_x = ref_xy_new[0] + 1.
    ref_shift_y = ref_xy_new[1] + 1.

    sechdr['CRPIX1'] = ref_shift_x
    sechdr['CRPIX2'] = ref_shift_y

    if correctiontype == 'full':
        sechdr['CD1_1'] = (-1) * A * 0.2 / 3600.
        sechdr['CD1_2'] = C * 0.2 / 3600.
        sechdr['CD2_1'] = (-1) * B * 0.2 / 3600.
        sechdr['CD2_2'] = D * 0.2 / 3600.

    sechdr.set('WCSNAME', wcsname)

    sechdr.extend(orig_hdr)

    if output_file == None:
        cube.writeto(path + '/' + input_fits + '.fits', overwrite=True)
    else:
        cube.writeto(path + '/' + output_file + '.fits', overwrite=True)


def pampelmuse_cat(ra, dec, mag, filter, idx=None, path=os.getcwd(),
                   sat=0., mag_sat=None, ifs_sat=None, mag_limit=None):

    '''

    This modules uses input parameters to create a catalog that is compatible
    with pampelmuse

    Args:

        ra : :obj:`float`
            RA coordinates of the stars

        dec : :obj:`float`
            Dec coordinates of the stars

        mag : :obj:`float`
            magnitudes of the stars

        filter : :obj:`float`
            filter used for the magnitudes

    Kwargs:

    idx : :obj:`float` (optional, default : counting up from 1)
        catalog index of the stars

    sat : :obj`float` (optional, default: 0)
        value assigned to saturated sources in the catalog

    mag_sat : :obj:`float` (optional, default: :obj:`None`)
        magnitude that replaces saturated sources in the catalog

    path : :obj:`str` (optional, default: current directory)
        path of the output file
    '''

    if idx == None:
        id = np.arange(len(ra)) + 1
    else:
        id = idx

    if mag_limit != None:
        mag_lim_id = np.where(mag <= mag_limit)
        id = id[mag_lim_id]
        ra = ra[mag_lim_id]
        dec = dec[mag_lim_id]
        mag = mag[mag_lim_id]

    if mag_sat != None:
        sat_source = np.where(mag == sat)

    if ifs_sat != None:
        ifs_sat_id = id[np.where((mag < ifs_sat) & (mag != sat))]
        f_ifs_sat_id = open(path + '/ifs_sat_id.list', 'w')
        f_ifs_sat_id.write('[')
        for i in ifs_sat_id:
            f_ifs_sat_id.write(str(i) + ',')

        f_ifs_sat_id.write(']')
        f_ifs_sat_id.close()

    tab = Table([id, ra, dec, mag], names=('id', 'ra', 'dec',\
    str(filter).lower()), dtype=('i4', 'f8', 'f8', 'f8'))

    if mag_sat != None:
        tab[str(filter).lower()][sat_source] = mag_sat

    tab.write(path + '/' + str(filter).upper(),\
    format='ascii.basic', delimiter=',', overwrite=True)


def linemaps(input_fits, path=os.getcwd(), elements=None, wavelengths=None):
    '''
    This module is intended to create linemaps of specified lines/elements

    Args:
        input_fits : :obj:`str`
            The fully reduced datacube Pampelmuse has been run on

    Kwargs:
        path : :obj:`str` (optional, default: current directory)
            I/O path

        element : obj:`list` (optional)
            list of elements the linemaps shall be produced

        wavelength : :obj:`list` (optional)
            list of wavelength for givene elements, optional,
            must be given if `elements` is given

    '''

    #predefined elements and their wavelength
    if elements == None:
        wavelengths = [6562.80, 6716.47, 5006.84, 6583.41]
        elements = ['Ha', 'SII_6716', 'OIII_5007', 'NII_6583']

    if elements != None and wavelengths == None:
        print('Error: element but no wavelength given')
        sys.exit()

    if os.path.exists(path + 'temp/') == False:
        os.mkdir(path + 'temp/')

    cube = SpectralCube.read(path + '/' + input_fits, hdu=1, fomrat='fits')
    for wavelength, element in zip(wavelengths, elements):
        slab = cube.spectral_slab((wavelength - 3) * u.AA,\
        (wavelength + 3) * u.AA).sum(axis=0)
        slab.hdu.writeto(path + '/' + element + '.fits', overwrite=True)

    shutil.rmtree(path + 'temp/')


def mosaics(input_list, name, path=os.getcwd()):

    '''
    This module is intended to create mosaics of specified lines/elements.
    linemaps should have been created beforehand using the ``linemaps`` module

    Args:
    input_list : :obj:`list`
            The list of specific linemaps to be used to mosaic

        name : :obj:`str`
            Name of the created mosaic

    Kwargs:
        path: :obj:`str` (optional, default: current directory)
            I/O path

    '''

    if os.path.exists(path + '/temp/') == False:
        os.mkdir(path + '/temp/')
    if os.path.exists(path + '/mosaics/') == False:
        os.mkdir(path + '/mosaics/')

    for idx, f in enumerate(input_list):
        shutil.copy(f, path + '/temp/' + str(idx) + '.fits')

    montage.mosaic(path + '/temp/', path + '/mosaic_temp/',\
    background_match=True, exact_size=True, cleanup=True)
    shutil.copy(path + '/mosaic_temp/mosaic_area.fits',\
    path + '/mosaics/exp_' + name + '.fits')

    shutil.copy(path + '/mosaic_temp/mosaic.fits',\
    path + '/mosaics/' + name + '.fits')
    shutil.rmtree(path + '/mosaic_temp/')
    shutil.rmtree(path + '/temp/')

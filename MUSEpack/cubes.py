#!/usr/bin/env python

__version__ = '0.3.0'

__revision__ = '20260812'

import sys
import os
import glob
import shutil
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clip
from spectral_cube import SpectralCube
import montage_wrapper as montage
from mpdaf.drs import PixTable
from pathlib import Path
from reproject import reproject_interp

''' internal modules'''
from MUSEpack.utils import ABtoVega
from MUSEpack.utils import _any_bit_in_number
from MUSEpack.utils import _make_circle_mask
from MUSEpack.utils import _make_rectangle_mask
from MUSEpack.utils import _make_diagnostic_plot

def wcs_cor(input_fits, offset_input, path=None, offset_path=None,
            output_file=None, output_path=None, out_frame=None, in_frame=None,
            correct_flux=False, spec_folder='stars', spec_path=None,
            hst_filter=['ACS', 'F814W'], AB_zpt=None, Vega_zpt=None,
            correctiontype='shift', save_output=True, debug=False, spec_snr=5.0,
            dqs = [2,4,8], man_multiplier=None):

    '''
    Args:
        inputfits : :obj:`str`
            The fully reduced datacube, whose WCS has to be corrected

        offset_input : :obj:`str`
            The prm file produced by Pampelmuse
            or the OFFSET_LIST.fits file from the ``exp_align`` routine 

    Kwargs:
        path : :obj:`str` (optional, Default: current directory)
            I/O path

        offset_path : :obj:`str` (optional, default: current directory)
            I/O path of prm file

        output_file : :obj:`str` (optional, default: input file name +_cor)
            outputfile name

        output_path : :obj:`str` (optional, default: current directory)
            output path. If not provided the I/O folder is identical

        output_frame : :obj:`str` (optional, default: input frame)
            coordinate frame of the output cube in case one want to
            change. Is always set to default input frame if corrections
            are not based on a .prm file 

        in_frame : :obj:`str` (optional, default: input frame)
            coordinate frame of the output cube in case it cannot be
            determined from the header information or it has to be manually
            changed

        correct_flux : :obj:`bool` (optional, default: :obj:`False`)
            If set :obj:`True` the fluxes of the data cube will be corrected
            to match the input catalog to correct for calibration offsets.
            If the input file is a prm-file:
            This step is only recommended if the input fluxes can be trusted.
            CUBEFIT and GETSPECTRA have to executed again using the corrected
            data cube to correct the prm file and to extract the corrected
            spectra.
            If the input file is a ESO OFFSET-LIST.fits file:
            The fluxes will be scaled based on the 'FLUX_SCALE' entries.

        save_output : :obj:`bool` (optional, default: :obj:`True`)
            If set :obj:`False` the corrected cube will not be saved. This can
            be used to just obtain the correction factors

        spec_folder : :obj:`str` (optional, default: ``stars``)
            The folder name, in which the the extracted stellar spectra are
            stored. This is only needed if correct_flux=:obj:`True` and the
            corrections are based on a .prm file

        spec_path : :obj:`str` (optional, default: current directory)
            I/O path of the ``spec_folder``. This keyword is only need if
            the corrections are based on a .prm file

        hst_filter : :obj:`list' (optional, default: ['ACS', 'F814W'])
            The HST filter used when convolving the spectra to measure the
            fluxes. It is used with STSYNPHOT to obtain the correct zeropoint
            conversion.

        AB_zpt : :obj:`float` (optional, default: None)
            AB zeropoint for any filter that is not part of SYNPHOT. If provided
            the Vega_zpt must be provided as well. hst_filter will be overwritten

        Vega_zpt : :obj:`float` (optional, default: None)
            Vega zeropoint for any filter that is not part of SYNPHOT. If provided
            the AB_zpt must be provided as well. hst_filter will be overwritten

        correctiontype : :obj:`str` (optional, default: ``shift``)
            the type of distortion correction

            ``full``: the full 2D CD matrix, only possible if based on a 
            .prm file

            ``shift``: shift in XY only.

        debug : :obj:`bool` (optional, default: :obj:`False`)
            If set :obj:`True` tadditional information will be printed

        spec_snr : :obj:`float` (optional, default: 5.0)
            the S/R limit when a spectra is used to calculate the flux correction.
            This limit might need to be lowered for low S/R obsverations

        dqs : :obj:`list` (optional, default: [2,4,8])
            list of DQ flags where source as dicarted.

        man_multiplier : :obj:`float` (optional, default: None)
            If the manual multiplier is set this factor will be applied for the flux calibration
            The automatic input will be omitted

    '''

    if path == None:
        path = os.getcwd()

    if offset_path == None:
        offset_path = os.getcwd()

    if spec_path == None:
        spec_path = os.getcwd()

    cube = fits.open(path + '/' + input_fits + '.fits')
    offset = fits.open(offset_path + '/' + offset_input)

    print('processing observation: ' + path + '/' + input_fits + '.fits')

    if len(offset) == 6:
        print('using prm file: ' + offset_path + '/' + offset_input)
        offset_type = 'prm'
    if len(offset) == 2:
        print('using OFFSET file: ' + offset_path + '/' + offset_input)
        offset_type = 'eso'
        spec_folder = None
        spec_path = None,
        correctiontype = 'shift'
        output_frame = None
        in_frame = None

    assert len(offset) != 6 or len(offset) != 2,\
    'This offset file is currently not supported: Please check'

    assert (len(cube) != 3 or len(cube) != 1),\
    'fits file has currently unsupported extensions: Please check'

    if offset_type == 'prm':

        if len(cube) == 3:
            prihdr = cube[0].header
            sechdr = cube[1].header
            assert prihdr['INSTRUME'].strip() == 'MUSE',\
            ' This is not a MUSE cube. Please check'

            if not in_frame:
                if 'RADESYS' in prihdr:
                    in_frame = prihdr['RADESYS'].lower()
                elif 'RADECSYS' in prihdr:
                    in_frame = prihdr['RADECSYS'].lower()
            print('MUSE cube detected')

        if len(cube) != 3:
            print('No MUSE cube header info assumed to be in one extension')
            prihdr = cube[0].header
            sechdr = cube[0].header
            if 'RADESYS' in prihdr:
                in_frame = prihdr['RADESYS'].lower()
            elif 'RADECSYS' in prihdr:
                in_frame = prihdr['RADECSYS'].lower()
            if correct_flux:
                print('Flux correction currently only supported'\
                + ' for MUSE cubes')

        assert in_frame is not None, 'No WCS frame provided'

        print(' Input WCS frame: ', in_frame)
        if out_frame is None:
            print('Output WCS frame: ', in_frame)
        else:
            print('Output WCS frame: ', out_frame)
        print('')

        A = np.nanmedian(offset[4].data[0][1])
        B = np.nanmedian(offset[4].data[1][1])
        C = np.nanmedian(offset[4].data[2][1])
        D = np.nanmedian(offset[4].data[3][1])
        x0 = np.nanmedian(offset[4].data[4][1])
        y0 = np.nanmedian(offset[4].data[5][1])

        CD = np.array([[A, C], [B, D]])
        r = np.array([[x0, 0.], [0., y0]])

        #### coord sys change
        if out_frame != None:
            ref_ra = sechdr['CRVAL1']
            ref_dec = sechdr['CRVAL2']
            ref_coord = SkyCoord(ra=ref_ra * u.degree,\
            dec=ref_dec * u.degree, frame=in_frame)
            trans_ref_coord = ref_coord.transform_to(out_frame)

            sechdr['CRVAL1'] = trans_ref_coord.ra.value
            sechdr['CRVAL2'] = trans_ref_coord.dec.value
            if 'RADESYS' in sechdr:
                sechdr.set('RADESYS', out_frame.upper())
            elif 'RADECSYS' in sechdr:
                sechdr.set('RADECSYS', out_frame.upper())

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

        #### flux correction
        if correct_flux and len(cube) == 3:

            print('using spec path: ' + spec_path)

            if AB_zpt != None and Vega_zpt != None:
                print('AB zeropoint and Vega zeropoint provided. hst_filter will'
                      ' be overwritten')
                aboffset = Vega_zpt - AB_zpt
            elif AB_zpt != None and Vega_zpt == None:
                print('both AB zeropoint and Vega zeropoint must be provided.')
                sys.exit()
            elif AB_zpt == None and Vega_zpt != None:
                print('both AB zeropoint and Vega zeropoint must be provided.')
                sys.exit()
            else:
                print('HST zeropoints used with STSYNPHOT')
                aboffset = ABtoVega(hst_filter[0], hst_filter[1])

            speclist = glob.glob(spec_path + '/' + spec_folder + '/specid*')

            flux_cor_tab = Table()
            cat_mag_list = []
            del_mag_list = []
            qltflag_list = []
            snrspec_list = []
            starid_list = []

            for i, temp_sp in enumerate(speclist):

                spec_hdu = fits.open(temp_sp)
                spec_head = spec_hdu[0].header

                qltflag = spec_head['HIERARCH SPECTRUM QLTFLAG']
                snrspec = spec_head['HIERARCH SPECTRUM FITSNR']

                if (_any_bit_in_number(dqs, qltflag) & (snrspec > spec_snr)):

                    del_mag_list = np.append(del_mag_list,\
                    -(spec_head['HIERARCH SPECTRUM MAG DELTA'] + 50. + aboffset))
                    cat_mag_list = np.append(cat_mag_list, spec_head['HIERARCH STAR MAG'])
                    qltflag_list = np.append(qltflag_list, qltflag)
                    snrspec_list = np.append(snrspec_list, snrspec)
                    starid_list = np.append(starid_list, spec_head['HIERARCH STAR ID'])

                muse_mag_list = np.array(cat_mag_list) - np.array(del_mag_list)
            if len (del_mag_list) == 0:
                print()
                sys.exit('No sources were slected. Please adjust SNR or make sure any sources are in the Cube'.upper())

            flux_cor_tab = Table([starid_list, muse_mag_list, cat_mag_list, del_mag_list, snrspec_list, qltflag_list],
                                 names=('spec_id', 'muse_mag', 'cat_mag', 'dmag', 'spec_snr', 'DQ'))
            flux_cor_tab['spec_id'].info.format = '5.0f'
            flux_cor_tab['muse_mag'].info.format = '5.2f'
            flux_cor_tab['cat_mag'].info.format = '5.2f'
            flux_cor_tab['dmag'].info.format = '5.2f'
            flux_cor_tab['spec_snr'].info.format = '5.1f'
            flux_cor_tab['DQ'].info.format = '3.0f'

            clippend_del_mag = sigma_clip(del_mag_list, sigma=2, cenfunc = np.ma.median)
            if man_multiplier == None:
                fmultipl = 10 ** ((-1) * 0.4 * np.ma.median(clippend_del_mag))
            else:
                fmultipl = man_multiplier

            flux_cor_tab.add_column(np.ma.getmask(clippend_del_mag), name='masked')
            flux_cor_tab.sort('spec_id')
            flux_cor_tab.write(os.path.join(output_path, 'flux_cor_info.tab'), format='ascii.rst', overwrite=True)

            if debug == True:
                print(flux_cor_tab)
                print()

            print('The magnitude difference:')
            print('cat - MUSE [mag]: ',\
            '{:.2f}'.format(np.ma.median(clippend_del_mag)))
            print('sigma cat - MUSE [mag]: ',\
            '{:.2f}'.format(np.ma.std(clippend_del_mag)))
            if man_multiplier != None:
                print("WARNING: using manual multiplyer")
            print('The flux multiplicator [f_catalog / f_MUSE]: ',\
            '{:.2f}'.format(fmultipl))
            print('number of sources: ',\
            '{:.0f}'.format(clippend_del_mag.count()))

            with open(os.path.join(output_path, 'flux_cor_info.tab'),"a+") as f:
                f.write("\n")
                f.write('cat - MUSE [mag]: ' + '{:.2f}'.format(np.ma.median(clippend_del_mag)) + '\n')
                f.write('sigma cat - MUSE [mag]: ' + '{:.2f}'.format(np.ma.std(clippend_del_mag)) + '\n')
                if man_multiplier != None:
                    f.write(("WARNING: using manual multiplyer \n"))
                f.write('The flux multiplicator [f_catalog / f_MUSE]: ' + '{:.2f}'.format(fmultipl) + '\n')
                f.write('N sources: ' + '{:.0f}'.format(clippend_del_mag.count()))

            cube['DATA'].data *= fmultipl
            cube['STAT'].data *= fmultipl ** 2

            if save_output == True:
                if output_file == None:
                    offset[0].header['HIERARCH PAMPELMUSE global prefix'] = offset_input[:-9] + '_cor'
                    offset.writeto(offset_path + '/' + offset_input[:-9] + '_cor.prm.fits', overwrite=True)

                else:
                    shutil.copyfile(offset_path + '/' + offset_input,\
                    offset_path + '/' + output_file)

    if offset_type == 'eso':

        obs = cube[0].header['MJD-OBS']
        indobs = np.where(obs == offset[1].data['MJD_OBS'])
        dRA = offset[1].data['RA_OFFSET'][indobs]
        dDEC = offset[1].data['DEC_OFFSET'][indobs]
        fscale = offset[1].data['FLUX_SCALE'][indobs]

        print('The RA offset is [arcsec]: ',\
        '{:.4f}'.format(dRA[0] * 3600.))
        print('The Dec offset is [arcsec]: ',\
        '{:.4f}'.format(dDEC[0] * 3600.))
        print('The flux scale is: ',\
        '{:.4f}'.format(fscale[0]))

        cube[1].header['CRVAL1'] -= dRA[0]
        cube[1].header['CRVAL2'] -= dDEC[0]

        if correct_flux and not np.isnan(fscale[0]):
            cube[1].data *= fscale[0]

    if save_output==True:
        if output_path == None:
            if output_file == None:
                cube.writeto(path + '/' + input_fits + '_cor.fits', overwrite=True)
            else:
                cube.writeto(path + '/' + output_file + '.fits', overwrite=True)
        else:
            if output_file == None:
                cube.writeto(os.path.join(output_path, input_fits + '_cor.fits'), overwrite=True)
            else:
                cube.writeto(os.path.join(output_path, output_file ), overwrite=True)


def pampelmuse_cat(ra, dec, mag, filter, idx=None, path=None,
                   sat=0., mag_sat=None, ifs_sat=None, mag_limit=None,
                   regsize=0.5):

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

        filter : :obj:`str`
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

    mag_limit : :obj:`float` (optional, default: :obj:`None`)
        the magnitude at which the output catalog should be truncated

    regsize : :obj:`float` (optional, default: :float:0.5)
        the size of the regions in arcsec
    '''

    if not path:
        path = os.getcwd()

    if not idx.all():
        id = np.arange(len(ra)) + 1
    else:
        id = idx

    if mag_limit:
        mag_lim_id = np.where(mag <= mag_limit)
        id = id[mag_lim_id]
        ra = ra[mag_lim_id]
        dec = dec[mag_lim_id]
        mag = mag[mag_lim_id]

    if mag_sat:
        sat_source = np.where(mag == sat)

    if ifs_sat:
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


    #write ds9 regions file
    regf = open(path + '/' + str(filter).upper() + '.reg', 'w')

    regf.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=0 delete=1 include=1 source=1\n')
    regf.write("fk5\n")

    for rai, deci in zip(tab['ra'], tab['dec']):
        regf.write("circle(" + str(rai) + ", " + str(deci) + ', '
                   + '{:.1f}'.format(regsize) +  '")\n')

    regf.close()


def linemaps(input_fits, path=None, elements=None, wavelengths=None, dlambda=3):
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

        dlambda : :obj:`float` (optional, default: 3)
            plus/minus lambda range in Angstrom around line centroid to construct
            the linemap

    '''

    if path == None:
        path = os.getcwd()

    #predefined elements and their wavelength
    if elements == None:
        wavelengths = [6562.80, 6716.47, 5006.84, 6583.41]
        elements = ['Ha', 'SII_6716', 'OIII_5007', 'NII_6583']

    if elements != None and wavelengths.any() == None:
        print('Error: element but no wavelength given')
        sys.exit()

    if os.path.exists(path + 'temp/') == False:
        os.mkdir(path + 'temp/')

    cube = SpectralCube.read(path + '/' + input_fits, hdu=1, format='fits')
    cube.allow_huge_operations = True
    for wavelength, element in zip(wavelengths, elements):
        slab = cube.spectral_slab((wavelength - dlambda) * u.AA,\
        (wavelength + dlambda) * u.AA).sum(axis=0)
        slab.hdu.writeto(path + '/' + element + '_' +input_fits, overwrite=True)

    shutil.rmtree(path + 'temp/')


def mosaics(input_list, name, path=None, background_match=True):

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

        background_match: :obj:`bool` (optional, default: `true`)
            turn on and off teh background matching


    '''

    if path == None:
        path = os.getcwd()

    if os.path.exists(path + '/temp/') == False:
        os.mkdir(path + '/temp/')
    if os.path.exists(path + '/mosaics/') == False:
        os.mkdir(path + '/mosaics/')

    for idx, f in enumerate(input_list):
        shutil.copy(f, path + '/temp/' + str(idx) + '.fits')

    montage.mosaic(path + '/temp/', path + '/mosaic_temp/',\
    background_match=background_match, exact_size=True, cleanup=True)
    shutil.copy(path + '/mosaic_temp/mosaic_area.fits',\
    path + '/mosaics/exp_' + name + '.fits')

    shutil.copy(path + '/mosaic_temp/mosaic.fits',\
    path + '/mosaics/' + name + '.fits')
    shutil.rmtree(path + '/mosaic_temp/')
    shutil.rmtree(path + '/temp/')


def apply_leak_mask(leak_mask, pixtable, flag_bad_pix = 1):

    '''
    This module is intended to apply external leak masks to the PIXTABLES of each
    individual dither. This allows the user to mask any spilling generated by bright
    pointsources that heavily saturated.

    Args:
        leak_mask : :obj:`str`
                The leak mask that should be applied. The format needs to be identical to
                 the sky mask

        pixtable : :obj:`str`
                The PIXTABLE the leak mask should be applied to.

    Kwargs:
        flag_bad_pix : :obj:`int` (optional, default: 1)
                This is the flag that will be set for a masked pixel. The default
                sets it to BAD_PIXEL

    '''

    print('   ... loading the mask and pixel table')
    pixtable_filename = Path(pixtable).name
    pixtable_path = Path(pixtable).parent

    with fits.open(leak_mask) as mask_hdu:
        mask_data = mask_hdu[1].data.astype(bool)
        mask_wcs = WCS(mask_hdu[1].header)

    pix = PixTable(pixtable)

    print('   ... creating the coordinate transforms')

    ra, dec = pix.get_pos_sky()

    x_pix, y_pix = mask_wcs.world_to_pixel_values(ra, dec)
    x_idx = np.round(x_pix).astype(int)
    y_idx = np.round(y_pix).astype(int)

    ### the next step is needed in case rounding puts the pixel outside the valid range

    valid_bounds = ((x_idx >= 0) & (x_idx < mask_data.shape[1]) &
                    (y_idx >= 0) & (y_idx < mask_data.shape[0]))

    print('   ... masking the PIXTABLE')

    pix_sel = np.ones(len(ra), dtype=bool)
    pix_sel[valid_bounds] = ~mask_data[y_idx[valid_bounds], x_idx[valid_bounds]]

    dq = pix.get_dq()
    dq[pix_sel] |= flag_bad_pix
    pix.set_dq(dq)

    pix_tab_out = os.path.join(pixtable_path, pixtable_filename.replace('.fits', '_MASKED.fits'))
    print('   ... saving the PIXTABLE: ', pix_tab_out)
    pix.write(pix_tab_out)


def build_leak_mask(whitelight_image, ra_cent, dec_cent, bleed_width, bleend_length, r_source=None, mask_source=False, offset_list=None):

    '''
    This module is intended to create the external leak mask.

    Args:
        whitelight_image : :obj:`str`
                The whitelight image with the same WCS as teh data cube

        bleed_width : :obj:`np.array`
                numpy array with bleed width in arcsec of one or multiple masks

        bleed_length : :obj:`np.array`
                numpy array with bleed length in arcsec of one or multiple masks

    Kwargs:
        r_source : :obj:`np.array` (optional, default: `None')
                The source radius to be masked. If `None' the source radius matches the bleed_mask

        mask_source : :obj:`bool` (option, default: `False')
                If true the source is also masked.

        offset_list : :obj:`str` (optional, default: `None')
                The path to the OFFSET_LIST.fits file created in the cube align step of teh MUSE pipeline.
                If provided it will correct the mask location for possiple shifts between different dithers
    '''

    print('   >>> load whitelight image')
    whitelight_path = Path(whitelight_image).parent
    img_hdu = fits.open(whitelight_image)
    img = img_hdu[1].data
    img_header = img_hdu[1].header


    wcs2d = WCS(img_header)
    bleed_pa = img_hdu[0].header['HIERARCH ESO INS DROT POSANG']
    date_obs = img_hdu[0].header['DATE-OBS']

    shape = img.data.shape
    print(f"Image shape: {shape}")

    if offset_list != None:
        offset_list_hdu = fits.open(offset_list)
        offset_tab = Table(offset_list_hdu[1].data)

        offset_select = np.where(date_obs == offset_tab['DATE_OBS'])
        offset_ra = offset_tab['RA_OFFSET'][offset_select].value[0]
        offset_dec = offset_tab['DEC_OFFSET'][offset_select].value[0]

    else:
        offset_ra = 0.
        offset_dec = 0.

    print('   >>> initiate the mask')
    mask = np.ones(shape, dtype=np.int32)

    src_geom_dict = {}
    src_mask_dict = {}
    bleed_geom_dict = {}
    bleed_mask_dict = {}

    for imask, (ra, dec, width, length) in enumerate(zip(ra_cent, dec_cent, bleed_width, bleend_length)):

        print(ra)
        ra += offset_ra * np.cos(dec/180.*np.pi)
        dec += offset_dec
        print(ra)

        if r_source==None:
            print('   ... no radius is given, source will not be masked')
            r= width / np.sqrt(2)
        else:
            if isinstance(r_source, float):
                r = r_source
            else:
                r = r_source[imask]

        src_mask, src_geom = _make_circle_mask(shape, wcs2d, ra, dec, r)
        bleed_mask, bleed_geom = _make_rectangle_mask(shape, wcs2d, ra, dec, bleed_pa, length, width)

        src_geom_dict['m'+str(imask)] = src_geom
        src_mask_dict['m'+str(imask)] = src_mask
        bleed_geom_dict['m'+str(imask)] = bleed_geom
        bleed_mask_dict['m'+str(imask)] = bleed_mask

        print('   >>> apply mask ',imask+1, 'of', len(ra_cent))
        mask[bleed_mask] = 0

    for mask_dict_key in src_mask_dict.keys():
        if mask_source:
            mask[src_mask_dict[mask_dict_key]] = 0
        else:
            mask[src_mask_dict[mask_dict_key]] = 1

    n_total = mask.size
    n_sky   = int(mask.sum())
    n_excl  = n_total - n_sky
    print(f"\nMask summary:")
    print(f"  total spaxels: {n_total}")
    print(f"  sky spaxels:   {n_sky}  ({100*n_sky/n_total:.1f}%)")
    print(f"  excluded:      {n_excl} ({100*n_excl/n_total:.1f}%)")
    if n_sky / n_total < 0.5:
        print("  WARNING: less than 50% of the FOV is sky.")
        print("  Auto-calibration may be noisy. Consider a tighter mask.")

    print('   >>> create LEAK_MASK.fits')

    img_header["EXTNAME"] = "DATA"
    img_header["BUNIT"]   = ""

    primary = fits.PrimaryHDU()
    image_hdu = fits.ImageHDU(data=mask, header=img_header, name="DATA")
    hdul = fits.HDUList([primary, image_hdu])
    leak_mask_file = whitelight_image.replace('IMAGE_FOV', 'LEAK_MASK')
    hdul.writeto(leak_mask_file, overwrite=True)
    print(f"\nWrote {leak_mask_file}")

    _make_diagnostic_plot(img, mask, src_geom_dict, bleed_geom_dict, whitelight_image.replace('IMAGE_FOV', 'LEAK_MASK').replace('fits','png'))


def unmask_overlapping_leak_mask(prime_mask_file, sec_mask_file ,mask_val=0, unmask_val=1):

    '''
    This module is intended to unmask parts of a leak mask to avoid multiple masked bleeds removing all data.
    sec_mask is the mask that will be unmasked.

    Args:
        prime_mask_file : :obj:`str`
                The LEAK_MASK which will be the dominant mask

        prime_mask_file : :obj:`str`
                The LEAK_MASK where the overlapping regions will be unmasked

    Kwargs:
        mask_val : :obj:`int` (optional, default: 0)
                The value which is used to mask data

        unmask_val : :obj:`int` (optional, default: 0)
                The value which is unused to mask data
    '''

    prime_mask_hdu = fits.open(prime_mask_file)
    sec_mask_hdu = fits.open(sec_mask_file)

    prime_mask = prime_mask_hdu[1].data
    prime_mask_header = prime_mask_hdu[1].header
    prime_mask_wcs = WCS(prime_mask_header)

    sec_mask = sec_mask_hdu[1].data
    sec_mask_header = sec_mask_hdu[1].header
    sec_mask_wcs = WCS(sec_mask_header)

    prime_mask_reprojected, footprint = reproject_interp(
        (prime_mask.astype(float), prime_mask_wcs), sec_mask_wcs,
        shape_out=sec_mask.shape,
        order='nearest-neighbor'
    )

    prime_mask_on_sec_mask = (prime_mask_reprojected == mask_val) & ~np.isnan(prime_mask_reprojected)
    sec_mask_bool = (sec_mask == mask_val)

    # 3. Identify where BOTH masks are masked
    both_masked = prime_mask_on_sec_mask & sec_mask_bool

    # 4. Unmask those pixels in mask2
    updated_sec_mask = sec_mask.copy()
    updated_sec_mask[both_masked] = unmask_val

    sec_mask_hdu[1].data = updated_sec_mask
    sec_mask_hdu.writeto(sec_mask_file, overwrite=True)

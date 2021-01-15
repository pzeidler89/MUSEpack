#!/usr/bin/env python

__version__ = '0.1'

__revision__ = '20210111'

import string
from astropy.io import fits
import logging
from astropy.coordinates import SkyCoord
from astropy import units as u


class bck_sub:

    '''
    Args:


    Kwargs:
        debug : :obj:`bool`, (optional), default: :obj:`None`
        :obj:`True`: :class:`MUSEreduce.bck_sub` runs in debug mode and
        all modules will be executed in single-thread mode.

        loglevel : :obj:`str` (optional, default: ``INFO``)
        ``DEBUG``: all functions run on a single core to obtain proper
        output in the correct order for proper debugging.

    '''


    def __init__(self, input_fits, loglevel="INFO", file_id=None):

        self.debug = False
        self.loglevel = loglevel
        self.file_id = file_id

        hdu_cube = fits.open(input_fits)
        header_data = hdu_cube[1].header

        if not file_id:
            RA = hdu_cube[0].header['RA']
            DEC = hdu_cube[0].header['DEC']
            EXPTIME = hdu_cube[0].header['EXPTIME']

            c = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree,\
            frame='fk5').to_string('hmsdms', sep='',\
            precision=0).translate(str.maketrans('', '', string.whitespace))

            self.file_id = c + '_' + str(int(EXPTIME)).rjust(4, '0')

        self.logger = logging.getLogger(self.file_id)

        if self.loglevel == "INFO":
            self.logger.setLevel(logging.INFO)
        if self.loglevel == "DEBUG":
            self.logger.setLevel(logging.DEBUG)

        # setting up the logger

        log = logging.FileHandler(str(self.file_id) + '.log', mode='w')

        log.setLevel(self.loglevel)
        formatter = logging.Formatter('%(asctime)s [%(levelname)8s]: ' \
                                      + '%(message)s [%(funcName)s]', datefmt="%Y-%m-%d %H:%M:%S")
        log.setFormatter(formatter)

        self.logger.addHandler(log)
        self.logger.info('bck_sub vers. ' + str(__version__) \
                         + ' rev. ' + str(__revision__))
        self.logger.info('Initiate instance of bck_sub for: ' \
                         + str(self.file_id))
        self.logger.debug('DEBUG mode => all modules run on one core')


        self.logger.info('Loading data')

        self.data = hdu_cube[1].data
        self.noise = hdu_cube[2].data


    def set_debug(self, dbg_mode):
        '''
        Args:
            dbg_mode : :obj:`bool`,
            to deactivate the debug mode set to :obj:`False`
            to activate the debug mode set to :obj:`True`

        '''

        self.debug = dbg_mode
        if self.loglevel == "INFO":
            self.logger.setLevel(logging.INFO)
        if self.loglevel == "DEBUG":
            self.logger.setLevel(logging.DEBUG)
#!/usr/bin/env python

'''
vers. 0.1.0: Executes the MUSE reduction pileline in the correct order
vers. 0.1.1: introducing only one calibration file folder per OB
vers. 0.1.2: choosing the illumination file closest to the observation
vers. 0.1.3: selecting the files for the different master file creations
vers. 0.1.4: minor corrections
vers. 0.1.5: looping order changed. each module loops by itself
vers. 0.1.6: always checking if calibration files already exist
vers. 0.1.7: Choosing between ESO calibrations and own calibrations
vers. 0.1.8: User can choose specific exposure time to reduce
vers. 0.2.0: exposures spread via different OBs for one pointing is supported.
             To do so, run the script normally for each OB including scipost.
             After all are processed then run exp_align and exp_combine 
             AO observations are supported
             reduction of multiple OBs in one run supported
             set rootpath manually
             Extend the toggles
vers. 0.2.1: user can choose the number of CPUs used
vers. 0.2.2: sky subtraction can be modified, so that individual elements can
             be excluded
vers. 0.2.3: use json file as input
vers. 0.2.4: additional parameters added: skyreject, skysubtraction
             set parameters and modules will be shown
vers. 0.2.5: bug fixes
vers. 0.3.0: using the correct ILLUM file for the STD reduction in the
             SCI_BASIC routine, SCI_BASIC now separated for STD and OBJECT
             reduction
vers. 0.3.1: one can select if the sof file are created automatically or
             provided by the user
vers. 0.4.0: supports now pipeline 2.4.2 and the NFM-AO
             added: pipeline_path
             choosing if darks may be used
             only reduces STD once per OB
             general use of external SKY fields
             collecting the files for exp_combine in an independent step
             exp_align is an independent step now
vers. 0.4.1: new file names to correct a problem where data gets replaced
             in the scipost routine if you reduce the data with and
             without sky
vers. 0.4.2: one can now change the ignore and fraction parameters in the
             JSON file
vers. 0.4.3: one can auto remove and rewrite the statics
vers. 0.4.4: changed the sky subtraction keyword
             the user can give now individual names for the different dither
             exposures: Does currently not work with multiple OBs or multiple
             pointings per OB

'''

__version__ = '0.4.4'

__revision__ = '20181223'

class MUSEreduce:
    def __init__(self):
        self
        

import sys,shutil,os,subprocess,glob,string,filecmp
import numpy as np
from astropy.io import ascii,fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import datetime,timedelta
import time
import json




def get_filelist(data_dir,rootpath,filename_wildcard):
    os.chdir(data_dir)
    raw_data_list=glob.glob(filename_wildcard)
    os.chdir(rootpath)
    return raw_data_list
    
def call_esorex(exec_dir,rootpath,esorex_cmd,n_CPU):
    os.chdir(exec_dir)
    os.system('export OMP_NUM_THREADS='+str(n_CPU))
    print('esorex '+ esorex_cmd)
    os.system('esorex '+ esorex_cmd)
    os.chdir(rootpath)

def sort_data(rootpath,raw_data_dir,working_dir,ESO_calibration_dir):
        
    file_list=get_filelist(raw_data_dir,rootpath,'*.fits*')
    science_files=np.array([])
    calibration_files=np.array([])
    science_type=np.array([])
    calibration_type=np.array([])
    ESO_calibration_files=np.array([])
    ESO_calibration_type=np.array([])
    
    for files in range(len(file_list)):
        
        hdu = fits.open(raw_data_dir+file_list[files])
        
        dprcatg_exist=hdu[0].header.get('HIERARCH ESO DPR CATG',False)
        procatg_exist=hdu[0].header.get('HIERARCH ESO PRO CATG',False)
        
        if dprcatg_exist:
            
            dprcatg=(hdu[0].header['HIERARCH ESO DPR CATG'])
            dprtype=(hdu[0].header['HIERARCH ESO DPR TYPE'])
                
            if dprcatg == 'SCIENCE':
                science_files = np.append(science_files,file_list[files])
                science_type = np.append(science_type,dprtype)
            
            if dprcatg == 'CALIB':
                calibration_files = np.append(calibration_files,file_list[files])
                if dprtype == 'FLAT,LAMP': calibration_type = np.append(calibration_type,'FLAT')
                elif dprtype == 'FLAT,SKY': calibration_type = np.append(calibration_type,'SKYFLAT')
                elif dprtype == 'WAVE': calibration_type = np.append(calibration_type,'ARC')
                elif dprtype == 'WAVE,MASK': calibration_type = np.append(calibration_type,'MASK')
                elif dprtype == 'FLAT,LAMP,ILLUM': calibration_type = np.append(calibration_type,'ILLUM')
                else: calibration_type = np.append(calibration_type,dprtype)
                
        if procatg_exist:
            
            procatg=(hdu[0].header['HIERARCH ESO PRO CATG'])
                        
            ESO_calibration_files = np.append(ESO_calibration_files,file_list[files])
            if procatg == 'MASTER_BIAS':
                ESO_calibration_type = np.append(ESO_calibration_type,'MASTER_BIAS')
                if os.path.isfile(ESO_calibration_dir+'MASTER_BIAS.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'MASTER_BIAS.fits')
            if procatg == 'MASTER_DARK':
                ESO_calibration_type = np.append(ESO_calibration_type,'MASTER_DARK')
                if os.path.isfile(ESO_calibration_dir+'MASTER_DARK.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'MASTER_DARK.fits')
            if procatg == 'MASTER_FLAT':
                ESO_calibration_type = np.append(ESO_calibration_type,'MASTER_FLAT')
                if os.path.isfile(ESO_calibration_dir+'MASTER_FLAT.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'MASTER_FLAT.fits')
            if procatg == 'TRACE_TABLE':
                ESO_calibration_type = np.append(ESO_calibration_type,'TRACE_TABLE')
                if os.path.isfile(ESO_calibration_dir+'TRACE_TABLE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'TRACE_TABLE.fits')
            if procatg == 'WAVECAL_TABLE':
                ESO_calibration_type = np.append(ESO_calibration_type,'WAVECAL_TABLE')
                if os.path.isfile(ESO_calibration_dir+'WAVECAL_TABLE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'WAVECAL_TABLE.fits')
            if procatg == 'LSF_PROFILE':
                ESO_calibration_type = np.append(ESO_calibration_type,'LSF_PROFILE')
                if os.path.isfile(ESO_calibration_dir+'LSF_PROFILE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'LSF_PROFILE.fits')
            if procatg == 'TWILIGHT_CUBE':
                ESO_calibration_type = np.append(ESO_calibration_type,'TWILIGHT_CUBE')
                if os.path.isfile(ESO_calibration_dir+'TWILIGHT_CUBE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'TWILIGHT_CUBE.fits')
            if procatg == 'FILTER_LIST':
                ESO_calibration_type = np.append(ESO_calibration_type,'FILTER_LIST')
                if os.path.isfile(ESO_calibration_dir+'FILTER_LIST.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'FILTER_LIST.fits')
            if procatg == 'EXTINCT_TABLE':
                ESO_calibration_type = np.append(ESO_calibration_type,'EXTINCT_TABLE')
                if os.path.isfile(ESO_calibration_dir+'EXTINCT_TABLE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'EXTINCT_TABLE.fits')
            if procatg == 'STD_FLUX_TABLE':
                ESO_calibration_type = np.append(ESO_calibration_type,'STD_FLUX_TABLE')
                if os.path.isfile(ESO_calibration_dir+'STD_FLUX_TABLE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'STD_FLUX_TABLE.fits')
            if procatg == 'SKY_LINES':
                ESO_calibration_type = np.append(ESO_calibration_type,'SKY_LINES')
                if os.path.isfile(ESO_calibration_dir+'SKY_LINES.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'SKY_LINES.fits')
            if procatg == 'GEOMETRY_TABLE':
                ESO_calibration_type = np.append(ESO_calibration_type,'GEOMETRY_TABLE')
                if os.path.isfile(ESO_calibration_dir+'GEOMETRY_TABLE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'GEOMETRY_TABLE.fits')
            if procatg == 'ASTROMETRY_WCS':
                ESO_calibration_type = np.append(ESO_calibration_type,'ASTROMETRY_WCS')
                if os.path.isfile(ESO_calibration_dir+'ASTROMETRY_WCS.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'ASTROMETRY_WCS.fits')
            if procatg == 'STD_RESPONSE':
                ESO_calibration_type = np.append(ESO_calibration_type,'STD_RESPONSE')
                if os.path.isfile(ESO_calibration_dir+'STD_RESPONSE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'STD_RESPONSE.fits')
            if procatg == 'STD_TELLURIC':
                ESO_calibration_type = np.append(ESO_calibration_type,'STD_TELLURIC')
                if os.path.isfile(ESO_calibration_dir+'STD_TELLURIC.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'STD_TELLURIC.fits')
            if procatg == 'ASTROMETRY_REFERENCE':
                ESO_calibration_type = np.append(ESO_calibration_type,'ASTROMETRY_REFERENCE')
                if os.path.isfile(ESO_calibration_dir+'ASTROMETRY_REFERENCE.fits') == False: shutil.copy(raw_data_dir+file_list[files],ESO_calibration_dir+'ASTROMETRY_REFERENCE.fits')
                
    rot_angles=np.zeros(len(science_files))
    points=np.zeros(len(science_files),dtype=object)
    rot_angles_ident=np.zeros(len(science_files))
    
    for files in range(len(science_files)):
        
        hdu = fits.open(raw_data_dir+science_files[files])[0]
        RA=hdu.header['RA']
        DEC=hdu.header['DEC']
        EXPTIME=hdu.header['EXPTIME']
        
        points[files] = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5').to_string('hmsdms',sep='',precision=0).translate(str.maketrans('', '', string.whitespace))+'_'+str(int(EXPTIME)).rjust(4, '0')
        rot_angles[files]=hdu.header['HIERARCH ESO INS DROT POSANG']
        
    for idx in range(len(rot_angles)):
        ident = 0
        for idx2 in range(len(rot_angles)):
            if rot_angles[idx] == rot_angles[idx2] and points[idx] == points[idx2]:
                rot_angles_ident[idx2] = ident
                ident +=1
                
    for files in range(len(science_files)):
        
        hdu = fits.open(raw_data_dir+science_files[files])[0]
        RA=hdu.header['RA']
        DEC=hdu.header['DEC']
        EXPTIME=hdu.header['EXPTIME']
        ROT=hdu.header['HIERARCH ESO INS DROT POSANG']
        DATE=hdu.header['MJD-OBS']
        
        c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5').to_string('hmsdms',sep='',precision=0).translate(str.maketrans('', '', string.whitespace))
        filelist_science=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_SCI.list'
        filelist_sky=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_SKY.list'
        filelist_dark=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_DAR.list'
        filelist_twilight=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_TWI.list'
        
        dark_date=np.array([])
        twilight_date=np.array([])
                
        for calibfiles in range(len(calibration_files)):
            if calibration_type[calibfiles] == 'DARK':
                hdu_temp=fits.open(raw_data_dir+calibration_files[calibfiles])[0]
                temp_date=hdu_temp.header['MJD-OBS']
                dark_date=np.append(dark_date,temp_date)
            if calibration_type[calibfiles] == 'SKYFLAT':
                hdu_temp=fits.open(raw_data_dir+calibration_files[calibfiles])[0]
                temp_date=hdu_temp.header['MJD-OBS']
                twilight_date=np.append(twilight_date,temp_date)
        
        dark_date=np.unique(dark_date)
        twilight_date=np.unique(twilight_date)
        
        if science_type[files] == 'OBJECT': f_science = open(working_dir+filelist_science, 'w')
        if science_type[files] == 'SKY': f_science = open(working_dir+filelist_sky, 'w')
        f_dark = open(working_dir+filelist_dark, 'w')
        f_twilight = open(working_dir+filelist_twilight, 'w')
        
        f_science.write(raw_data_dir+science_files[files]+'  '+science_type[files]+'\n')
            
        for calibfiles in range(len(calibration_files)):
            temp_date=fits.open(raw_data_dir+calibration_files[calibfiles])[0].header['MJD-OBS']
            if abs(temp_date - DATE) <= 1.: f_science.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            if calibration_type[calibfiles] == 'STD' and abs(temp_date - DATE) > 1.:
                f_science.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
                if calibration_type[calibfiles] == 'ILLUM' and abs(temp_date - DATE) > 1.: f_science.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
                print('WARNING: STD star and SCI observation more than 24 hours apart')
            if (abs(temp_date - dark_date) <= 1.).all(): f_dark.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            if (abs(temp_date - twilight_date) <= 1.).all(): f_twilight.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            
        f_science.close()
        f_dark.close()
        f_twilight.close()

### CALIBRATION PRE-PROCESSING ###

def bias(rootpath,working_dir,exposure_list,exposure_list_DARK,exposure_list_TWILIGHT,dark,calibration_dir,creating_sof,n_CPU=24):
    
    print('... Creating the MASTER BIAS')
    
    esorex_cmd = '--log-file=bias.log --log-level=debug muse_bias --nifu=-1 --merge bias.sof'

    if creating_sof:
        
        if os.path.exists(calibration_dir+'SCIENCE/bias.sof') == True: os.remove(calibration_dir+'SCIENCE/bias.sof')
        if dark == True and os.path.exists(calibration_dir+'DARK/bias.sof') == True: os.remove(calibration_dir+'DARK/bias.sof')
        if os.path.exists(calibration_dir+'TWILIGHT/bias.sof') == True: os.remove(calibration_dir+'TWILIGHT/bias.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')
    
            raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
            if dark: raw_data_list_DARK = ascii.read(exposure_list_DARK[exposure_ID],format='no_header')
            raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
            f_science = open(calibration_dir+'SCIENCE/bias_temp.sof', 'w')
            if dark: f_dark = open(calibration_dir+'DARK/bias_temp.sof', 'w')
            f_twilight = open(calibration_dir+'TWILIGHT/bias_temp.sof', 'w')

            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'BIAS':
                    f_science.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            f_science.close()

            if dark:
                for i in range(len(raw_data_list_DARK[1][:])):
                    if raw_data_list_DARK[i][1] == 'BIAS':
                        f_dark.write(raw_data_list_DARK[i][0]+'  '+raw_data_list_DARK[i][1]+'\n')
                f_dark.close()

            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'BIAS':
                    f_twilight.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            f_twilight.close()

            if os.path.isfile(calibration_dir+'SCIENCE/bias.sof'):
                if filecmp.cmp(calibration_dir+'SCIENCE/bias.sof',calibration_dir+'SCIENCE/bias_temp.sof'):
                    print('calibration files are identical: proceed without error')
                    os.remove(calibration_dir+'SCIENCE/bias_temp.sof')
                else:
                    sys.exit('CAUTION CALIBRATION FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'SCIENCE/bias_temp.sof',calibration_dir+'SCIENCE/bias.sof')
                call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
            
            if os.path.isfile(calibration_dir+'TWILIGHT/bias.sof'):
                if filecmp.cmp(calibration_dir+'TWILIGHT/bias.sof',calibration_dir+'TWILIGHT/bias_temp.sof'):
                    print('TWILIGHT files are identical: proceed without error')
                    os.remove(calibration_dir+'TWILIGHT/bias_temp.sof')
                else:
                    sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'TWILIGHT/bias_temp.sof',calibration_dir+'TWILIGHT/bias.sof')
                call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)
        
            if dark:
                if os.path.isfile(calibration_dir+'DARK/bias.sof'):
                    if filecmp.cmp(calibration_dir+'DARK/bias.sof',calibration_dir+'DARK/bias_temp.sof'):
                        print('dark files are identical: proceed without error')
                        os.remove(calibration_dir+'DARK/bias_temp.sof')
                    else:
                        sys.exit('CAUTION DARK FILES ARE DIFFERENT: PLEASE CHECK')
                else:
                    os.rename(calibration_dir+'DARK/bias_temp.sof',calibration_dir+'DARK/bias.sof')
                    call_esorex(calibration_dir+'DARK/',rootpath,esorex_cmd,n_CPU)
                
    if creating_sof == False:
        call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
        call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)
        if dark: call_esorex(calibration_dir+'DARK/',rootpath,esorex_cmd,n_CPU)

def dark(rootpath,working_dir,exposure_list,exposure_list_DARK,calibration_dir,creating_sof,n_CPU=24):
    
    print('... Creating the MASTER DARK')
    
    esorex_cmd = '--log-file=dark.log --log-level=debug muse_dark --nifu=-1 --merge dark.sof'
    
    if creating_sof:
        
        if os.path.exists(calibration_dir+'DARK/dark.sof') == True: os.remove(calibration_dir+'DARK/dark.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')
    
            raw_data_list = ascii.read(exposure_list_DARK[exposure_ID],format='no_header')
    
            f = open(calibration_dir+'DARK/dark_temp.sof', 'w')
            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'DARK': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            f.write(calibration_dir+'DARK/MASTER_BIAS.fits MASTER_BIAS\n')
            f.close()            
            
            if os.path.isfile(calibration_dir+'DARK/dark.sof'):
                if filecmp.cmp(calibration_dir+'DARK/dark.sof',calibration_dir+'DARK/dark_temp.sof'):
                    print('dark files are identical: proceed without error')
                    os.remove(calibration_dir+'DARK/dark_temp.sof')
                else:
                    sys.exit('CAUTION DARK FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'DARK/dark_temp.sof',calibration_dir+'DARK/dark.sof')
                call_esorex(calibration_dir+'DARK/',rootpath,esorex_cmd,n_CPU)
                
    if creating_sof == False: call_esorex(calibration_dir+'DARK/',rootpath,esorex_cmd,n_CPU)
        
def flat(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,creating_sof,n_CPU=24):
    print('... Creating the MASTER FLAT')
    
    esorex_cmd = '--log-file=flat.log --log-level=debug muse_flat --samples=true --nifu=-1 --merge flat.sof'
    
    if creating_sof:
        
        if os.path.exists(calibration_dir+'SCIENCE/flat.sof') == True: os.remove(calibration_dir+'SCIENCE/flat.sof')
        if os.path.exists(calibration_dir+'TWILIGHT/flat.sof') == True: os.remove(calibration_dir+'TWILIGHT/flat.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')      
        
            raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
            raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
            f = open(calibration_dir+'SCIENCE/flat_temp.sof', 'w')
            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'FLAT': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            f.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
            if dark: f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.close()

            if os.path.isfile(calibration_dir+'SCIENCE/flat.sof'):
                if filecmp.cmp(calibration_dir+'SCIENCE/flat.sof',calibration_dir+'SCIENCE/flat_temp.sof'):
                    print('science files are identical: proceed without error')
                    os.remove(calibration_dir+'SCIENCE/flat_temp.sof')
                else:
                    sys.exit('CAUTION SCIENCE FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'SCIENCE/flat_temp.sof',calibration_dir+'SCIENCE/flat.sof')
                call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
        
            f = open(calibration_dir+'TWILIGHT/flat_temp.sof', 'w')
            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'FLAT': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
            if dark: f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.close()
    
            if os.path.isfile(calibration_dir+'TWILIGHT/flat.sof'):
                if filecmp.cmp(calibration_dir+'TWILIGHT/flat.sof',calibration_dir+'TWILIGHT/flat_temp.sof'):
                    print('twilight files are identical: proceed without error')
                    os.remove(calibration_dir+'TWILIGHT/flat_temp.sof')
                else:
                    sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'TWILIGHT/flat_temp.sof',calibration_dir+'TWILIGHT/flat.sof')
                call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)

    if creating_sof == False:
        call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
        call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)
        
def wavecal(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,n_CPU=24):
    print('... Creating the WAVELENGTH CALIBRATION')
    
    esorex_cmd = '--log-file=wavecal.log --log-level=debug muse_wavecal --nifu=-1 --residuals --merge wavecal.sof'

    if creating_sof:
        
        if os.path.exists(calibration_dir+'SCIENCE/wavecal.sof') == True: os.remove(calibration_dir+'SCIENCE/wavecal.sof')
        if os.path.exists(calibration_dir+'TWILIGHT/wavecal.sof') == True: os.remove(calibration_dir+'TWILIGHT/wavecal.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')

            raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
            raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')

            f = open(calibration_dir+'SCIENCE/wavecal_temp.sof', 'w')
            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'ARC': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
            f.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'SCIENCE/TRACE_TABLE.fits TRACE_TABLE\n')
            if dark: f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.close()

            if os.path.isfile(calibration_dir+'SCIENCE/wavecal.sof'):
                if filecmp.cmp(calibration_dir+'SCIENCE/wavecal.sof',calibration_dir+'SCIENCE/wavecal_temp.sof'):
                    print('science files are identical: proceed without error')
                    os.remove(calibration_dir+'SCIENCE/wavecal_temp.sof')
                else:
                    sys.exit('CAUTION SCIENCE FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'SCIENCE/wavecal_temp.sof',calibration_dir+'SCIENCE/wavecal.sof')
                call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)

            f = open(calibration_dir+'TWILIGHT/wavecal_temp.sof', 'w')
            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'ARC': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
            if dark: f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.close()

            if os.path.isfile(calibration_dir+'TWILIGHT/wavecal.sof'):
                if filecmp.cmp(calibration_dir+'TWILIGHT/wavecal.sof',calibration_dir+'TWILIGHT/wavecal_temp.sof'):
                    print('twilight files are identical: proceed without error')
                    os.remove(calibration_dir+'TWILIGHT/wavecal_temp.sof')
                else:
                    sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'TWILIGHT/wavecal_temp.sof',calibration_dir+'TWILIGHT/wavecal.sof')
                call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)
                
    if creating_sof == False:
        call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
        call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)
    
def lsf(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,n_CPU=24):
    print('... Creating the LINE SPREAD FUNCTION')
    
    esorex_cmd = '--log-file=lsf.log --log-level=debug muse_lsf --nifu=-1 --merge --save_subtracted lsf.sof'
    
    if creating_sof:
        
        if os.path.exists(calibration_dir+'SCIENCE/lsf.sof') == True: os.remove(calibration_dir+'SCIENCE/lsf.sof')
        if os.path.exists(calibration_dir+'TWILIGHT/lsf.sof') == True: os.remove(calibration_dir+'TWILIGHT/lsf.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')
    
            raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
            raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
            f = open(calibration_dir+'SCIENCE/lsf_temp.sof', 'w')
            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'ARC': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
            f.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'SCIENCE/TRACE_TABLE.fits TRACE_TABLE\n')
            f.write(calibration_dir+'SCIENCE/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            if dark: f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.write(calibration_dir+'SCIENCE/MASTER_FLAT.fits MASTER_FLAT\n')
            f.close()

            if os.path.isfile(calibration_dir+'SCIENCE/lsf.sof'):
                if filecmp.cmp(calibration_dir+'SCIENCE/lsf.sof',calibration_dir+'SCIENCE/lsf_temp.sof'):
                    print('science files are identical: proceed without error')
                    os.remove(calibration_dir+'SCIENCE/lsf_temp.sof')
                else:
                    sys.exit('CAUTION SCIENCE FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'SCIENCE/lsf_temp.sof',calibration_dir+'SCIENCE/lsf.sof')
                call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
    
            f = open(calibration_dir+'TWILIGHT/lsf_temp.sof', 'w')
            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'ARC': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
            f.write(calibration_dir+'TWILIGHT/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            if dark: f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_FLAT.fits MASTER_FLAT\n')
            f.close()

            if os.path.isfile(calibration_dir+'TWILIGHT/lsf.sof'):
                if filecmp.cmp(calibration_dir+'TWILIGHT/lsf.sof',calibration_dir+'TWILIGHT/lsf_temp.sof'):
                    print('twilight files are identical: proceed without error')
                    os.remove(calibration_dir+'TWILIGHT/lsf_temp.sof')
                else:
                    sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'TWILIGHT/lsf_temp.sof',calibration_dir+'TWILIGHT/lsf.sof')
                call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)

    if creating_sof == False:
        call_esorex(calibration_dir+'SCIENCE/',rootpath,esorex_cmd,n_CPU)
        call_esorex(calibration_dir+'TWILIGHT/',rootpath,esorex_cmd,n_CPU)

def twilight(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,mode,n_CPU=24):
    print('... Creating the TWILIGHT FLAT')
    
    esorex_cmd = '--log-file=twilight.log --log-level=debug muse_twilight twilight.sof'

    if creating_sof:
        
        if os.path.exists(calibration_dir+'TWILIGHT/twilight.sof') == True: os.remove(calibration_dir+'TWILIGHT/twilight.sof')
        
        for exposure_ID in range(len(exposure_list)):
            print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
            print('>>> processing: '+exposure_list[exposure_ID])
            print(' ')   
               
            raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
            raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
            
            MJDsillum=np.array([])
            MJDsskyflat=np.array([])
            illum_index=np.array([])

            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'SKYFLAT':
                    skyflathdu=fits.open(raw_data_list_TWILIGHT[i][0])
                    MJDskyflat=skyflathdu[0].header['MJD-OBS']
                    MJDsskyflat=np.append(MJDsskyflat,MJDskyflat)
                    
                if raw_data_list_TWILIGHT[i][1] == 'ILLUM':
                    illumhdu=fits.open(raw_data_list_TWILIGHT[i][0])
                    MJDillum=illumhdu[0].header['MJD-OBS']
                    MJDsillum=np.append(MJDsillum,MJDillum)
                    illum_index=np.append(illum_index,i)
                
            choosen_illum=int(illum_index[np.argmin(np.abs(MJDsillum-np.min(MJDskyflat)))])

            f = open(calibration_dir+'TWILIGHT/twilight_temp.sof', 'w')
            for i in range(len(raw_data_list_TWILIGHT[1][:])):
                if raw_data_list_TWILIGHT[i][1] == 'SKYFLAT': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            f.write(raw_data_list_TWILIGHT[choosen_illum][0]+'  '+raw_data_list_TWILIGHT[choosen_illum][1]+'\n')
            if mode == 'WFM-AO' or mode == 'WFM-NOAO': f.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
            if mode == 'NFM-AO':f.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'TWILIGHT/MASTER_FLAT.fits MASTER_FLAT\n')
            if dark: f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
            f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
            f.write(calibration_dir+'TWILIGHT/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            if MJDsskyflat.all() < 57823.5: f.write(static_calibration_dir+'vignetting_mask.fits VIGNETTING_MASK\n')
        
            f.close()

            if os.path.isfile(calibration_dir+'TWILIGHT/twilight.sof'):
                if filecmp.cmp(calibration_dir+'TWILIGHT/twilight.sof',calibration_dir+'TWILIGHT/twilight_temp.sof'):
                    print('twilight files are identical: proceed without error')
                    os.remove(calibration_dir+'TWILIGHT/twilight_temp.sof')
                else:
                    sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
            else:
                os.rename(calibration_dir+'TWILIGHT/twilight_temp.sof',calibration_dir+'TWILIGHT/twilight.sof')
                call_esorex(calibration_dir+'TWILIGHT',rootpath,esorex_cmd,n_CPU)
                
    if creating_sof == False:
        call_esorex(calibration_dir+'TWILIGHT',rootpath,esorex_cmd,n_CPU)

###OBSERVATION PRE-PROCESSING###

def science_pre(rootpath,working_dir,exposure_list,dark,calibration_dir,ESO_calibration_dir,static_calibration_dir,skyreject,using_ESO_calibration,creating_sof,mode,n_CPU=24):
    print('... Science PREPROCESSING')
    
    esorex_cmd = '--log-file=sci_basic_object.log --log-level=debug muse_scibasic --nifu=-1 --resample --saveimage=true --skyreject='+skyreject+' --merge  sci_basic_object.sof'
    esorex_cmd_std = '--log-file=sci_basic_std.log --log-level=debug muse_scibasic --nifu=-1 --resample --saveimage=true --skyreject=15.,15.,1 --merge  sci_basic_std.sof'
    
    
    if os.path.exists(working_dir+'std/sci_basic_std.sof') == True: os.remove(working_dir+'std/sci_basic_std.sof')
        
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')

        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        exposure_dir=exposure_list[exposure_ID][:-9]+'/'

        MJDsillum=np.array([])
        illum_index=np.array([])

        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'OBJECT' or raw_data_list[i][1] == 'SKY':
                objecthdu=fits.open(raw_data_list[i][0])
                MJDobject=objecthdu[0].header['MJD-OBS']
                
            if raw_data_list[i][1] == 'STD':
                stdhdu=fits.open(raw_data_list[i][0])
                MJDstd=stdhdu[0].header['MJD-OBS']
                
            if raw_data_list[i][1] == 'ILLUM':
                illumhdu=fits.open(raw_data_list[i][0])
                MJDillum=illumhdu[0].header['MJD-OBS']
                MJDsillum=np.append(MJDsillum,MJDillum)
                illum_index=np.append(illum_index,i)
                
        choosen_illum_object=int(illum_index[np.argmin(np.abs(MJDsillum-MJDobject))])
        choosen_illum_std=int(illum_index[np.argmin(np.abs(MJDsillum-MJDstd))])
            
        f_std = open(working_dir+'std/sci_basic_std_temp.sof', 'w')
        
        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'STD': f_std.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')

        f_std.write(raw_data_list[choosen_illum_std][0]+'  '+raw_data_list[choosen_illum_std][1]+'\n')
        
        if using_ESO_calibration == False:
            f_std.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
            f_std.write(calibration_dir+'SCIENCE/MASTER_FLAT.fits MASTER_FLAT\n')
            f_std.write(calibration_dir+'TWILIGHT/TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
            f_std.write(calibration_dir+'SCIENCE/TRACE_TABLE.fits TRACE_TABLE\n')
            f_std.write(calibration_dir+'SCIENCE/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            if dark: f_std.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')

        if using_ESO_calibration == True:
            f_std.write(ESO_calibration_dir+'MASTER_BIAS.fits MASTER_BIAS\n')
            f_std.write(ESO_calibration_dir+'MASTER_FLAT.fits MASTER_FLAT\n')
            f_std.write(ESO_calibration_dir+'TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
            f_std.write(ESO_calibration_dir+'TRACE_TABLE.fits TRACE_TABLE\n')
            f_std.write(ESO_calibration_dir+'WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            if dark: f_std.write(ESO_calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')

        if mode == 'WFM-AO' or mode == 'WFM-NOAO': f_std.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
        if mode == 'NFM-AO': f_std.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')

        f_std.write(static_calibration_dir+'badpix_table.fits BADPIX_TABLE\n')
        f_std.close()
        
        if creating_sof:
            
            if os.path.exists(exposure_dir+'sci_basic_object.sof') == True: os.remove(exposure_dir+'sci_basic_object.sof')
            f_object = open(exposure_dir+'sci_basic_object.sof', 'w')

            for i in range(len(raw_data_list[1][:])):
                if raw_data_list[i][1] == 'OBJECT': f_object.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
                if raw_data_list[i][1] == 'SKY': f_object.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')

            f_object.write(raw_data_list[choosen_illum_object][0]+'  '+raw_data_list[choosen_illum_object][1]+'\n')
            
            if using_ESO_calibration == False:
                f_object.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
                f_object.write(calibration_dir+'SCIENCE/MASTER_FLAT.fits MASTER_FLAT\n')
                f_object.write(calibration_dir+'TWILIGHT/TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
                f_object.write(calibration_dir+'SCIENCE/TRACE_TABLE.fits TRACE_TABLE\n')
                f_object.write(calibration_dir+'SCIENCE/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
                if dark: f_object.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')

            if using_ESO_calibration == True:
                f_object.write(ESO_calibration_dir+'MASTER_BIAS.fits MASTER_BIAS\n')            
                f_object.write(ESO_calibration_dir+'MASTER_FLAT.fits MASTER_FLAT\n')
                f_object.write(ESO_calibration_dir+'TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
                f_object.write(ESO_calibration_dir+'TRACE_TABLE.fits TRACE_TABLE\n')
                f_object.write(ESO_calibration_dir+'WAVECAL_TABLE.fits WAVECAL_TABLE\n')
                if dark: f_object.write(ESO_calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')

            if mode == 'WFM-AO' or mode == 'WFM-NOAO': f_object.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
            if mode == 'NFM-AO': f_object.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
            f_object.write(static_calibration_dir+'badpix_table.fits BADPIX_TABLE\n')
    
            f_object.close()
    
        call_esorex(exposure_dir,rootpath,esorex_cmd,n_CPU)
        if os.path.isfile(working_dir+'std/sci_basic_std.sof'):
            if filecmp.cmp(working_dir+'std/sci_basic_std.sof',working_dir+'std/sci_basic_std_temp.sof'):
                os.remove(working_dir+'std/sci_basic_std_temp.sof')
            else:
                sys.exit('CAUTION DIFFERENT STD STARS FOR VARIOUS FIELDS: PLEASE CHECK')
        else:
            os.rename(working_dir+'std/sci_basic_std_temp.sof',working_dir+'std/sci_basic_std.sof')
    
    call_esorex(working_dir+'std/',rootpath,esorex_cmd_std,n_CPU)

### OBSERVATION POST-PROCESSING ###

def std_flux(rootpath,working_dir,exposure_list,calibration_dir,static_calibration_dir,creating_sof,n_CPU=24):
    print('... FLUX CALIBRATION')
    
    esorex_cmd = ' --log-file=std_flux.log --log-level=debug muse_standard --filter=white std_flux.sof'
    
    PIXTABLE_STD_list=get_filelist(working_dir+'std/',rootpath,'PIXTABLE_STD*.fits')
    if creating_sof:
        
        if os.path.exists(working_dir+'std/'+'std_flux.sof') == True: os.remove(working_dir+'std/'+'std_flux.sof')
        
        f = open(working_dir+'std/'+'std_flux.sof', 'w')
        for i in range(len(PIXTABLE_STD_list)): f.write(working_dir+'std/'+PIXTABLE_STD_list[i]+' PIXTABLE_STD\n')
        f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
        f.write(static_calibration_dir+'std_flux_table.fits STD_FLUX_TABLE\n')
        f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
        f.close()
    
    call_esorex(working_dir+'std/',rootpath,esorex_cmd,n_CPU)
    
def sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,creating_sof,skyfield,skyignore,skyfraction,n_CPU=24):
    print('... SKY CREATION')
    
    sky = np.zeros_like(exposure_list,dtype=bool)
    sci = np.zeros_like(exposure_list,dtype=bool)
    
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SKY': sky[idx] = True
        if exposure[-8:-5] == 'SCI': sci[idx] = True

    if skyfield == 'auto' and (sky == True).any():
        exposure_list_sky = np.array(exposure_list)[sky]
    else:
        exposure_list_sky = np.array(exposure_list)
    
    for exposure_ID in range(len(exposure_list_sky)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list_sky)))
        print('>>> processing: '+exposure_list_sky[exposure_ID])
        print(' ')
    
        raw_data_list = ascii.read(exposure_list_sky[exposure_ID],format='no_header')
        exposure_dir=exposure_list_sky[exposure_ID][:-9]+'/'

        if skyfield == 'auto' and (sky == True).any():
            PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_SKY*.fits')
        else:
            PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_OBJECT*.fits')

        if creating_sof:
            
            if os.path.exists(exposure_dir+'sky.sof') == True: os.remove(exposure_dir+'sky.sof')
            
            f = open(exposure_dir+'sky.sof', 'w')
            for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
            if using_ESO_calibration == False:
                f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
            if using_ESO_calibration == True:
                f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
            f.write(working_dir+'std/'+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
            f.write(working_dir+'std/'+'STD_TELLURIC_0001.fits STD_TELLURIC\n')

            f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
            f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
            f.close()
        esorex_cmd = "--log-file=sky.log --log-level=debug muse_create_sky --fraction="+str(skyfraction)+" --ignore="+str(skyignore)+" sky.sof"
        if skyfield == 'auto' and (sky == True).any(): call_esorex(exposure_dir,rootpath,esorex_cmd,n_CPU)
        else: call_esorex(exposure_dir,rootpath,'--log-file=sky.log --log-level=debug muse_create_sky --fraction='+str(skyfraction)+' --ignore='+str(skyignore)+' sky.sof',n_CPU)
        
    if skyfield == 'auto' and (sky == True).any():
        skydate = np.ones_like(exposure_list_sky,dtype=float)
        
        for idx,exps in enumerate(exposure_list_sky):
            skydate[idx] = fits.open(exps[:-9]+'/PIXTABLE_SKY_0001-01.fits')[0].header['MJD-OBS']
        for idx,exps in enumerate(np.array(exposure_list)[sci]):
            scidate = fits.open(exps[:-9]+'/PIXTABLE_OBJECT_0001-01.fits')[0].header['MJD-OBS']
            
            ind = np.argmin(abs(skydate - scidate))
            flist = glob.glob(exposure_list_sky[ind][:-9]+'/SKY_*.fits')
            for f in flist: shutil.copy(f,exps[:-9]+'/.')
            
def modified_sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,creating_sof,skyfield,skyignore,skyfraction,n_CPU=24):

    print('... SKY CREATION MODIFIED')
    
    sky = np.zeros_like(exposure_list,dtype=bool)
    sci = np.zeros_like(exposure_list,dtype=bool)
    
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SKY': sky[idx] = True
        if exposure[-8:-5] == 'SCI': sci[idx] = True

    if skyfield == 'auto' and (sky == True).any():
        exposure_list_sky = np.array(exposure_list)[sky]
    else:
        exposure_list_sky = np.array(exposure_list)
    
    for exposure_ID in range(len(exposure_list_sky)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list_sky)))
        print('>>> processing: '+exposure_list_sky[exposure_ID])
        print(' ')
    
        raw_data_list = ascii.read(exposure_list_sky[exposure_ID],format='no_header')
        exposure_dir=exposure_list_sky[exposure_ID][:-9]+'/'

        if skyfield == 'auto' and (sky == True).any():
            PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_SKY*.fits')
        else:
            PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_OBJECT*.fits')
        
        if creating_sof:
            
            if os.path.exists(exposure_dir+'sky.sof') == True: os.remove(exposure_dir+'sky.sof')
            
            f = open(exposure_dir+'sky.sof', 'w')
            for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
            if using_ESO_calibration == False:
                f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
            if using_ESO_calibration == True:
                f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
            f.write(working_dir+'std/'+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
            f.write(working_dir+'std/'+'STD_TELLURIC_0001.fits STD_TELLURIC\n')

            f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
            f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
            f.close()

        if skyfield == 'auto' and (sky == True).any(): call_esorex(exposure_dir,rootpath,'--log-file=sky.log --log-level=debug muse_create_sky --fraction='+str(skyfraction)+' --ignore='+str(skyignore)+' sky.sof',n_CPU)
        else: call_esorex(exposure_dir,rootpath,'--log-file=sky.log --log-level=debug muse_create_sky --fraction='+str(skyfraction)+' --ignore='+str(skyignore)+' sky.sof',n_CPU)
        
        ### continuum set to zero
        
        os.chdir(exposure_dir)
        sky_cont_hdu=fits.open('SKY_CONTINUUM.fits',checksum=True)
        sky_cont=sky_cont_hdu[1].data
        for i in range(len(sky_cont)):
            sky_cont[i][1] = 0.
        sky_cont_hdu[1].data = sky_cont
        sky_cont_hdu.writeto('SKY_CONTINUUM_zero.fits',overwrite = True,checksum=True)
        os.chdir(rootpath)
        
        if creating_sof:
            
            if os.path.exists(exposure_dir+'sky.sof') == True: os.remove(exposure_dir+'sky.sof')
            
            f = open(exposure_dir+'sky.sof', 'w')
            for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
            if using_ESO_calibration == False:
                f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
            if using_ESO_calibration == True:
                f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
            f.write(working_dir+'std/'+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
            f.write(working_dir+'std/'+'STD_TELLURIC_0001.fits STD_TELLURIC\n')
        
            f.write(exposure_dir+'SKY_CONTINUUM_zero.fits SKY_CONTINUUM\n')

            f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
            f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
            f.close()
        
        if skyfield == 'auto' and (sky == True).any(): call_esorex(exposure_dir,rootpath,'--log-file=sky.log --log-level=debug muse_create_sky --fraction='+str(skyfraction)+' --ignore='+str(skyignore)+' sky.sof',n_CPU)
        else: call_esorex(exposure_dir,rootpath,'--log-file=sky.log --log-level=debug muse_create_sky --fraction='+str(skyfraction)+' --ignore='+str(skyignore)+' sky.sof',n_CPU)
        
        os.chdir(exposure_dir)
        print('SKY_CONTINUUM_zero.fits ==> SKY_CONTINUUM.fits')
        shutil.copy('SKY_CONTINUUM_zero.fits','SKY_CONTINUUM.fits')
        hdu=fits.open('SKY_LINES.fits',checksum=True)
        data=hdu[1].data

        lines_to_keep=[i for i,sky_line in enumerate(data) if (sky_line[0][:2] == 'O2' or sky_line[0][:2] == 'OH')]
        data=data[lines_to_keep]

        hdu[1].data = data
        hdu.writeto('SKY_LINES.fits', overwrite = True, checksum=True)
        os.chdir(rootpath)
        
    if skyfield == 'auto' and (sky == True).any():
        skydate = np.ones_like(exposure_list_sky,dtype=float)
    
        for idx,exps in enumerate(exposure_list_sky):
            skydate[idx] = fits.open(exps[:-9]+'/PIXTABLE_SKY_0001-01.fits')[0].header['MJD-OBS']
        for idx,exps in enumerate(np.array(exposure_list)[sci]):
            scidate = fits.open(exps[:-9]+'/PIXTABLE_OBJECT_0001-01.fits')[0].header['MJD-OBS']
        
            ind = np.argmin(abs(skydate - scidate))
            flist = glob.glob(exposure_list_sky[ind][:-9]+'/SKY_*.fits')
            for f in flist: shutil.copy(f,exps[:-9]+'/.')
        
### SCIENCE POST-PROCESSING ###

def scipost(rootpath,working_dir,static_calibration_dir,exposure_list,calibration_dir,ESO_calibration_dir,using_ESO_calibration,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,OB,creating_sof,raman,mode,n_CPU=24):
    print('... SCIENCE POST-PROCESSING')
    
    unique_pointings=np.array([])
    unique_tester=' '
    
    sci = np.zeros_like(exposure_list,dtype=bool)
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SCI': sci[idx] = True
    exposure_list = np.array(exposure_list)[sci]
    
    for expnum in range(len(exposure_list)):
        if unique_tester.find(exposure_list[expnum][:-16])==-1:
            unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-16])
            unique_tester=unique_tester+exposure_list[expnum][:-16]
    
    for unique_pointing_num in range(len(unique_pointings)):
        
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')
        
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]
        exp_list=glob.glob(sec+'*SCI.list')

        for exp_num in range(len(exp_list)):
            
            print(' ')
            print('>>> processing exposure: '+str(exp_num+1)+'/'+str(len(exp_list)))
            print(' ')
            
            PIXTABLE_OBJECT_list=get_filelist(exp_list[exp_num][:-9],rootpath,'PIXTABLE_OBJECT*.fits')
            
            if creating_sof:

                if os.path.exists(exp_list[exp_num][:-9]+'/scipost.sof') == True: os.remove(exp_list[exp_num][:-9]+'/scipost.sof')
                
                f = open(exp_list[exp_num][:-9]+'/scipost.sof', 'w')
                for i in range(len(PIXTABLE_OBJECT_list)): f.write(exp_list[exp_num][:-9]+'/'+PIXTABLE_OBJECT_list[i]+' PIXTABLE_OBJECT\n')
            
                if using_ESO_calibration == False: f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
                if using_ESO_calibration == True: f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
                
                f.write(working_dir+'std/'+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
                f.write(working_dir+'std/'+'STD_TELLURIC_0001.fits STD_TELLURIC\n')
                if skysub: f.write(exp_list[exp_num][:-9]+'/'+'SKY_LINES.fits SKY_LINES\n')
                if skysub: f.write(exp_list[exp_num][:-9]+'/'+'SKY_CONTINUUM.fits SKY_CONTINUUM\n')
                if mode =='WFM-AO' or mode =='WFM-NOAO': f.write(static_calibration_dir+'astrometry_wcs_wfm.fits ASTROMETRY_WCS\n')
                if mode =='NFM-AO': print("CURRENTLY NO ASTRONOMY_WCS_NFM AVAILABLE !!!!")

                f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
                f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                if raman == True: f.write(static_calibration_dir+'raman_lines.fits RAMAN_LINES\n')
                    
                f.close()
            if withrvcorr:
                if skysub:
                    print('without sky ...')
                    
                    if raman == True: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual,raman --skymethod=subtract-model --filter=white scipost.sof',n_CPU)
                    if raman == False: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=subtract-model --filter=white scipost.sof',n_CPU)
                    
                    os.chdir(exp_list[exp_num][:-9])
                    os.rename('DATACUBE_FINAL.fits','DATACUBE_FINAL_wosky.fits')
                    os.rename('IMAGE_FOV_0001.fits','IMAGE_FOV_0001_wosky.fits')
                    os.rename('PIXTABLE_REDUCED_0001.fits','PIXTABLE_REDUCED_0001_wosky.fits')
                    os.chdir(rootpath)
                    
                if not skysub:
                    print('with sky ...')
                    
                    if raman == True: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual,raman --skymethod=none --filter=white scipost.sof',n_CPU)
                    if raman == False: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=none --filter=white scipost.sof',n_CPU)
                    
                    os.chdir(exp_list[exp_num][:-9])
                    os.rename('DATACUBE_FINAL.fits','DATACUBE_FINAL_wsky.fits')
                    os.rename('IMAGE_FOV_0001.fits','IMAGE_FOV_0001_wsky.fits')
                    os.rename('PIXTABLE_REDUCED_0001.fits','PIXTABLE_REDUCED_0001_wsky.fits')
                    os.chdir(rootpath)
            else:
                
                if raman == True: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual,raman --skymethod=none --rvcorr=none --filter=white scipost.sof',n_CPU)
                if raman == False: call_esorex(exp_list[exp_num][:-9],rootpath,'--log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=none --rvcorr=none --filter=white scipost.sof',n_CPU)
                
                os.chdir(exp_list[exp_num][:-9])
                os.rename('DATACUBE_FINAL.fits','DATACUBE_FINAL_wskynorvcorr.fits')
                os.rename('IMAGE_FOV_0001.fits','IMAGE_FOV_0001_wskynorvcorr.fits')
                os.rename('PIXTABLE_REDUCED_0001.fits','PIXTABLE_REDUCED_0001_wskynorvcorr.fits')
                os.chdir(rootpath)

def dither_collect(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,user_list,OB):
    
    print('... COLLECT DITHER POSTITIONS')
        
    unique_pointings=np.array([])
    unique_tester=' '
    
    sci = np.zeros_like(exposure_list,dtype=bool)
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SCI': sci[idx] = True
    exposure_list = np.array(exposure_list)[sci]
        
    print(' ')
    print('>>> Copying files:')
    print(' ')

    if len(user_list) == 0:
        for expnum in range(len(exposure_list)):
                if unique_tester.find(exposure_list[expnum][:-16]) == -1:
                    unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-16])
                    unique_tester=unique_tester+exposure_list[expnum][:-16]
    
    if len(user_list) > 0: unique_pointings = np.array([user_list[0]+'_usr'])
            
    for unique_pointing_num in range(len(unique_pointings)):
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]
            
            
        if len(user_list) == 0: exp_list=glob.glob(sec+'*SCI.list')
        if len(user_list) > 0: exp_list = user_list+'_SCI.list'
        
        if dithering_multiple_OBs:
            if withrvcorr:
                combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
        
        if dithering_multiple_OBs == False:
            if withrvcorr:
                combining_exposure_dir_withoutsky=sec+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=sec+'/withsky_withrvcorr'
            else: combining_exposure_dir=sec+'/withsky_withoutrvcorr'
                
        if not os.path.exists(combining_exposure_dir_withoutsky): os.makedirs(combining_exposure_dir_withoutsky)
        if not os.path.exists(combining_exposure_dir_withsky): os.makedirs(combining_exposure_dir_withsky)
        
        if withrvcorr:
            files = glob.glob(combining_exposure_dir_withoutsky+'/*FOV_0001*')
            if len(files) > 0:
                for f in files: os.remove(f)
            files = glob.glob(combining_exposure_dir_withsky+'/*FOV_0001*')
            if len(files) > 0:
                for f in files: os.remove(f)
        if not withrvcorr:
            files = glob.glob(combining_exposure_dir+'/*FOV_0001*')
            if len(files) > 0:
                for f in files: os.remove(f)
            
            
    for unique_pointing_num in range(len(unique_pointings)):
    
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]

        if len(user_list) == 0: exp_list=glob.glob(sec+'*SCI.list')
        if len(user_list) > 0: exp_list = user_list+'_SCI.list'
        
        if dithering_multiple_OBs:
            if withrvcorr:
                combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
        
        if dithering_multiple_OBs == False:
            if withrvcorr:
                combining_exposure_dir_withoutsky=sec+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=sec+'/withsky_withrvcorr'
            else: combining_exposure_dir=sec+'/withsky_withoutrvcorr'
        
        for exp_num in range(len(exp_list)):
            if withrvcorr:
                if skysub:
                    if dithering_multiple_OBs:
                        print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wosky.fits',combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wosky.fits',combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+OB+'.fits')
                        print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wosky.fits',combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                    else:
                        print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wosky.fits',combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wosky.fits',combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wosky.fits ==> '+combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wosky.fits',combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')

                if not skysub:
                    if dithering_multiple_OBs:
                        print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wsky.fits ==> '+combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wsky.fits',combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wsky.fits ==> '+combining_exposure_dir_withsky+'/IMAGE_FOV_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wsky.fits',combining_exposure_dir_withsky+'/IMAGE_FOV_'+OB+'.fits')
                        print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wsky.fits ==> '+combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wsky.fits',combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                        
                    else:
                        print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wsky.fits ==> '+combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wsky.fits',combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wsky.fits ==> '+combining_exposure_dir_withsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wsky.fits',combining_exposure_dir_withsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wsky.fits ==> '+combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wsky.fits',combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
            else:
                
                if dithering_multiple_OBs:
                    print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wskynorvcorr.fits ==> '+combining_exposure_dir+'/DATACUBE_FINAL_'+OB+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wskynorvcorr.fits',combining_exposure_dir+'/DATACUBE_FINAL_'+OB+'.fits')
                    print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wskynorvcorr.fits ==> '+combining_exposure_dir+'/IMAGE_FOV_'+OB+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wskynorvcorr.fits',combining_exposure_dir+'/IMAGE_FOV_'+OB+'.fits')
                    print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wskynorvcorr.fits ==> '+combining_exposure_dir+'/PIXTABLE_REDUCED_'+OB+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wskynorvcorr.fits',combining_exposure_dir+'/PIXTABLE_REDUCED_'+OB+'.fits')
                else:
                    print(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wskynorvcorr.fits ==> '+combining_exposure_dir+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/DATACUBE_FINAL_wskynorvcorr.fits',combining_exposure_dir+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    print(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wskynorvcorr.fits ==> '+combining_exposure_dir+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/IMAGE_FOV_0001_wskynorvcorr.fits',combining_exposure_dir+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    print(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wskynorvcorr.fits ==> '+combining_exposure_dir+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.copy(exp_list[exp_num][:-9]+'/PIXTABLE_REDUCED_0001_wskynorvcorr.fits',combining_exposure_dir+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
    
def exp_align(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,creating_sof,user_list,OB,n_CPU=24):
    print('... CUBE ALIGNMENT')
    
    esorex_cmd = '--log-file=exp_align.log --log-level=debug muse_exp_align exp_align.sof'
    
    unique_pointings=np.array([])
    unique_tester=' '
    
    sci = np.zeros_like(exposure_list,dtype=bool)
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SCI': sci[idx] = True
    exposure_list = np.array(exposure_list)[sci]
    
    if len(user_list) == 0:
        for expnum in range(len(exposure_list)):
                if unique_tester.find(exposure_list[expnum][:-16]) == -1:
                    unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-16])
                    unique_tester=unique_tester+exposure_list[expnum][:-16]
    
    if len(user_list) > 0: unique_pointings = np.array([user_list[0]+'_usr'])
    
            
    for unique_pointing_num in range(len(unique_pointings)):
    
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]
        
        
        if len(user_list) == 0: exp_list=glob.glob(sec+'*SCI.list')
        if len(user_list) > 0: exp_list = user_list+'_SCI.list'
        
        
        if dithering_multiple_OBs:
            if withrvcorr:
                combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
        
        if dithering_multiple_OBs == False:
            if withrvcorr:
                combining_exposure_dir_withoutsky=sec+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=sec+'/withsky_withrvcorr'
            else: combining_exposure_dir=sec+'/withsky_withoutrvcorr'
        
    
    for unique_pointing_num in range(len(unique_pointings)):
    
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')

        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]
    
        if dithering_multiple_OBs:
            if withrvcorr:
                print(unique_pointings_ID)
                if skysub: combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                if not skysub: combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
        
        if dithering_multiple_OBs == False:
            if withrvcorr:
                if skysub: combining_exposure_dir_withoutsky=sec+'/withoutsky_withrvcorr'
                if not skysub: combining_exposure_dir_withsky=sec+'/withsky_withrvcorr'
            else: combining_exposure_dir=sec+'/withsky_withoutrvcorr'
    
        if withrvcorr:
            if skysub:
                exp_list=get_filelist(combining_exposure_dir_withoutsky,rootpath,'IMAGE_FOV_*.fits')
                if creating_sof:
                    if os.path.exists(combining_exposure_dir_withoutsky+'/exp_align.sof') == True: os.remove(combining_exposure_dir_withoutsky+'/exp_align.sof')
                
                    f = open(combining_exposure_dir_withoutsky+'/exp_align.sof', 'w')
                    for i in range(len(exp_list)): f.write(combining_exposure_dir_withoutsky+'/'+exp_list[i]+' IMAGE_FOV\n')
                    f.close()
                
                call_esorex(combining_exposure_dir_withoutsky,rootpath,esorex_cmd,n_CPU)
            
            if not skysub:
                exp_list=get_filelist(combining_exposure_dir_withsky,rootpath,'IMAGE_FOV_*.fits')
                if creating_sof:
                
                    if os.path.exists(combining_exposure_dir_withsky+'/exp_align.sof') == True: os.remove(combining_exposure_dir_withsky+'/exp_align.sof')
                
                    f = open(combining_exposure_dir_withsky+'/exp_align.sof', 'w')
                    for i in range(len(exp_list)): f.write(combining_exposure_dir_withsky+'/'+exp_list[i]+' IMAGE_FOV\n')
                    f.close()
                
                call_esorex(combining_exposure_dir_withsky,rootpath,esorex_cmd,n_CPU)
        
        else:
            exp_list=get_filelist(combining_exposure_dir,rootpath,'IMAGE_FOV_*.fits')
            if creating_sof:
                if os.path.exists(combining_exposure_dir+'/exp_align.sof') == True: os.remove(combining_exposure_dir+'/exp_align.sof')
            
                f = open(combining_exposure_dir+'/exp_align.sof', 'w')
                for i in range(len(exp_list)): f.write(combining_exposure_dir+'/'+exp_list[i]+' IMAGE_FOV\n')
                f.close()
        
            call_esorex(combining_exposure_dir,rootpath,esorex_cmd,n_CPU)
          
def exp_combine(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,static_calibration_dir,combining_OBs_dir,creating_sof,user_list,n_CPU=24):
    print('... EXPOSURE COMBINATION')
    
    esorex_cmd = '--log-file=exp_combine.log --log-level=debug muse_exp_combine --filter=white --save=cube --crsigma=5. exp_combine.sof'
    
    unique_pointings=np.array([])
    unique_tester=' '
    
    sci = np.zeros_like(exposure_list,dtype=bool)
    for idx,exposure in enumerate(exposure_list):
        if exposure[-8:-5] == 'SCI': sci[idx] = True
    exposure_list = np.array(exposure_list)[sci]
    
    if len(user_list) == 0:
        for expnum in range(len(exposure_list)):
                if unique_tester.find(exposure_list[expnum][:-16]) == -1:
                    unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-16])
                    unique_tester=unique_tester+exposure_list[expnum][:-16]
    
    if len(user_list) > 0: unique_pointings = np.array([user_list[0]+'_usr'])
    
    for unique_pointing_num in range(len(unique_pointings)):
        
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')
        
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        sec=unique_pointings[unique_pointing_num]
                
        if dithering_multiple_OBs:
            if withrvcorr:
                if skysub: combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                if not skysub: combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
            
        if dithering_multiple_OBs == False:
            if withrvcorr:
                if skysub: combining_exposure_dir_withoutsky=sec+'/withoutsky_withrvcorr'
                if not skysub: combining_exposure_dir_withsky=sec+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=sec+'/withsky_withoutrvcorr'
        
        if withrvcorr:
            if skysub:
                pixtable_list=get_filelist(combining_exposure_dir_withoutsky,rootpath,'PIXTABLE_REDUCED_*.fits')
                if creating_sof:
                    if os.path.exists(combining_exposure_dir_withoutsky+'/exp_combine.sof') == True: os.remove(combining_exposure_dir_withoutsky+'/exp_combine.sof')
                    
                    f = open(combining_exposure_dir_withoutsky+'/exp_combine.sof', 'w')
                    for i in range(len(pixtable_list)): f.write(combining_exposure_dir_withoutsky+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
                    f.write(combining_exposure_dir_withoutsky+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
                    f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                    f.close()
                
                call_esorex(combining_exposure_dir_withoutsky,rootpath,esorex_cmd,n_CPU)

            if not skysub:
                pixtable_list=get_filelist(combining_exposure_dir_withsky,rootpath,'PIXTABLE_REDUCED_*.fits')
                if creating_sof:
                    if os.path.exists(combining_exposure_dir_withsky+'/exp_combine.sof') == True: os.remove(combining_exposure_dir_withsky+'/exp_combine.sof')
                    
                    f = open(combining_exposure_dir_withsky+'/exp_combine.sof', 'w')
                    for i in range(len(pixtable_list)): f.write(combining_exposure_dir_withsky+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
                    f.write(combining_exposure_dir_withsky+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
                    f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                    f.close()
                
                call_esorex(combining_exposure_dir_withsky,rootpath,esorex_cmd,n_CPU)
            
        else:
            pixtable_list=get_filelist(combining_exposure_dir,rootpath,'PIXTABLE_REDUCED_*.fits')
            if creating_sof:
                if os.path.exists(combining_exposure_dir+'/exp_combine.sof') == True: os.remove(combining_exposure_dir+'/exp_combine.sof')
                
                f = open(combining_exposure_dir+'/exp_combine.sof', 'w')
                for i in range(len(pixtable_list)): f.write(combining_exposure_dir+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
                f.write(combining_exposure_dir+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
                f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                f.close()
                
            call_esorex(combining_exposure_dir,rootpath,esorex_cmd,n_CPU)
    
# if __name__=='__main__':
def musereduce(configfile=None):
    
    
    print('##########################################################################################')
    print('##########################################################################################')
    print('### #### ##     ##   ####   #### ### #### ### ###### #### ###### ##     ##################')
    print('###  ##  ## ###### ##  ## ##  ### # ###### # #######  ##  #####  ## ######################')
    print('### #  # ##   ####    ###    ##### ######## ###   ## #### #### # ##     ##################')
    print('### #### ## ###### ## ### ## ##### ####### # ####### #### ###    ###### ##################')
    print('### #### ##     ## ### ## ### #### ###### ### ###### #### ## ### ##     ##################')
    print('##########################################################################################')
    print('##########################################################################################')
    print('##########################################################################################')
    
    time.sleep(5)
    print(' ')
    print('##########################################################################################')
    print('#####                                                                                #####')
    print('#####                   MUSE data reduction pipeline wrapper                         #####')
    print('#####  This package is meant to be used together with ESORex and ESO MUSE pipeline   #####')
    print('#####    ftp://ftp.eso.org/pub/dfs/pipelines/muse/muse-pipeline-manual-2.4.2.pdf     #####')
    print('#####                 author: Peter Zeidler (zeidler@stsci.edu)                      #####')
    print('#####                               Dec 23, 2018                                     #####')
    print('#####                              Version: 0.4.4                                    #####')
    print('#####                                                                                #####')
    print('##########################################################################################')
    print(' ')

    ###################################################################################################
    ############################   Toggles and directories: User input   ##############################
    ###################################################################################################
    
    if configfile == None: configfile=os.path.dirname(__file__)+"/config.json"
    
    with open(configfile, "r") as read_file:
        config = json.load(read_file)
        
    withrvcorr = config['global']['withrvcorr'] #Set True if you want sky subtraction and heliocentric correction (recommended)
    OB_list=np.array(config['global']['OB_list'])#,'long_2b','long_2c']) #Name of the OB folder
    dithername=config['global']['OB'] #Name of the pointing, for which multiple OBs are used to dither
    manual_rootpath=config['global']['rootpath'] #set if you don't want to use the location of this script
    mode = config['global']['mode'] #sets the mode the OBs were taken: default: ; str: WFM-NOAO, WFM-AO, NFM-AO
    auto_sort_data = config['global']['auto_sort_data'] #Set True if the script shall sort the raw data automatically (recommended)
    using_specific_exposure_time = config['global']['using_specific_exposure_time'] #Set specific exposure time if needed, otherwise set False
    dithering_multiple_OBs = config['global']['dither_multiple_OBs'] # Set if individual dither positions are distributed via multiple OBs
    n_CPU=config['global']['n_CPU'] #Set the number of Threads used to reduce the data manually. Default and recommended is 24.
    
    using_ESO_calibration = config['calibration']['using_ESO_calibration'] #Set True if you want to use the ESO provided calibration files (recommended)
    dark = config['calibration']['dark'] #sets if the dark will be reduce (usually not necessary): default: false, bool
    renew_statics = config['calibration']['renew_statics'] #sets if the dark will be reduce (usually not necessary): default: false, bool
    
    skyreject=config['sci_basic']['skyreject'] #sets to control the sigma clipping to detect skylines in SCI_BASIC: default: 15,15,1
    
    skyfield = config['sky']['sky_field'] #sets if a skyfield is used (if available): default: auto; str: auto, object
    skyfraction = config['sky']['fraction'] # Sets the skyfraction to be used
    skyignore = config['sky']['ignore'] #Sets how much sky will be ignored
    
    skysub=config['sci_post']['subtract_sky'] #toggles whether you want skysubtraction or not in the case of the wavelength calibrated cube
    raman=config['sci_post']['raman'] #sets control if raman lines are masked.
    
    user_list = np.array(config['dither_collect']['user_list'],dtype=str)
    
    ###################################################################################################
    ########################################   END: User input   ######################################
    ###################################################################################################
    
    
    print('... Checking various necessary variables')
    assert sys.version_info > (3,0), 'YOU ARE NOT USING PYTHON 3.x'
    assert config['global']['pipeline_path'], 'NO PIPELINE PATH DEFINED'
    assert config['global']['rootpath'], 'NO ROOTPATH DEFINED'
    assert config['global']['OB'], 'NO OBs given'
    
    print('... Perfect, everything checks out')
    
        
    startime=time.time()
    print('... Settings for data reduction')
    print('>>> The data reduction will run on '+str(n_CPU)+' cores')
    
    if withrvcorr: print('>>> correct for bariocentric movement')
    else: print('>>> NO sky subtraction and NO correction for bariocentric movement')
    if skyfield == 'auto':
        print('>>> Sky runs on "automatic" and will detect if a skyfield was taken')
    if skyfield == 'object':
        print('>>> Sky runs on "OBJECT": the sky will be determined on the science field')
    print('')
    
    print('>>> The observation mode was: '+mode)
    if mode != 'NFM-AO': raman = False

    if config['calibration']['excecute'] == True:
        if using_ESO_calibration:
            print('>>> Using ESO calibration files')
            if dark: print('>>> DARK will be used')
            if not dark: print('>>> DARK will not be used')
        else:
            print('>>> Using self-processed calibration files')
            if dark: print('>>> DARK will be reduced and used')
            if not dark: print('>>> DARK will not be reduced and used')
            
    if raman == True: print('>>> Raman lines are estimated and removed (recommended for an empty field only)')
    else: print('>>> Raman lines remain in the spectrum')
    
    if dithering_multiple_OBs:
        print('>>> Exposures per pointing are spread over multiple OBs')
        print('==> The pointing name is: '+dithername)
    else:
        print('>>> All exposures are located in one OB')
        OB_list=np.array([dithername])
        if len(user_list) > 0:
            print('>>> Dithering exposures from manual input list')
            print('==> Input list: ',user_list)
    
    if dithering_multiple_OBs and len(user_list) > 0:
        print('Currently a user list of dithers cannot be provided with multiple OBs')
        sys.exit()

        
    if len(OB_list) > 1:
        print('>>> Reducing more than one OB: ',str(len(OB_list)))
    print(' ')
    print(' ')
    print('... The following modules will be excecuted')
    if config['calibration']['excecute'] == True and using_ESO_calibration == False:
        print('>>> BIAS')
        if dark: print('>>> DARK')
        print('>>> FLAT')
        print('>>> WAVECAL')
        print('>>> LSF')
        print('>>> TWILIGHT')
    if config['sci_basic']['excecute'] == True: print('>>> SCIBASIC')
    if config['std_flux']['excecute'] == True: print('>>> STANDARD')
    if config['sky']['excecute'] == True: print('>>> CREATE_SKY')
    if config['sci_post']['excecute'] == True: print('>>> SCI_POST')
    if config['exp_combine']['excecute'] == True:
        print('>>> EXP_ALIGN')
        print('>>> EXP_COMBINE')
        print(' ')
        print(' ')

    for OB in OB_list:
        print(' ')
        print('... Creating directories')
        print('>>> for OB: '+OB)
        print(' ')
        
        if manual_rootpath: rootpath=manual_rootpath
        else: rootpath=os.getcwd()+'/' #rootpath: normally set to the location of this script
        
        if dithering_multiple_OBs:
            raw_data_dir=rootpath+'raw/'+dithername+'/'+OB+'/' #path of the raw data
            working_dir=rootpath+'reduced/'+dithername+'/'+OB+'/' #path of the working directory
            combining_OBs_dir=rootpath+'reduced/'+dithername+'/'
            if not os.path.exists(combining_OBs_dir): os.mkdir(combining_OBs_dir)
        else:
            raw_data_dir=rootpath+'raw/'+OB+'/' #path of the raw data
            working_dir=rootpath+'reduced/'+OB+'/' #path of the working directory
            combining_OBs_dir = None
            
        calibration_dir=working_dir+'calibrations/' #path of the calibration file directory (in case of self prepared calibrations)
        ESO_calibration_dir=working_dir+'ESO_calibrations/' #path of the ESO calibration file directory
        static_calibration_dir=working_dir+'static_calibration_files/' #path of the static calibration file directory
        
        if renew_statics and os.path.exists(static_calibration_dir): os.rmdir(static_calibration_dir)
        if renew_statics and os.path.exists(ESO_calibration_dir): os.rmdir(ESO_calibration_dir)
        
        if not os.path.exists(rootpath+'reduced/'): os.mkdir(rootpath+'reduced/')
        if not os.path.exists(working_dir): os.mkdir(working_dir)
        if not os.path.exists(working_dir+'std/'): os.mkdir(working_dir+'std/')
        if not os.path.exists(calibration_dir): os.mkdir(calibration_dir)
        if not os.path.exists(ESO_calibration_dir): os.mkdir(ESO_calibration_dir)
        if not os.path.exists(calibration_dir+'DARK/'): os.mkdir(calibration_dir+'DARK/')
        if not os.path.exists(calibration_dir+'TWILIGHT/'): os.mkdir(calibration_dir+'TWILIGHT/')
        if not os.path.exists(calibration_dir+'SCIENCE/'): os.mkdir(calibration_dir+'SCIENCE/')
        if os.path.exists(static_calibration_dir): shutil.rmtree(static_calibration_dir)
        os.mkdir(static_calibration_dir)
        for itername in glob.glob(config['global']['pipeline_path']+'calib/muse*/cal/*.*'): shutil.copy(itername,static_calibration_dir+'.')
        
        print('... Sorting the data')
        
        if auto_sort_data:
            print('>>> The raw data will be sorted according to their header information')
            sort_data(rootpath,raw_data_dir,working_dir,ESO_calibration_dir)
        else:
            print('>>> MANUAL INTERACTION NEEDED')
    
        if using_specific_exposure_time == False:
            exposure_list=sorted(np.concatenate([glob.glob(working_dir+'*_SCI.list'),glob.glob(working_dir+'*_SKY.list')]))
            exposure_list_DARK=sorted(glob.glob(working_dir+'*_DAR.list'))
            exposure_list_TWILIGHT=sorted(glob.glob(working_dir+'*_TWI.list'))

        if using_specific_exposure_time:
            exposure_list=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_SCI.list'))
            exposure_list=sorted(np.concatenate([glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_SCI.list'),\
                                                 glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_SKY.list')]))
            exposure_list_DARK=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_DAR.list'))
            exposure_list_TWILIGHT=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_TWI.list'))
        
        for exposure_ID in range(len(exposure_list)):
            exposure_dir=exposure_list[exposure_ID][:-9]+'/'
            if os.path.exists(exposure_dir)==False: os.mkdir(exposure_dir)
            
    ###################################################################################################
    ############################   END FILE AND DIRECTORY PREPARARTION   ##############################
    ###################################################################################################
    
        print('>>>>>> reducing OB: '+OB+' <<<<<<')
        print(' ')
        ### CALIBRATION PRE-PROCESSING ###
        
        if config['calibration']['excecute'] == True:
            
            creating_sof = config['calibration']['creating_sof']
            
            if using_ESO_calibration == False:
                bias(rootpath,working_dir,exposure_list,exposure_list_DARK,exposure_list_TWILIGHT,dark,calibration_dir,creating_sof,n_CPU=n_CPU)
                if dark: dark(rootpath,working_dir,exposure_list,exposure_list_DARK,calibration_dir,creating_sof,n_CPU=n_CPU)
                flat(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,creating_sof,n_CPU=n_CPU)
                wavecal(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,n_CPU=n_CPU)
                lsf(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,n_CPU=n_CPU)
                twilight(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,dark,calibration_dir,static_calibration_dir,creating_sof,mode,n_CPU=n_CPU)
                
                
        ### OBSERVATION PRE-PROCESSING ###
        if config['sci_basic']['excecute'] == True:
            
            creating_sof = config['sci_basic']['creating_sof']
            science_pre(rootpath,working_dir,exposure_list,dark,calibration_dir,ESO_calibration_dir,static_calibration_dir,skyreject,using_ESO_calibration,creating_sof,mode,n_CPU=n_CPU)
            
        ### OBSERVATION POST-PROCESSING ###
        if config['std_flux']['excecute'] == True:
            
            creating_sof = config['std_flux']['creating_sof']
            std_flux(rootpath,working_dir,exposure_list,calibration_dir,static_calibration_dir,creating_sof,n_CPU=n_CPU)
            
        if config['sky']['excecute'] == True:
            
            creating_sof = config['sky']['creating_sof']
            if config['sky']['modified'] == False: sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,creating_sof,skyfield,skyignore,skyfraction,n_CPU=n_CPU)
            if config['sky']['modified'] == True: modified_sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,creating_sof,skyfield,skyignore,skyfraction,n_CPU=n_CPU)
        
        ###s SCIENCE POST-PROCESSING ###
        if config['sci_post']['excecute'] == True:
            
            creating_sof = config['sci_post']['creating_sof']
            scipost(rootpath,working_dir,static_calibration_dir,exposure_list,calibration_dir,ESO_calibration_dir,using_ESO_calibration,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,OB,creating_sof,raman,mode,n_CPU=n_CPU)
            
        if config['dither_collect']['excecute']: dither_collect(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,user_list,OB)
    
    if config['exp_align']['excecute']:
        creating_sof = config['exp_align']['creating_sof']
        exp_align(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,creating_sof,user_list,OB,n_CPU=24)
    if config['exp_combine']['excecute']:
        
        creating_sof = config['exp_combine']['creating_sof']
        exp_combine(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,static_calibration_dir,combining_OBs_dir,creating_sof,user_list,n_CPU=n_CPU)
    
    endtime=time.time()
    print('>>> The total execution time of the script was: ',timedelta(seconds=endtime-startime))
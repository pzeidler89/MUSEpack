'''
vers. 1.0: Excecutes the MUSE reduction pileline in the correct order
vers. 1.1: introducing only one calibration file folder per OB
vers. 1.2: choosing the illumination file closest to the observation
vers. 1.3: selecting the files for the different master file creations
vers. 1.3.2: minor corrections
vers. 1.3.3: looping order changed. each module loops by itself
vers. 1.4: always checking if calibration files already exist
vers. 1.5: Choosing between ESO calibrations and own calibrations
vers. 1.5.1: User can choose specific exposure time to reduce
vers. 2.0.0: exposures spread via different OBs for one pointing is supported.
             To do so, run the script normally for each OB inscluding scipost.
             After all are processed then run exp_align and exp_combine 
             AO observations are supported
             reduction of multiple OBs in one run supported
             set rootpath manually
             Extend the toggles
vers. 2.1.0: user can choose the number of CPUs used
vers. 2.2.0: sky subtraction can be modified, so that individual elements can be excluded
vers. 2.3.0: use json file as input
vers. 2.3.1: additional parameters added: skyreject, skysubtraction
             set parameters and modules will be shown
'''                
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
import argparse


def get_filelist(data_dir,rootpath,filename_wildcard):
    os.chdir(data_dir)
    raw_data_list=glob.glob(filename_wildcard)
    os.chdir(rootpath)
    return raw_data_list


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
                
                ### test if multiple dither positions have the same rotation angle
                
                
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
                
                ################
    for files in range(len(science_files)):
        
        hdu = fits.open(raw_data_dir+science_files[files])[0]
        RA=hdu.header['RA']
        DEC=hdu.header['DEC']
        EXPTIME=hdu.header['EXPTIME']
        ROT=hdu.header['HIERARCH ESO INS DROT POSANG']
        DATE=datetime.strptime(hdu.header['DATE-OBS'][0:10],"%Y-%m-%d")
        
        c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5').to_string('hmsdms',sep='',precision=0).translate(str.maketrans('', '', string.whitespace))
        filelist_science=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_SCIENCE.list'
        filelist_dark=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_DARK.list'
        filelist_twilight=c+'_'+str(int(EXPTIME)).rjust(4, '0')+'_'+str(int(ROT)).rjust(3, '0')+'_'+str(int(rot_angles_ident[files])).rjust(2, '0')+'_TWILIGHT.list'
        
        dark_date=np.array([])
        twilight_date=np.array([])
        
        for calibfiles in range(len(calibration_files)):
            if calibration_type[calibfiles] == 'DARK':
                hdu_temp=fits.open(raw_data_dir+calibration_files[calibfiles])[0]
                temp_date=datetime.strptime(hdu_temp.header['DATE-OBS'][0:10],"%Y-%m-%d")
                dark_date=np.append(dark_date,temp_date)
            if calibration_type[calibfiles] == 'SKYFLAT':
                hdu_temp=fits.open(raw_data_dir+calibration_files[calibfiles])[0]
                temp_date=datetime.strptime(hdu_temp.header['DATE-OBS'][0:10],"%Y-%m-%d")
                twilight_date=np.append(twilight_date,temp_date)
                
        twilight_date=np.unique(twilight_date)
        dark_date=np.unique(dark_date)
        
        f_science = open(working_dir+filelist_science, 'w')
        f_dark = open(working_dir+filelist_dark, 'w')
        f_twilight = open(working_dir+filelist_twilight, 'w')
        
        f_science.write(raw_data_dir+science_files[files]+'  '+science_type[files]+'\n')
        
        for calibfiles in range(len(calibration_files)):
            hdu_temp=fits.open(raw_data_dir+calibration_files[calibfiles])[0]
            temp_date=datetime.strptime(hdu_temp.header['DATE-OBS'][0:10],"%Y-%m-%d")
            if abs(temp_date - DATE) <= timedelta(days=1): f_science.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            if abs(temp_date - dark_date) <= timedelta(days=1): f_dark.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            if abs(temp_date - twilight_date) <= timedelta(days=1): f_twilight.write(raw_data_dir+calibration_files[calibfiles]+'  '+calibration_type[calibfiles]+'\n')
            
        f_science.close()
        f_dark.close()
        f_twilight.close()

### CALIBRATION PRE-PROCESSING ###

def bias(rootpath,working_dir,exposure_list,exposure_list_DARK,exposure_list_TWILIGHT,calibration_dir,n_CPU=24):
    
    print('... Creating the MASTER BIAS')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')
    
        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        raw_data_list_DARK = ascii.read(exposure_list_DARK[exposure_ID],format='no_header')
        raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
    
        f_science = open(calibration_dir+'SCIENCE/bias_temp.sof', 'w')
        f_dark = open(calibration_dir+'DARK/bias_temp.sof', 'w')
        f_twilight = open(calibration_dir+'TWILIGHT/bias_temp.sof', 'w')

        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'BIAS':
                f_science.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
        f_science.close()

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

            os.chdir(calibration_dir+'SCIENCE/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=bias.log --log-level=debug muse_bias --nifu=-1 --merge bias.sof')
            os.chdir(rootpath)

        if os.path.isfile(calibration_dir+'TWILIGHT/bias.sof'):
            if filecmp.cmp(calibration_dir+'TWILIGHT/bias.sof',calibration_dir+'TWILIGHT/bias_temp.sof'):
                print('TWILIGHT files are identical: proceed without error')
                os.remove(calibration_dir+'TWILIGHT/bias_temp.sof')
            else:
                sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'TWILIGHT/bias_temp.sof',calibration_dir+'TWILIGHT/bias.sof')
            
            os.chdir(calibration_dir+'TWILIGHT/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=bias.log --log-level=debug muse_bias --nifu=-1 --merge bias.sof')
            os.chdir(rootpath)
            
        if os.path.isfile(calibration_dir+'DARK/bias.sof'):
            if filecmp.cmp(calibration_dir+'DARK/bias.sof',calibration_dir+'DARK/bias_temp.sof'):
                print('dark files are identical: proceed without error')
                os.remove(calibration_dir+'DARK/bias_temp.sof')
            else:
                sys.exit('CAUTION DARK FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'DARK/bias_temp.sof',calibration_dir+'DARK/bias.sof')

            os.chdir(calibration_dir+'DARK/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=bias.log --log-level=debug muse_bias --nifu=-1 --merge bias.sof')
            os.chdir(rootpath)

def dark(rootpath,working_dir,exposure_list,exposure_list_DARK,calibration_dir,n_CPU=24):
    
    print('... Creating the MASTER DARK')
    
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
        
            os.chdir(calibration_dir+'DARK/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=dark.log --log-level=debug muse_dark --nifu=-1 --merge dark.sof')
            os.chdir(rootpath)

def flat(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,n_CPU=24):
    print('... Creating the MASTER FLAT')

    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')      
        
        
        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
        f = open(calibration_dir+'SCIENCE/flat_temp.sof', 'w')
        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'FLAT': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            #optional
            #if raw_data_list[i][1] == 'BADPIX_TABLE': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
        f.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
        #optional
        # f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
        f.close()

        if os.path.isfile(calibration_dir+'SCIENCE/flat.sof'):
            if filecmp.cmp(calibration_dir+'SCIENCE/flat.sof',calibration_dir+'SCIENCE/flat_temp.sof'):
                print('science files are identical: proceed without error')
                os.remove(calibration_dir+'SCIENCE/flat_temp.sof')
            else:
                sys.exit('CAUTION SCIENCE FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'SCIENCE/flat_temp.sof',calibration_dir+'SCIENCE/flat.sof')

            os.chdir(calibration_dir+'SCIENCE/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=flat.log --log-level=debug muse_flat --samples=true --nifu=-1 --merge flat.sof')
            os.chdir(rootpath)
        
        
        
        f = open(calibration_dir+'TWILIGHT/flat_temp.sof', 'w')
        for i in range(len(raw_data_list_TWILIGHT[1][:])):
            if raw_data_list_TWILIGHT[i][1] == 'FLAT': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
            #optional
            #if raw_data_list_TWILIGHT[i][1] == 'BADPIX_TABLE': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
        f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
        #optional
        # f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
        f.close()
    
        if os.path.isfile(calibration_dir+'TWILIGHT/flat.sof'):
            if filecmp.cmp(calibration_dir+'TWILIGHT/flat.sof',calibration_dir+'TWILIGHT/flat_temp.sof'):
                print('twilight files are identical: proceed without error')
                os.remove(calibration_dir+'TWILIGHT/flat_temp.sof')
            else:
                sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'TWILIGHT/flat_temp.sof',calibration_dir+'TWILIGHT/flat.sof')

            os.chdir(calibration_dir+'TWILIGHT/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=flat.log --log-level=debug muse_flat --samples=true --nifu=-1 --merge flat.sof')
            os.chdir(rootpath)

def wavecal(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=24):
    print('... Creating the WAVELENGTH CALIBRATION')

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
        #optional
        # f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
        # f.write(calibration_dir+'SCIENCE/MASTER_FLAT.fits MASTER_FLAT\n')
        f.close()

        if os.path.isfile(calibration_dir+'SCIENCE/wavecal.sof'):
            if filecmp.cmp(calibration_dir+'SCIENCE/wavecal.sof',calibration_dir+'SCIENCE/wavecal_temp.sof'):
                print('science files are identical: proceed without error')
                os.remove(calibration_dir+'SCIENCE/wavecal_temp.sof')
            else:
                sys.exit('CAUTION SCIENCE FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'SCIENCE/wavecal_temp.sof',calibration_dir+'SCIENCE/wavecal.sof')

            os.chdir(calibration_dir+'SCIENCE/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=wavecal.log --log-level=debug muse_wavecal --nifu=-1 --residuals --merge wavecal.sof')
            os.chdir(rootpath)

        f = open(calibration_dir+'TWILIGHT/wavecal_temp.sof', 'w')
        for i in range(len(raw_data_list_TWILIGHT[1][:])):
            if raw_data_list_TWILIGHT[i][1] == 'ARC': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
        f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
        f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
        f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
        #optional
        # f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
        # f.write(calibration_dir+'TWILIGHT/MASTER_FLAT.fits MASTER_FLAT\n')
        f.close()

        if os.path.isfile(calibration_dir+'TWILIGHT/wavecal.sof'):
            if filecmp.cmp(calibration_dir+'TWILIGHT/wavecal.sof',calibration_dir+'TWILIGHT/wavecal_temp.sof'):
                print('twilight files are identical: proceed without error')
                os.remove(calibration_dir+'TWILIGHT/wavecal_temp.sof')
            else:
                sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'TWILIGHT/wavecal_temp.sof',calibration_dir+'TWILIGHT/wavecal.sof')

            os.chdir(calibration_dir+'TWILIGHT/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=wavecal.log --log-level=debug muse_wavecal --nifu=-1 --residuals --merge wavecal.sof')
            os.chdir(rootpath)
    
def lsf(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=24):
    print('... Creating the LINE SPREAD FUNCTION')
    
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
        #optional
        # f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
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

            os.chdir(calibration_dir+'SCIENCE/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=lsf.log --log-level=debug muse_lsf --nifu=-1 --merge --save_subtracted lsf.sof')
            os.chdir(rootpath)
    
        f = open(calibration_dir+'TWILIGHT/lsf_temp.sof', 'w')
        for i in range(len(raw_data_list_TWILIGHT[1][:])):
            if raw_data_list_TWILIGHT[i][1] == 'ARC': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
        f.write(static_calibration_dir+'line_catalog.fits LINE_CATALOG\n')
        f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
        f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
        f.write(calibration_dir+'TWILIGHT/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
        #optional
        # f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
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

            os.chdir(calibration_dir+'TWILIGHT/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=lsf.log --log-level=debug muse_lsf --nifu=-1 --merge --save_subtracted lsf.sof')
            os.chdir(rootpath)

def twilight(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=24):
    print('... Creating the TWILIGHT FLAT')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')   
               
        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        raw_data_list_TWILIGHT = ascii.read(exposure_list_TWILIGHT[exposure_ID],format='no_header')
    
        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'OBJECT': MJDobject=fits.open(raw_data_list[i][0])[0].header['MJD-OBS']
    
        f = open(calibration_dir+'TWILIGHT/twilight_temp.sof', 'w')
        for i in range(len(raw_data_list_TWILIGHT[1][:])):
            if raw_data_list_TWILIGHT[i][1] == 'SKYFLAT': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')

        f.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
        f.write(calibration_dir+'TWILIGHT/MASTER_BIAS.fits MASTER_BIAS\n')
        f.write(calibration_dir+'TWILIGHT/MASTER_FLAT.fits MASTER_FLAT\n')
        # f.write(exposure_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')
        f.write(calibration_dir+'TWILIGHT/TRACE_TABLE.fits TRACE_TABLE\n')
        f.write(calibration_dir+'TWILIGHT/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
        #optional
        if MJDobject < 57823.5: f.write(static_calibration_dir+'vignetting_mask.fits VIGNETTING_MASK\n')
        if raw_data_list_TWILIGHT[i][1] == 'ILLUM': f.write(raw_data_list_TWILIGHT[i][0]+'  '+raw_data_list_TWILIGHT[i][1]+'\n')
        f.close()

        if os.path.isfile(calibration_dir+'TWILIGHT/twilight.sof'):
            if filecmp.cmp(calibration_dir+'TWILIGHT/twilight.sof',calibration_dir+'TWILIGHT/twilight_temp.sof'):
                print('twilight files are identical: proceed without error')
                os.remove(calibration_dir+'TWILIGHT/twilight_temp.sof')
            else:
                sys.exit('CAUTION TWILIGHT FILES ARE DIFFERENT: PLEASE CHECK')
        else:
            os.rename(calibration_dir+'TWILIGHT/twilight_temp.sof',calibration_dir+'TWILIGHT/twilight.sof')

            os.chdir(calibration_dir+'TWILIGHT/')
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=twilight.log --log-level=debug muse_twilight twilight.sof')
            os.chdir(rootpath)

###OBSERVATION PRE-PROCESSING###

def science_pre(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,skyreject,using_ESO_calibration,n_CPU=24):
    print('... Science PREPROCESSING')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')

        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        exposure_dir=exposure_list[exposure_ID][:-13]+'/'

        MJDsillum=np.array([])
        illum_index=np.array([])

        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'OBJECT':
                objecthdu=fits.open(raw_data_list[i][0])
                MJDobject=objecthdu[0].header['MJD-OBS']
            if raw_data_list[i][1] == 'ILLUM':
                illumhdu=fits.open(raw_data_list[i][0])
                MJDillum=illumhdu[0].header['MJD-OBS']
                MJDsillum=np.append(MJDsillum,MJDillum)
                illum_index=np.append(illum_index,i)

        choosen_illum=int(illum_index[np.argmin(np.abs(MJDsillum-MJDobject))])

        f = open(exposure_dir+'sci_pre.sof', 'w')

        for i in range(len(raw_data_list[1][:])):
            if raw_data_list[i][1] == 'OBJECT': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            if raw_data_list[i][1] == 'STD': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            # if raw_data_list[i][1] == 'SKY': f.write(raw_data_list[i][0]+'  '+raw_data_list[i][1]+'\n')
            # if raw_data_list[i][1] == 'OBJECT': f.write(raw_data_list[i][0]+'  SKY\n')

        f.write(raw_data_list[choosen_illum][0]+'  '+raw_data_list[choosen_illum][1]+'\n')
        
        if using_ESO_calibration == False:
            f.write(calibration_dir+'SCIENCE/MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(calibration_dir+'SCIENCE/MASTER_FLAT.fits MASTER_FLAT\n')
            f.write(calibration_dir+'TWILIGHT/TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
            f.write(calibration_dir+'SCIENCE/TRACE_TABLE.fits TRACE_TABLE\n')
            f.write(calibration_dir+'SCIENCE/WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            # f.write(calibration_dir+'DARK/MASTER_DARK.fits MASTER_DARK\n')

        if using_ESO_calibration == True:
            f.write(ESO_calibration_dir+'MASTER_BIAS.fits MASTER_BIAS\n')
            f.write(ESO_calibration_dir+'MASTER_FLAT.fits MASTER_FLAT\n')
            f.write(ESO_calibration_dir+'TWILIGHT_CUBE.fits TWILIGHT_CUBE\n')
            f.write(ESO_calibration_dir+'TRACE_TABLE.fits TRACE_TABLE\n')
            f.write(ESO_calibration_dir+'WAVECAL_TABLE.fits WAVECAL_TABLE\n')
            # f.write(ESO_calibration_dir+'MASTER_DARK.fits MASTER_DARK\n')
            
        f.write(static_calibration_dir+'geometry_table_wfm.fits GEOMETRY_TABLE\n')
        f.write(static_calibration_dir+'badpix_table.fits BADPIX_TABLE\n')
        #optinal
        f.close()

        os.chdir(exposure_dir)
        os.system('export OMP_NUM_THREADS='+str(n_CPU))
        os.system('esorex --log-file=sci_pre.log --log-level=debug muse_scibasic --nifu=-1 --resample --saveimage=true --skyreject='+skyreject+' --merge  sci_pre.sof')
        os.chdir(rootpath)

### OBSERVATION POST-PROCESSING ###

def std_flux(rootpath,working_dir,exposure_list,calibration_dir,static_calibration_dir,n_CPU=24):
    print('... FLUX CALIBRATION')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')      
    

        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        exposure_dir=exposure_list[exposure_ID][:-13]+'/'
    
        PIXTABLE_STD_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_STD*.fits')
    
        f = open(exposure_dir+'std_flux.sof', 'w')
        for i in range(len(PIXTABLE_STD_list)): f.write(exposure_dir+PIXTABLE_STD_list[i]+' PIXTABLE_STD\n')
        f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
        f.write(static_calibration_dir+'std_flux_table.fits STD_FLUX_TABLE\n')
        f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
        f.close()

        os.chdir(exposure_dir)
        os.system('export OMP_NUM_THREADS='+str(n_CPU))
        os.system('esorex --log-file=std_flux.log --log-level=debug muse_standard --filter=white std_flux.sof')
        os.chdir(rootpath)
    
def sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,n_CPU=24):
    print('... SKY CREATION')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')
    

        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        exposure_dir=exposure_list[exposure_ID][:-13]+'/'

        PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_OBJECT*.fits')


        f = open(exposure_dir+'sky.sof', 'w')
        for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
        if using_ESO_calibration == False:
            f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
        if using_ESO_calibration == True:
            f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
        f.write(exposure_dir+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
        f.write(exposure_dir+'STD_TELLURIC_0001.fits STD_TELLURIC\n')

        f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
        f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
        f.close()

        os.chdir(exposure_dir)
        os.system('export OMP_NUM_THREADS='+str(n_CPU))
        os.system('esorex --log-file=sky.log --log-level=debug muse_create_sky --fraction=0.005 --ignore=0.017 sky.sof')
        
        os.chdir(rootpath)
        
def modified_sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,n_CPU=24):

    print('... SKY CREATION MODIFIED')
    
    for exposure_ID in range(len(exposure_list)):
        print('>>> processing exposure: '+str(exposure_ID+1)+'/'+str(len(exposure_list)))
        print('>>> processing: '+exposure_list[exposure_ID])
        print(' ')
    
        raw_data_list = ascii.read(exposure_list[exposure_ID],format='no_header')
        exposure_dir=exposure_list[exposure_ID][:-13]+'/'

        PIXTABLE_SKY_list=get_filelist(exposure_dir,rootpath,'PIXTABLE_OBJECT*.fits')

        f = open(exposure_dir+'sky.sof', 'w')
        for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
        if using_ESO_calibration == False:
            f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
        if using_ESO_calibration == True:
            f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
        f.write(exposure_dir+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
        f.write(exposure_dir+'STD_TELLURIC_0001.fits STD_TELLURIC\n')

        f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
        f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
        f.close()

        os.chdir(exposure_dir)
        os.system('export OMP_NUM_THREADS='+str(n_CPU))
        os.system('esorex --log-file=sky.log --log-level=debug muse_create_sky --fraction=0.05 --ignore=0.05 sky.sof')
                
        os.chdir(rootpath)
        
        ### continuum set to zero
        
        os.chdir(exposure_dir)
        sky_cont_hdu=fits.open('SKY_CONTINUUM.fits',checksum=True)
        sky_cont=sky_cont_hdu[1].data
        for i in range(len(sky_cont)):
            sky_cont[i][1] = 0.
        sky_cont_hdu[1].data = sky_cont
        sky_cont_hdu.writeto('SKY_CONTINUUM_zero.fits',overwrite = True,checksum=True)
        os.chdir(rootpath)
        
        f = open(exposure_dir+'sky.sof', 'w')
        for i in range(len(PIXTABLE_SKY_list)): f.write(exposure_dir+PIXTABLE_SKY_list[i]+' PIXTABLE_SKY\n')
        if using_ESO_calibration == False:
            f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
        
        if using_ESO_calibration == True:
            f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
        
        f.write(exposure_dir+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
        f.write(exposure_dir+'STD_TELLURIC_0001.fits STD_TELLURIC\n')
        
        f.write(exposure_dir+'SKY_CONTINUUM_zero.fits SKY_CONTINUUM\n')

        f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
        f.write(static_calibration_dir+'sky_lines.fits SKY_LINES\n')
        
        f.close()

        os.chdir(exposure_dir)
        os.system('export OMP_NUM_THREADS='+str(n_CPU))
        os.system('esorex --log-file=sky.log --log-level=debug muse_create_sky --fraction=0.05 --ignore=0.05 sky.sof')
        print('SKY_CONTINUUM_zero.fits ==> SKY_CONTINUUM.fits')
        shutil.copy('SKY_CONTINUUM_zero.fits','SKY_CONTINUUM.fits')
        os.chdir(rootpath)
        
        
        os.chdir(exposure_dir)
        hdu=fits.open('SKY_LINES.fits',checksum=True)
        data=hdu[1].data

        lines_to_keep=[i for i,sky_line in enumerate(data) if (sky_line[0][:2] == 'O2' or sky_line[0][:2] == 'OH')]
        data=data[lines_to_keep]

        hdu[1].data = data
        hdu.writeto('SKY_LINES.fits', overwrite = True, checksum=True)
        os.chdir(rootpath)
        
### SCIENCE POST-PROCESSING ###

def scipost(rootpath,working_dir,static_calibration_dir,exposure_list,calibration_dir,ESO_calibration_dir,using_ESO_calibration,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,OB,n_CPU=24):
    print('... SCIENCE POST-PROCESSING')
    
    unique_pointings=np.array([])
    unique_tester=' '
    
    for expnum in range(len(exposure_list)):
        if unique_tester.find(exposure_list[expnum][:-20])==-1:
            unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-20])
            unique_tester=unique_tester+exposure_list[expnum][:-20]
    
    for unique_pointing_num in range(len(unique_pointings)):
        
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')
        
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]
        combined_exposure_dir=unique_pointings[unique_pointing_num]
        exp_list=glob.glob(combined_exposure_dir+'*SCIENCE.list')
        
        if dithering_multiple_OBs:
            if withrvcorr:
                combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
            
        if dithering_multiple_OBs == False:
            if withrvcorr:
                combining_exposure_dir_withoutsky=combined_exposure_dir+'/withoutsky_withrvcorr'
                combining_exposure_dir_withsky=combined_exposure_dir+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combined_exposure_dir+'/withsky_withoutrvcorr'
                    
        if os.path.exists(combining_exposure_dir_withoutsky)==False: os.makedirs(combining_exposure_dir_withoutsky)
        if os.path.exists(combining_exposure_dir_withsky)==False: os.makedirs(combining_exposure_dir_withsky)

        for exp_num in range(len(exp_list)):
            
            print(' ')
            print('>>> processing exposure: '+str(exp_num+1)+'/'+str(len(exp_list)))
            print(' ')
            
            PIXTABLE_OBJECT_list=get_filelist(exp_list[exp_num][:-13],rootpath,'PIXTABLE_OBJECT*.fits')
        
            f = open(exp_list[exp_num][:-13]+'/scipost.sof', 'w')
            for i in range(len(PIXTABLE_OBJECT_list)): f.write(exp_list[exp_num][:-13]+'/'+PIXTABLE_OBJECT_list[i]+' PIXTABLE_OBJECT\n')
            
            if using_ESO_calibration == False:

                f.write(calibration_dir+'SCIENCE/LSF_PROFILE.fits LSF_PROFILE\n')
                
            if using_ESO_calibration == True:

                f.write(ESO_calibration_dir+'LSF_PROFILE.fits LSF_PROFILE\n')
                
            f.write(exp_list[exp_num][:-13]+'/'+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
            f.write(exp_list[exp_num][:-13]+'/'+'STD_TELLURIC_0001.fits STD_TELLURIC\n')
            f.write(exp_list[exp_num][:-13]+'/'+'SKY_LINES.fits SKY_LINES\n')
            
            f.write(exp_list[exp_num][:-13]+'/'+'SKY_MASK.fits SKY_MASK\n')
            f.write(exp_list[exp_num][:-13]+'/'+'SKY_CONTINUUM.fits SKY_CONTINUUM\n')
            
            ## ESO STD files
            # f.write(ESO_calibration_dir+'STD_RESPONSE_0001.fits STD_RESPONSE\n')
            # f.write(ESO_calibration_dir+'STD_TELLURIC_0001.fits STD_TELLURIC\n')
            # f.write(ESO_calibration_dir+'SKY_LINES.fits SKY_LINES\n')
            
            f.write(static_calibration_dir+'astrometry_wcs_wfm.fits ASTROMETRY_WCS\n')
            f.write(static_calibration_dir+'extinct_table.fits EXTINCT_TABLE\n')
            f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
            f.close()

            os.system('export OMP_NUM_THREADS='+str(n_CPU))

            if withrvcorr:

                if skysub != 'exclude':
                    print('without sky ...')
                    os.chdir(exp_list[exp_num][:-13])
                    os.system('esorex --log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=subtract-model --filter=white scipost.sof')
                    os.chdir(rootpath)

                    print('copying files ...')
                    if dithering_multiple_OBs:
                        print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+OB+'.fits')
                        print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                    else:
                        print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir_withoutsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir_withoutsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir_withoutsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')

                if skysub != 'include':
                    print('with sky ...')

                    os.chdir(exp_list[exp_num][:-13])
                    os.system('esorex --log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=none --filter=white scipost.sof')
                    os.chdir(rootpath)

                    print('copying files ...')
                    if dithering_multiple_OBs:
                        print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+OB+'.fits')
                        print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir_withsky+'/IMAGE_FOV_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir_withsky+'/IMAGE_FOV_'+OB+'.fits')
                        print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+OB+'.fits')
                    else:
                        print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir_withsky+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir_withsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir_withsky+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                        shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir_withsky+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')

            else:

                os.chdir(exp_list[exp_num][:-13])
                os.system('esorex --log-file=scipost.log --log-level=debug muse_scipost --save=cube,skymodel,individual --skymethod=none --rvcorr=none --filter=white scipost.sof') #with sky without_rvcorr
                os.chdir(rootpath)

                print('copying files ...')
                if dithering_multiple_OBs:
                    print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir+'/DATACUBE_FINAL_'+OB+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir+'/DATACUBE_FINAL_'+OB+'.fits')
                    print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir+'/IMAGE_FOV_'+OB+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir+'/IMAGE_FOV_'+OB+'.fits')
                    print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir+'/PIXTABLE_REDUCED_'+OB+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir+'/PIXTABLE_REDUCED_'+OB+'.fits')
                else:
                    print(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits ==> '+combining_exposure_dir+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/DATACUBE_FINAL.fits',combining_exposure_dir+'/DATACUBE_FINAL_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    print(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits ==> '+combining_exposure_dir+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/IMAGE_FOV_0001.fits',combining_exposure_dir+'/IMAGE_FOV_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    print(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits ==> '+combining_exposure_dir+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
                    shutil.move(exp_list[exp_num][:-13]+'/PIXTABLE_REDUCED_0001.fits',combining_exposure_dir+'/PIXTABLE_REDUCED_'+str(exp_num+1).rjust(2, '0')+'.fits')
    
def exp_align(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,n_CPU=24):
    print('... CUBE ALIGNMENT')
    
    unique_pointings=np.array([])
    unique_tester=' '
    print(exposure_list)
    for expnum in range(len(exposure_list)):
        if unique_tester.find(exposure_list[expnum][:-20])==-1:
            unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-20])
            unique_tester=unique_tester+exposure_list[expnum][:-20]
    
    for unique_pointing_num in range(len(unique_pointings)):
        
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')
        print(unique_pointings)
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]#was -21 but didn't work for multilpe OBs ??????
        combined_exposure_dir=unique_pointings[unique_pointing_num]
        
        if dithering_multiple_OBs:
            if withrvcorr:
                print(unique_pointings_ID)
                if skysub != 'exclude': combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                if skysub != 'include': combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
            
        if dithering_multiple_OBs == False:
            if withrvcorr:
                if skysub != 'exclude': combining_exposure_dir_withoutsky=combined_exposure_dir+'/withoutsky_withrvcorr'
                if skysub != 'include': combining_exposure_dir_withsky=combined_exposure_dir+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combined_exposure_dir+'/withsky_withoutrvcorr'
        
        if withrvcorr:
            if skysub != 'exclude':
                exp_list=get_filelist(combining_exposure_dir_withoutsky,rootpath,'IMAGE_FOV_*.fits')
            
                f = open(combining_exposure_dir_withoutsky+'/exp_align.sof', 'w')
                for i in range(len(exp_list)): f.write(combining_exposure_dir_withoutsky+'/'+exp_list[i]+' IMAGE_FOV\n')
                f.close()
        
                os.chdir(combining_exposure_dir_withoutsky)
                os.system('export OMP_NUM_THREADS='+str(n_CPU))
                os.system('esorex --log-file=exp_align.log --log-level=debug muse_exp_align exp_align.sof')
                os.chdir(rootpath)
            if skysub != 'include':
                exp_list=get_filelist(combining_exposure_dir_withsky,rootpath,'IMAGE_FOV_*.fits')

                f = open(combining_exposure_dir_withsky+'/exp_align.sof', 'w')
                for i in range(len(exp_list)): f.write(combining_exposure_dir_withsky+'/'+exp_list[i]+' IMAGE_FOV\n')
                f.close()

                os.chdir(combining_exposure_dir_withsky)
                os.system('export OMP_NUM_THREADS='+str(n_CPU))
                os.system('esorex --log-file=exp_align.log --log-level=debug muse_exp_align exp_align.sof')
                os.chdir(rootpath)
            
        else:
            exp_list=get_filelist(combining_exposure_dir,rootpath,'IMAGE_FOV_*.fits')
            
            f = open(combining_exposure_dir+'/exp_align.sof', 'w')
            for i in range(len(exp_list)): f.write(combining_exposure_dir+'/'+exp_list[i]+' IMAGE_FOV\n')
            f.close()
        
            os.chdir(combining_exposure_dir)
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=exp_align.log --log-level=debug muse_exp_align exp_align.sof')
            os.chdir(rootpath)
          
def exp_combine(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,static_calibration_dir,combining_OBs_dir,n_CPU=24):
    print('... EXPOSURE COMBINATION')
    
    unique_pointings=np.array([])
    unique_tester=' '
    2,
    for expnum in range(len(exposure_list)):
        if unique_tester.find(exposure_list[expnum][:-20])==-1:
            unique_pointings=np.append(unique_pointings,exposure_list[expnum][:-20])
            unique_tester=unique_tester+exposure_list[expnum][:-20]
    
    for unique_pointing_num in range(len(unique_pointings)):
        
        print(' ')
        print('>>> processing pointing: '+str(unique_pointing_num+1)+'/'+str(len(unique_pointings)))
        print(' ')
        
        unique_pointings_ID=unique_pointings[unique_pointing_num][-18:]#was -21 but didn't work for multilpe OBs ??????
        combined_exposure_dir=unique_pointings[unique_pointing_num]
                
        if dithering_multiple_OBs:
            if withrvcorr:
                if skysub != 'exclude': combining_exposure_dir_withoutsky=combining_OBs_dir+unique_pointings_ID+'/withoutsky_withrvcorr'
                if skysub != 'inlcude': combining_exposure_dir_withsky=combining_OBs_dir+unique_pointings_ID+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combining_OBs_dir+unique_pointings_ID+'/withsky_withoutrvcorr'
            
        if dithering_multiple_OBs == False:
            if withrvcorr:
                if skysub != 'exclude': combining_exposure_dir_withoutsky=combined_exposure_dir+'/withoutsky_withrvcorr'
                if skysub != 'include': combining_exposure_dir_withsky=combined_exposure_dir+'/withsky_withrvcorr'
                
            else: combining_exposure_dir=combined_exposure_dir+'/withsky_withoutrvcorr'
        
        if withrvcorr:
            if skysub != 'exclude':
                pixtable_list=get_filelist(combining_exposure_dir_withoutsky,rootpath,'PIXTABLE_REDUCED_*.fits')
                f = open(combining_exposure_dir_withoutsky+'/exp_combine.sof', 'w')
                for i in range(len(pixtable_list)): f.write(combining_exposure_dir_withoutsky+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
                f.write(combining_exposure_dir_withoutsky+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
                f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                f.close()
                
                
                os.chdir(combining_exposure_dir_withoutsky)
                os.system('export OMP_NUM_THREADS='+str(n_CPU))
                os.system('esorex --log-file=exp_combine.log --log-level=debug muse_exp_combine --filter=white --save=cube --crsigma=5. --weight=fwhm exp_combine.sof')
                os.chdir(rootpath)
            if skysub != 'include':
                pixtable_list=get_filelist(combining_exposure_dir_withsky,rootpath,'PIXTABLE_REDUCED_*.fits')
                f = open(combining_exposure_dir_withsky+'/exp_combine.sof', 'w')
                for i in range(len(pixtable_list)): f.write(combining_exposure_dir_withsky+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
                f.write(combining_exposure_dir_withsky+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
                f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
                f.close()
                
                os.chdir(combining_exposure_dir_withsky)
                os.system('export OMP_NUM_THREADS='+str(n_CPU))
                os.system('esorex --log-file=exp_combine.log --log-level=debug muse_exp_combine --filter=white --save=cube --crsigma=5. --weight=fwhm exp_combine.sof')
                os.chdir(rootpath)
            
        else:
            pixtable_list=get_filelist(combining_exposure_dir,rootpath,'PIXTABLE_REDUCED_*.fits')
            
            f = open(combining_exposure_dir+'/exp_combine.sof', 'w')
            for i in range(len(pixtable_list)): f.write(combining_exposure_dir+'/'+pixtable_list[i]+' PIXTABLE_REDUCED\n')
            f.write(combining_exposure_dir+'/'+'OFFSET_LIST.fits OFFSET_LIST\n')
            f.write(static_calibration_dir+'filter_list.fits FILTER_LIST\n')
            f.close()
        
        
            os.chdir(combining_exposure_dir)
            os.system('export OMP_NUM_THREADS='+str(n_CPU))
            os.system('esorex --log-file=exp_combine.log --log-level=debug muse_exp_combine --filter=white --save=cube --crsigma=5. --weight=fwhm exp_combine.sof')
            os.chdir(rootpath)
    
    
# if __name__=='__main__':
def musereduce(configfile=None):
    
    print(' ')
    print('##########################################################################################')
    print('#####                                                                                #####')
    print('#####                   MUSE data reduction pipeline wrapper                         #####')
    print('#####  This package is meant to be used together with ESORex and ESO MUSE pipeline   #####')
    print('#####     ftp://ftp.eso.org/pub/dfs/pipelines/muse/muse-pipeline-manual-2.2.pdf      #####')
    print('#####                 author: Peter Zeidler (zeidler@stsci.edu)                      #####')
    print('#####                               Aug 09, 2018                                     #####')
    print('#####                              Version: 2.3.1                                    #####')
    print('#####                                                                                #####')
    print('##########################################################################################')
    print(' ')

    ###################################################################################################
    ############################   Toggles and directories: User input   ##############################
    ###################################################################################################
    
    if configfile == None:
        configfile=os.path.dirname(__file__)+"/config.json"
    else:
        configfile=configfile+'.json'
    
    with open(configfile, "r") as read_file:
        config = json.load(read_file)
        
    withrvcorr = config['global']['withrvcorr'] #Set True if you want sky subtraction and heliocentric correction (recommended)
    using_ESO_calibration = config['calibration']['using_ESO_calibration'] #Set True if you want to use the ESO provided calibration files (recommended)
    auto_sort_data = config['global']['auto_sort_data'] #Set True if the script shall sort the raw data automatically (recommended)
    using_specific_exposure_time = config['global']['using_specific_exposure_time'] #Set specific exposure time if needed, otherwise set False
    dithering_multiple_OBs = config['global']['dither_multiple_OBs'] # Set if individual dither positions are distributed via multiple OBs
    n_CPU=config['global']['n_CPU'] #Set the number of Threads used to reduce the data manually. Default and recommended is 24.

    OB_list=np.array(config['global']['OB_list'])#,'long_2b','long_2c']) #Name of the OB folder

    dithername=config['global']['OB'] #Name of the pointing, for which multiple OBs are used to dither

    manual_rootpath=config['global']['rootpath'] #set if you don't want to use the location of this script
    
    skysub=config['sci_post']['skysub'] #toggles whether you want skysubtraction or not in the case of the wavelength calibrated cube
    
    skyreject=config['sci_basic']['skyreject'] #sets to control the sigma clipping to detect skylines in SCI_BASIC: default: 15,15,1
    
    ###################################################################################################
    ########################################   END: User input   ######################################
    ###################################################################################################
    
    
    print('... Checking the python version')
    if sys.version_info < (3,0):
        print('YOU ARE NOT USING PYTHON 3.x')
        sys.exit()
    else:
        print('>>> Perfect, you are using Python 3.5 or higher')
        print(' ')
        
    startime=time.time()
    print('... Settings for data reduction')
    print('>>> The data reduction will run on '+str(n_CPU)+' cores')
    if withrvcorr: print('>>> correct for bariocentric movement')
    else: print('>>> NO sky subtraction and NO correction for bariocentric movement')

    if config['calibration']['excecute'] == True:
        if using_ESO_calibration: print('>>> Using ESO calibration files')
        else: print('>>> Using self-processed calibration files')
    
    if dithering_multiple_OBs:
        print('>>> Exposures per pointing are spread over multiple OBs')
        print('==> The pointing name is: '+dithername)
    else:
        print('>>> All exposures are located in one OB')
        OB_list=np.array([dithername])
    
    if len(OB_list) > 1:
        print('>>> Reducing more than one OB: ',str(len(OB_list)))
    
    for OB in OB_list:
        print(' ')
        print('... Creating directories')
        print('>>> for OB: '+OB)
        print(' ')
        
        if manual_rootpath:
            rootpath=manual_rootpath
        else:
            rootpath=os.getcwd()+'/' #rootpath: normally set to the location of this script
        
        if dithering_multiple_OBs:
            raw_data_dir=rootpath+'raw/'+dithername+'/'+OB+'/' #path of the raw data
            working_dir=rootpath+'reduced/'+dithername+'/'+OB+'/' #path of the working directory
            combining_OBs_dir=rootpath+'reduced/'+dithername+'/'
            if os.path.exists(combining_OBs_dir) == False: os.mkdir(combining_OBs_dir)
        else:
            raw_data_dir=rootpath+'raw/'+OB+'/' #path of the raw data
            working_dir=rootpath+'reduced/'+OB+'/' #path of the working directory
            combining_OBs_dir = None
        calibration_dir=working_dir+'calibrations/' #path of the calibration file directory (in case of self prepared calibrations)
        ESO_calibration_dir=working_dir+'ESO_calibrations/' #path of the ESO calibration file directory
        static_calibration_dir=working_dir+'static_calibration_files/' #path of the static calibration file directory
        
        if os.path.exists(rootpath+'reduced/') == False: os.mkdir(rootpath+'reduced/')
        if os.path.exists(working_dir) == False: os.mkdir(working_dir)
        if os.path.exists(calibration_dir) == False: os.mkdir(calibration_dir)
        if os.path.exists(ESO_calibration_dir) == False: os.mkdir(ESO_calibration_dir)
        if os.path.exists(calibration_dir+'DARK/') == False: os.mkdir(calibration_dir+'DARK/')
        if os.path.exists(calibration_dir+'TWILIGHT/') == False: os.mkdir(calibration_dir+'TWILIGHT/')
        if os.path.exists(calibration_dir+'SCIENCE/') == False: os.mkdir(calibration_dir+'SCIENCE/')
        if os.path.exists(static_calibration_dir) == True: shutil.rmtree(static_calibration_dir)
        os.mkdir(static_calibration_dir)
        for itername in glob.glob('/usr/local/calib/muse*/cal/*.*'): shutil.copy(itername,static_calibration_dir+'.')
        
        print('... Sorting the data')
        
        if auto_sort_data:
            print('>>> The raw data will be sorted according to their header information')
            sort_data(rootpath,raw_data_dir,working_dir,ESO_calibration_dir)
        else:
            print('>>> MANUAL INTERACTION NEEDED')
    
        if using_specific_exposure_time == False:
            exposure_list=sorted(glob.glob(working_dir+'*_SCIENCE.list'))
            exposure_list_DARK=sorted(glob.glob(working_dir+'*_DARK.list'))
            exposure_list_TWILIGHT=sorted(glob.glob(working_dir+'*_TWILIGHT.list'))

        if using_specific_exposure_time:
            exposure_list=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_SCIENCE.list'))
            exposure_list_DARK=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_DARK.list'))
            exposure_list_TWILIGHT=sorted(glob.glob(working_dir+'*'+str('{:04d}'.format(using_specific_exposure_time))+'*_TWILIGHT.list'))
        
        for exposure_ID in range(len(exposure_list)):
            exposure_dir=exposure_list[exposure_ID][:-13]+'/'
            if os.path.exists(exposure_dir)==False: os.mkdir(exposure_dir)
            
    ###################################################################################################
    ############################   END FILE AND DIRECTORY PREPARARTION   ##############################
    ###################################################################################################
    print(' ')
    print(' ')
    print('##########################################################################################')
    print('################################# Starting the pipeline ##################################')
    print('##########################################################################################')
    print(' ')
    print(' ')
    print('... The following modules will be excecuted')
    if config['calibration']['excecute'] == True and using_ESO_calibration == False:
        print('>>> BIAS')
        print('>>> DARK')
        print('>>> FLAT')
        print('>>> WAVECAL')
        print('>>> LSF')
        print('>>> TWILIGHT')
    if config['sci_basic']['excecute'] == True:
        print('>>> SCIBASIC')

    if config['std_flux']['excecute'] == True:
        print('>>> STANDARD')

    if config['sky']['excecute'] == True:
        print('>>> CREATE_SKY')

    if config['sci_post']['excecute'] == True:
        print('>>> SCI_POST')

    if config['exp_combine']['excecute'] == True:
        print('>>> EXP_ALIGN')
        print('>>> EXP_COMBINE')
        
    print(' ')
    print(' ')
    for OB in OB_list:
        
        print('>>>>>> reducing OB: '+OB+' <<<<<<')
        print(' ')
        ### CALIBRATION PRE-PROCESSING ###
        
        if config['calibration']['excecute'] == True:
            
            if using_ESO_calibration == False:
                bias(rootpath,working_dir,exposure_list,exposure_list_DARK,exposure_list_TWILIGHT,calibration_dir,n_CPU=n_CPU)
                dark(rootpath,working_dir,exposure_list,exposure_list_DARK,calibration_dir,n_CPU=n_CPU)
                flat(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,n_CPU=n_CPU)
                wavecal(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=n_CPU)
                lsf(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=n_CPU)
                twilight(rootpath,working_dir,exposure_list,exposure_list_TWILIGHT,calibration_dir,static_calibration_dir,n_CPU=n_CPU)
                
                
        ### OBSERVATION PRE-PROCESSING ###
        if config['sci_basic']['excecute'] == True:
            science_pre(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,skyreject,using_ESO_calibration,n_CPU=n_CPU)
            
        ### OBSERVATION POST-PROCESSING ###
        if config['std_flux']['excecute'] == True:
            std_flux(rootpath,working_dir,exposure_list,calibration_dir,static_calibration_dir,n_CPU=n_CPU)
            
        if config['sky']['excecute'] == True:
            if config['sky']['modified'] == False:
                sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,n_CPU=n_CPU)
            if config['sky']['modified'] == True:
                modified_sky(rootpath,working_dir,exposure_list,calibration_dir,ESO_calibration_dir,static_calibration_dir,using_ESO_calibration,n_CPU=n_CPU)
                
        ###s SCIENCE POST-PROCESSING ###
        if config['sci_post']['excecute'] == True:
            scipost(rootpath,working_dir,static_calibration_dir,exposure_list,calibration_dir,ESO_calibration_dir,using_ESO_calibration,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,OB,n_CPU=n_CPU)
    if config['exp_combine']['excecute'] == True:
        exp_align(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,combining_OBs_dir,n_CPU=n_CPU)
        exp_combine(rootpath,exposure_list,withrvcorr,skysub,dithering_multiple_OBs,static_calibration_dir,combining_OBs_dir,n_CPU=n_CPU)
    
    endtime=time.time()
    print('>>> The total execution time of the script was: ',timedelta(seconds=endtime-startime))
    

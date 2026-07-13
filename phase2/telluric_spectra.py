from astropy.io.fits import open as openFits
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import maximum_filter1d
from scipy.signal import medfilt
from matplotlib import pyplot as pl
from sys import argv
import configparser
from glob import glob
import numpy as np
import sys
import os

import warnings
warnings.filterwarnings("ignore")

# read the config file
config = configparser.RawConfigParser(allow_no_value=True)
config_filename = argv[1]
with open(config_filename,'r') as config_file:
    config.read_file(config_file)
    
# load config params
if config.get('data', 'tmpl') != 'yes':
    tmpl_ans = config.get('data', 'tmpl')
else:
    tmpl_ans = 'y'
instrument    = config.get('data', 'instrument')
target        = config.get('data', 'target')
data_path     = config.get('data', 'data_path')+target+'/'
rv_data_ans   = config.get('data', 'rv_data')
tables_path   = config.get('data', 'tables_path')
table_berv    = config.get('data', 'table_berv')
rv_avg_val  = config.getfloat('values', 'rv_avg_val')
rv_ini      = config.getfloat('values', 'rv_ini')
rv_end      = config.getfloat('values', 'rv_end')
savename      = target+'_'+instrument+'_'+config.get('data', 'savename')
opt_tmpl_ans      = config.getboolean('template', 'opt_tmpl')
skip_bjd_tmpl_ans = config.getboolean('template', 'skip_bjd_tmpl')
allow_simul       = config.getboolean('template', 'allow_simul')
rv_avg            = config.getboolean('template', 'rv_avg')
flat_omc          = config.getboolean('template', 'flat_omc')
neglect_bjd          = config.getboolean('other', 'neglect_bjd')
plot_ans_tmpl        = config.getboolean('other', 'plot_tmpl')
save_ans             = config.getboolean('other', 'save')


if instrument == 'ESPRESSO':
    #spectra_list = glob(data_path+'reflex_end_products/*/ESPRE*/*S2D_BLAZE_A*fits')
    #spectra_list = glob(data_path+'*/ES_S2SA*fits')
    spectra_list = glob(data_path+'ADP*/ES_S2DA*fits')
    
else:
    spectra_list = glob(data_path+'*_e2ds_A.fits')
    #spectra_list = glob(data_path+'/*S2D_BLAZE_A*.fits')
    #spectra_list = glob(data_path+'/*/*_S2D_A.fits')
    
spectra_list.sort()

N_spectra = len(spectra_list)
print(' * Loading',N_spectra,' spectra in',data_path,'...')
    
trasmittance_file = 'trasmittance_curve_btsettl.dat'


if len(config.get('values','special_bjd')) > 0:
    special_bjd = np.asarray(config.get('values','special_bjd').split(','),dtype='float')
else:
    special_bjd = np.array([])


if len(config.get('values','special_bjd_interval')) > 0:
    special_bjd_interval = np.asarray(config.get('values','special_bjd_interval').split(','),dtype='float')
else:
    special_bjd_interval = np.array([])


    
# Use the instrument keyword on headers
if instrument == 'HARPS':
    instrument_key = 'ESO'

if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
    instrument_key = 'ESO QC'


c = 299792.458 #light speed in vacuum km/s

if skip_bjd_tmpl_ans:
    tmpl_ans = 'y'



class Spec:
    def __init__(self, data_path, file_list):
        self.data_path  = data_path
        self.file_list = file_list


    def avoid_bjd(self, current_bjd, bjd_val, bjd_interval, neglect_bjd):
        """ Avoid spectra given a bjd interval or matching bjd values (bjd-2450000) """
        avoid_i = False

        if len(bjd_interval) > 0:
            if (round(current_bjd-2450000.,5) > bjd_interval[0]) & (round(current_bjd-2450000.,5) < bjd_interval[1]) and neglect_bjd:
                avoid_i = True

            if (round(current_bjd-2450000.,5) < bjd_interval[0]) | (round(current_bjd-2450000.,5) > bjd_interval[1]) and not neglect_bjd:
                avoid_i = True


        if len(bjd_val) > 0:
            if not neglect_bjd:
                avoid_i = True

            for bjd_val_i in bjd_val:
                if round(bjd_val_i,5) == round(current_bjd-2450000,5) and neglect_bjd:
                    avoid_i = True
                    break

                if round(bjd_val_i,5) == round(current_bjd-2450000,5) and not neglect_bjd:
                    avoid_i = False
                    break

        return avoid_i


    def load_spec(self, tmpl_ans, rv_data_ans, special_bjd, special_bjd_interval, neglect_bjd, rv_avg, rv_avg_val):
        """ Load spectra and all parameters """
        read_file1, read_file2 = True, True
        self.N_spectra =  len(self.file_list)
        self.berv_Cor = np.ones(self.N_spectra, dtype='bool')

        with openFits(self.file_list[0]) as first_fits:
            if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                self.N_order, self.N_pix = np.shape(first_fits[1].data)
            else:
                self.N_order, self.N_pix = np.shape(first_fits[0].data)


        self.Flux = np.zeros((self.N_spectra, self.N_order, self.N_pix))
        self.Wavelength = np.zeros((self.N_spectra, self.N_order, self.N_pix))
        self.Blaze = np.zeros((self.N_spectra, self.N_order, self.N_pix))
        self.errFlux = np.zeros((self.N_spectra, self.N_order, self.N_pix))

        self.Bjd, self.Berv = np.zeros(self.N_spectra), np.zeros(self.N_spectra)
        self.Rvc, self.readout_noise = np.zeros(self.N_spectra), np.zeros(self.N_spectra)
        self.SN10, self.SN50, self.SN60 = np.zeros(self.N_spectra), np.zeros(self.N_spectra), np.zeros(self.N_spectra)
        self.Airmass, self.Seeing, self.Contrast = np.zeros(self.N_spectra), np.zeros(self.N_spectra), np.zeros(self.N_spectra)
        self.FWHM, self.Bis_span = np.zeros(self.N_spectra), np.zeros(self.N_spectra)

        self.root_file = []

        self.mask_reject = np.ones(self.N_spectra, dtype=bool)

        # Creating a copy to avoid any problem with the modification of the self.file_list
        # through removing spectra
        file_list_copy = self.file_list.copy()

        #width_display = 100

        for i, file_i in enumerate(file_list_copy):
            fits_file = openFits(file_i)

            # Wavelength solutions from HARPS Fabry-Perot have *e2ds_A.fits extension (like science files)
            try:
                if instrument == 'HARPS':

                    if fits_file[0].header['OBJECT']  == 'WAVE,WAVE,FP' or fits_file[0].header['OBJECT'] == 'WAVE,WAVE,THAR2':
                        self.file_list.remove(file_i)
                        self.mask_reject[i] = False
                        continue

                if instrument == 'HARPS-N':
                    if fits_file[0].header['OBS-TYPE']  != 'SCIENCE':
                        print(file_i,'is not a SCIENCE file')
                        self.file_list.remove(file_i)
                        self.mask_reject[i] = False
                        continue


                    # The wave files have e2ds in their name so to avoid them
                    # the keyword OBJECT is used
                    kw_obj = fits_file[0].header['OBJECT']
                    if kw_obj  == 'Calibration' or kw_obj  == 'Engineering':
                        self.file_list.remove(file_i)
                        self.mask_reject[i] = False
                        continue

            except:
                print(' * Check observation type or instrument entry')
                pass

            if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                root_file_aux = file_i[:(file_i).rfind('/')]
                #root_file_aux = file_i[:(file_i).rfind(':')+7]
                current_bjd = fits_file[0].header[instrument_key+' BJD']
                rv_header_key = instrument_key+' CCF RV'
                
            elif instrument == 'HARPS' or instrument == 'HARPS-N':
                root_file_aux = file_i[:(file_i).find('_e2ds_A.fits')]
                #root_file_aux = file_i[:(file_i).rfind(':')+7]
                current_bjd = fits_file[0].header[instrument_key+' DRS BJD']
                rv_header_key = instrument_key+' DRS CCF RVC'


            if instrument == 'HARPS' or instrument == 'HARPS-REFLEX':
                # Discard simultaneous wavelength calibration with the Fabry-Perot
                # but allows simultaneous wavelength calibration with the Thorium (HARPS_ech_obs_thosimult)
                if not allow_simul and fits_file[0].header['ESO TPL NAME'] == 'HARPS_ech_obs_wavesimult':
                    self.file_list.remove(file_i)
                    self.mask_reject[i] = False
                    continue


            if len(special_bjd) > 0 or len(special_bjd_interval) > 0:
                if self.avoid_bjd(current_bjd, special_bjd,special_bjd_interval,neglect_bjd):
                    self.file_list.remove(file_i)
                    self.mask_reject[i] = False
                    continue


            if rv_data_ans.find('.rdb') != -1:
                try:
                    if read_file1:
                        rvc_file = np.loadtxt(rv_data_ans,skiprows=2,usecols=(0,1))
                        read_file1 = False

                    pos_rv = np.searchsorted(rvc_file[:,0],int(current_bjd-2400000.))
                        
                    #pos_rv = np.where(rvc_file[:,0] == round(current_bjd-2400000.,6))[0]
                    
                    #if len(pos_rv) == 0:
                    if round(rvc_file[pos_rv,0],4) != round(current_bjd-2400000.,4):
                        try:
                            if instrument != 'ESPRESSO':
                                rvc_val = openFits(glob(root_file_aux+'_ccf_*_A.fits')[0])[0].header[rv_header_key]

                            elif instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                                rvc_val = fits_file[0].header[rv_header_key]
                                    
                        except:
                            rvc_val = openFits(glob((root_file_aux+'_ccf_*_A.fits').replace(':','_',2))[0])[0].header[rv_header_key]
                            print(' * RV not found for bjd '+str(round(current_bjd-2400000.,6))+'. Using RV from CCF header')

                    else:
                        rvc_val = float(rvc_file[pos_rv,1])
                        rvc_file = np.delete(rvc_file, pos_rv, axis=0)

                except:
                    print(' * RV file not found. Using RV from mask CCF')
                    rv_data_ans = ''


            elif rv_data_ans.find('.rdb') == -1 and not rv_avg:

                try:
                    if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                        rvc_val = fits_file[0].header[rv_header_key]

                    else:
                        try:
                            rvc_val = openFits(glob(root_file_aux+'_ccf_*_A.fits')[0])[0].header[rv_header_key]

                        except:
                            rvc_val = openFits(glob((root_file_aux+'_ccf_*_A.fits').replace(':','_',2))[0])[0].header[rv_header_key]

                except:
                    print(' * CCF file not found for e2ds file',file_i)
                    self.file_list.remove(file_i)
                    self.mask_reject[i] = False
                    continue

            if instrument != 'ESPRESSO' and instrument != 'HARPS-REFLEX':

                self.Flux[i] = fits_file[0].data
                self.errFlux[i] = np.ones((self.N_order, self.N_pix)) # error Flux, to be estimated for instr !ESPRESSO


                wavelength_file_name = fits_file[0].header[instrument_key+' DRS CAL TH FILE']
                    
                if wavelength_file_name.find('e2ds') != -1:
                    wavelength_file_name = wavelength_file_name[:-12]+'_wave_A.fits'

                if os.path.exists(data_path+wavelength_file_name):
                    with openFits(data_path+wavelength_file_name) as fits_wave:
                        wave_solution = fits_wave[0].data

                elif os.path.exists((data_path+wavelength_file_name).replace(':', '_', 2)):
                    with openFits((data_path+wavelength_file_name).replace(':', '_', 2)) as fits_wave:
                        wave_solution = fits_wave[0].data

                elif os.path.exists((data_path+wavelength_file_name).replace(':', '-', 2)):
                    with openFits((data_path+wavelength_file_name).replace(':', '-', 2)) as fits_wave:
                        wave_solution = fits_wave[0].data

                else:
                    print(' * Wavelength solution file', wavelength_file_name, 'not found for spectrum', file_i[len(self.data_path):])
                    try:
                        print('Computing wavelength solution')
                        wave_solution = self.get_wave(fits_file)[0]

                    except:
                        print('Discarding spectrum', file_i[len(self.data_path):])
                        #wave_solution = self.get_wave(fits_file)[0]
                        self.file_list.remove(file_i)
                        self.mask_reject[i] = False
                        continue
 
                    
                if instrument == 'HARPS':
                    blaze_file_name = fits_file[0].header[instrument_key+' DRS BLAZE FILE']


                if instrument == 'HARPS-N':
                    blaze_file_name = fits_file[0].header[instrument_key+' DRS BLAZE FILE']

                if os.path.exists(data_path+blaze_file_name):
                    with openFits(data_path+blaze_file_name) as fits_blaze:
                        self.Blaze[i] = fits_blaze[0].data

                elif os.path.exists((data_path+blaze_file_name).replace(':', '_', 2)):
                    with openFits((data_path+blaze_file_name).replace(':', '_', 2)) as fits_blaze:
                        self.Blaze[i] = fits_blaze[0].data

                elif os.path.exists((data_path+blaze_file_name).replace(':', '-', 2)):
                    with openFits((data_path+blaze_file_name).replace(':', '-', 2)) as fits_blaze:
                        self.Blaze[i] = fits_blaze[0].data

                    
                else:
                    print(' * Missed blaze file', blaze_file_name, 'for spectrum', file_i[len(self.data_path):])
                    #self.Blaze[i] = np.ones((self.N_order, self.N_pix))############ only synth spectra test
                    self.file_list.remove(file_i)
                    self.mask_reject[i] = False
                    continue


            elif instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                blaze_file_name = ''
                
                try:
                    wavelength_file_name = glob(root_file_aux+'/*_AIR_WAVE_MATRIX*')[0]
                    wave_solution = openFits(wavelength_file_name)[1].data # Air wavelength solution from wave is NOT berv corrected
                    
                except:
                    wave_solution = fits_file[5].data # Air wavelength solution from s2d is berv corrected
                    self.berv_Cor[i] = False
                    #print('AIR_WAVE_MATRIX file not found. Air wavelength solution is BERV corrected :'+file_i)#root_file_aux)


                self.Flux[i] = fits_file[1].data # Flux without quality control
                self.errFlux[i] = fits_file[2].data # error Flux
                
                #qualdata = np.where(fits_file[3].data != 0)[0]
                #flux_aux = fits_file[1].data # Flux without quality control
                #flux_aux[qualdata] = 0
                #self.Flux[i] = flux_aux # Flux with quality control



                try:
                    blaze_file_name = glob(root_file_aux+'/Calibration_BLAZE_A*')[0]
                    self.Blaze[i] = openFits(blaze_file_name)[1].data

                except:
                    try:
                        #blaze_file_name = glob(root_file_aux+'/ToO_calibrations_BLAZE_A*')[0]
                        blaze_file_name = glob(root_file_aux+'/*cal*_BLAZE_A*')[0]
                        self.Blaze[i] = openFits(blaze_file_name)[1].data
                        
                    except:
                        #self.Blaze[i] = np.ones((self.N_order, self.N_pix))############ only synth spectra test
                        print(' * Missed blaze file',blaze_file_name,'for spectrum',file_i[len(self.data_path):])
                        self.file_list.remove(file_i)
                        self.mask_reject[i] = False
                        continue

                
            self.root_file.append(root_file_aux)
            self.Wavelength[i] = wave_solution

            if not rv_avg:
                self.Rvc[i] = rvc_val

            if instrument == 'ESPRESSO':
                if fits_file[0].header['ESO PRO REC1 CAL18 NAME'].find('SINGLEHR_2x1') != -1:
                    ccd_sigdet = 3.0

                if fits_file[0].header['ESO PRO REC1 CAL18 NAME'].find('SINGLEHR_1x1') != -1:
                    ccd_sigdet = 8.0

                ccd_conad = 1.1 # e/ADU
                RO_val = 1#np.sqrt(6.)*float(hdr[' DRS CCD SIGDET'])/float(hdr[instrument_key+' DRS CCD CONAD']) #read-out-noise

            elif instrument == 'HARPS-REFLEX':
                RO_val = np.sqrt(6.)*float(fits_file[0].header[instrument_key+' EXT0 ROX0 ROY0 BIAS RON'])

            else:
                RO_val = np.sqrt(6.)*float(fits_file[0].header[instrument_key+' DRS CCD SIGDET']) / float(fits_file[0].header[instrument_key+' DRS CCD CONAD'])

            self.Bjd[i] = current_bjd
            self.readout_noise[i] = RO_val

            if table_berv.find('rdb') == -1:

                if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                    self.Berv[i] = fits_file[0].header[instrument_key+' BERV']

                else:
                    self.Berv[i] = fits_file[0].header[instrument_key+' DRS BERV']

            else:
                if read_file2:
                    berv_file = np.loadtxt(table_berv,skiprows=2)
                    print('BERV from',table_berv)
                    read_file2 = False

                pos_berv = np.searchsorted(berv_file[:,0],int(current_bjd-2400000.))
                #pos_berv = np.where(berv_file[:,0] == round(current_bjd-24e5,6))[0]

                #if len(pos_berv) > 0:
                if round(berv_file[pos_berv,0],4) == round(current_bjd-2400000.,4):
                    berv_val = float(berv_file[pos_berv,1])
                    self.Berv[i] = berv_val
                    berv_file = np.delete(berv_file, pos_berv, axis=0)

                else:
                    print(' ### Berv not found. Using from header')
                    self.Berv[i] = fits_file[0].header[instrument_key+' DRS BERV']



            # Store the SNR around 405nm, 555nm, 612nm
            if instrument == 'HARPS':
                self.SN10[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN10']
                self.SN50[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN50']
                self.SN60[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN60']

                self.Airmass[i] = round((fits_file[0].header[instrument_key+' TEL AIRM START']+fits_file[0].header[instrument_key+' TEL AIRM END'])/2.,2)
                self.Seeing[i] = round((fits_file[0].header[instrument_key+' TEL AMBI FWHM START']+fits_file[0].header[instrument_key+' TEL AMBI FWHM END'])/2.,2)

            if instrument == 'HARPS-N':
                self.SN10[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN10']
                self.SN50[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN50']
                self.SN60[i] = fits_file[0].header[instrument_key+' DRS SPE EXT SN60']

                self.Airmass[i] = fits_file[0].header['AIRMASS']
                self.Seeing[i] = 99

            if instrument == 'HARPS-REFLEX':
                self.SN10[i] = fits_file[0].header[instrument_key+' ORDER10 SNR']
                self.SN50[i] = fits_file[0].header[instrument_key+' ORDER50 SNR']
                self.SN60[i] = fits_file[0].header[instrument_key+' ORDER60 SNR']
                
            if instrument == 'ESPRESSO':
                self.SN10[i] = fits_file[0].header[instrument_key+' ORDER21 SNR']
                self.SN50[i] = fits_file[0].header[instrument_key+' ORDER104 SNR']
                self.SN60[i] = fits_file[0].header[instrument_key+' ORDER124 SNR']              
                
                try:
                    self.Airmass[i] = round((fits_file[0].header['ESO TEL1 AIRM START']+fits_file[0].header['ESO TEL1 AIRM END'])/2.,2)
                    self.Seeing[i] = round((fits_file[0].header['ESO TEL1 AMBI FWHM START']+fits_file[0].header['ESO TEL1 AMBI FWHM END'])/2.,2)

                except:
                    try:
                        self.Airmass[i] = round((fits_file[0].header['ESO TEL2 AIRM START']+fits_file[0].header['ESO TEL2 AIRM END'])/2.,2)
                        self.Seeing[i] = round((fits_file[0].header['ESO TEL2 AMBI FWHM START']+fits_file[0].header['ESO TEL2 AMBI FWHM END'])/2.,2)
                    except:
                        try:
                            self.Airmass[i] = round((fits_file[0].header['ESO TEL3 AIRM START']+fits_file[0].header['ESO TEL3 AIRM END'])/2.,2)
                            self.Seeing[i] = round((fits_file[0].header['ESO TEL3 AMBI FWHM START']+fits_file[0].header['ESO TEL3 AMBI FWHM END'])/2.,2)
                        except:
                            try:
                                self.Airmass[i] = round((fits_file[0].header['ESO TEL4 AIRM START']+fits_file[0].header['ESO TEL4 AIRM END'])/2.,2)
                                self.Seeing[i] = round((fits_file[0].header['ESO TEL4 AMBI FWHM START']+fits_file[0].header['ESO TEL4 AMBI FWHM END'])/2.,2)
                            except:
                                print('Airmass, seeing not found for',root_file_aux)



            # Close main fits file
            fits_file.close()

            if instrument == 'ESPRESSO' or instrument == 'HARPS-REFLEX':
                self.Contrast[i] = fits_file[0].header[instrument_key+' CCF CONTRAST']
                self.FWHM[i] = fits_file[0].header[instrument_key+' CCF FWHM']
                self.Bis_span[i] = 99

            else:
                try:
                    try:
                        hdr = openFits(glob(root_file_aux+'_bis_*_A.fits')[0])[0].header
                        self.Contrast[i] = hdr[instrument_key+' DRS CCF CONTRAST']
                        self.FWHM[i] = hdr[instrument_key+' DRS CCF FWHM']
                        self.Bis_span[i] = hdr[instrument_key+' DRS BIS SPAN']
                    except:
                        hdr = openFits(glob((root_file_aux+'_bis_*_A.fits').replace(':','_',2))[0])[0].header
                        self.Contrast[i] = hdr[instrument_key+' DRS CCF CONTRAST']
                        self.FWHM[i] = hdr[instrument_key+' DRS CCF FWHM']
                        self.Bis_span[i] = hdr[instrument_key+' DRS BIS SPAN']

                except:
                    #print(' * Bis file not found for root file',self.file_list[file_i][-41:-12])
                    try:
                        hdr = openFits(glob(root_file_aux+'_ccf_*_A.fits')[0])[0].header
                        self.Contrast[i] = hdr[instrument_key+' DRS CCF CONTRAST']
                        self.FWHM[i] = hdr[instrument_key+' DRS CCF FWHM']
                        self.Bis_span[i] = 99
                    except:
                        try:
                            hdr = openFits(glob((root_file_aux+'_ccf_*_A.fits').replace(':','_',2))[0])[0].header
                            self.Contrast[i] = hdr[instrument_key+' DRS CCF CONTRAST']
                            self.FWHM[i] = hdr[instrument_key+' DRS CCF FWHM']
                            self.Bis_span[i] = 99
                        except:
                            self.Contrast[i] = 99
                            self.FWHM[i] = 99
                            self.Bis_span[i] = 99


        # The mask is applied to all the arrays to select the non rejected spectra
        self.Flux = self.Flux[self.mask_reject]
        self.Wavelength = self.Wavelength[self.mask_reject]
        self.Blaze = self.Blaze[self.mask_reject]
        self.errFlux = self.errFlux[self.mask_reject]
        self.Bjd = self.Bjd[self.mask_reject]
        self.Berv = self.Berv[self.mask_reject]
        self.berv_Cor = self.berv_Cor[self.mask_reject]
        self.Rvc = self.Rvc[self.mask_reject]
        self.readout_noise = self.readout_noise[self.mask_reject]
        self.SN10 = self.SN10[self.mask_reject]
        self.SN50 = self.SN50[self.mask_reject]
        self.SN60 = self.SN60[self.mask_reject]
        self.Airmass = self.Airmass[self.mask_reject]
        self.Seeing = self.Seeing[self.mask_reject]
        self.Contrast = self.Contrast[self.mask_reject]
        self.FWHM = self.FWHM[self.mask_reject]
        self.Bis_span = self.Bis_span[self.mask_reject]


        self.N_spectra =  len(self.file_list)  # Update in case some spec were rejected

        sys.stdout.write('\n')

        # Add the option to use a single value (average rv_ccf or manual) of RV
        if rv_avg:
            if rv_avg_val == 99:
                self.Rvc = np.zeros(self.N_spectra) + np.average(self.Rvc)
            else:
                self.Rvc = np.zeros(self.N_spectra) + rv_avg_val



    def flux(self,spec_i,order_i):
        """ Return  the flux for a given order and spectra """
        return self.Flux[spec_i][order_i].astype(np.float64)

    def errflux(self,spec_i,order_i):
        """ Return  the flux error for a given order and spectra (only ESPRESSO or HARPS-REFLEX)"""
        return self.errFlux[spec_i][order_i].astype(np.float64)

    def wave(self,spec_i,order_i):
        """ Return  the wavelength solution for a given order and spectra """
        return self.Wavelength[spec_i][order_i].astype(np.float64)

    def wave_restFrame(self,spec_i,order_i):
        """ Shift the wavelength solution to the rest frame for a given order and spectra """
        if self.berv_cor(spec_i):
            zberv = (c+self.berv(spec_i))/(c**2-self.berv(spec_i)**2)**.5
        else:
            zberv = 1.0

        zvrad = (c+self.vrad(spec_i))/(c**2-self.vrad(spec_i)**2)**.5

        return self.wave(spec_i,order_i)*zberv/zvrad

    def wave_starFrame(self,wave,spec_i):
        """ Shift the wavelength solution to the star frame for a given order and spectra """
        if self.berv_cor(spec_i):
            zberv = (c+self.berv(spec_i))/(c**2-self.berv(spec_i)**2)**.5
        else:
            zberv = 1.0

        zvrad = (c+self.vrad(spec_i))/(c**2-self.vrad(spec_i)**2)**.5
        
        return wave*zvrad/zberv


    def blaze(self,spec_i,order_i):
        """ Return the blaze function for a given order and spectra """
        return self.Blaze[spec_i][order_i].astype(np.float64)


    def subBlaze(self,spec_i,order_i,used_pix=np.array([], dtype=bool)):
        """ Subtract the blaze function for a given order and spectra and normalize """
        
        toreturn = np.zeros(self.n_pix())
        nonzeros = np.where(self.blaze(spec_i,order_i) != 0)[0]
        #nonzeros = np.where(self.blaze(spec_i,order_i) >= 0.1)[0]

        if not(False in used_pix):
            coeff = self.scale_spec(self.blaze(spec_i,order_i), self.flux(spec_i,order_i), np.arange(self.n_pix()), maxfilter=2)[1]

        else:
            coeff = self.scale_spec(self.blaze(spec_i,order_i), self.flux(spec_i,order_i), used_pix, maxfilter=2)[1]
            
        toreturn[nonzeros] = coeff*self.flux(spec_i,order_i)[nonzeros]/self.blaze(spec_i,order_i)[nonzeros]       

        return toreturn


    def addBlaze(self, noBlaze, spec_i, order_i, used_pix=np.array([], dtype=bool)):
        """ Add the blaze function and scale to a given order and spectra """
        if not(False in used_pix):
            coeff = self.scale_spec(self.flux(spec_i,order_i), self.blaze(spec_i,order_i), np.arange(self.n_pix()), maxfilter=1)[1]

        else:
            coeff = self.scale_spec(self.flux(spec_i,order_i), self.blaze(spec_i,order_i), used_pix, maxfilter=1)[1]

        return coeff*noBlaze*self.blaze(spec_i,order_i)

    def scale_spec(self, spec1, spec2, used_pix, maxfilter=1):
        """ Scale spectra 2 to spectra 1 """
        maxwidth = int(self.n_pix()*0.1)
        if maxwidth % 2 == 0:
            maxwidth += 1

        if maxfilter == 1:
            scale_coeff = np.nanmedian(maximum_filter1d(spec1[used_pix],200))/np.nanmedian(spec2[used_pix])

        elif maxfilter == 2:
            scale_coeff = np.nanmedian(spec1[used_pix])/np.nanmedian(maximum_filter1d(spec2[used_pix],200))

        elif maxfilter == 3:
            scale_coeff = np.nanmedian(spec1[used_pix])/np.nanmedian(spec2[used_pix])

        return scale_coeff*spec2, scale_coeff

    def re_norm(self, wave, flux, size, trim=0):
        """ Re-normalize a given order to ensure the flatness in the continua"""
        #####
        #pl.step(wave, flux,color='k',alpha=0.5)
        flux_med = medfilt(flux,49)
        #####
        if trim == 0:
            flux_M = maximum_filter1d(flux_med, size)
            a, b = np.polyfit(wave, flux_M, deg=1)
        else:
            flux_M = maximum_filter1d(flux_med[trim:-trim], int(self.n_pix()/5.))#size)
            a, b = np.polyfit(wave[trim:-trim], flux_M, deg=1)
            #####
            #pl.step(wave, flux_med,color='g',alpha=0.5)
            #pl.step(wave[trim:-trim], flux_M,color='r',alpha=0.5)
            #pl.step(wave, a*wave+b,color='y',alpha=0.5)
            #pl.step(wave, flux/(a*wave+b),color='b',alpha=0.5)
            #pl.show()
            #####
        return flux/(a*wave+b)

    
    def re_norm2(self, spec_i, order_i, deblazed_spec):
        #deblazed_spec = self.flux(spec_i,order_i)/self.blaze(spec_i,order_i)
    
        if order_i > 0 and order_i < self.n_order()-1:

            w1 = np.searchsorted(self.wave(spec_i,order_i),self.wave(spec_i,order_i-1)[-1])
            w2 = np.searchsorted(self.wave(spec_i,order_i),self.wave(spec_i,order_i+1)[0])

            if w1 == 0:
                w1 = self.n_pix() - w2

            if w2 == self.n_pix():
                w2 = self.n_pix() - w1

            x_edge = np.append(np.arange(w1),np.arange(w2,self.n_pix()))
            y_edge = np.append(deblazed_spec[:w1],deblazed_spec[w2:])

        
        elif order_i == 0:
        
            w2 = np.searchsorted(self.wave(spec_i,order_i),self.wave(spec_i,order_i+1)[0])
            w1 = self.n_pix() - w2

            x_edge = np.append(np.arange(w1),np.arange(w2,self.n_pix()))
            y_edge = np.append(deblazed_spec[:w1],deblazed_spec[w2:])

        
        elif order_i == self.n_order()-1:

            w1 = np.searchsorted(self.wave(spec_i,order_i),self.wave(spec_i,order_i-1)[-1])
            w2 =  self.n_pix() - w1

            x_edge = np.append(np.arange(w1),np.arange(w2,self.n_pix()))
            y_edge = np.append(deblazed_spec[:w1],deblazed_spec[w2:])
        

        a, b = np.polyfit(x_edge,y_edge,deg=1)

        return deblazed_spec/(a*np.arange(self.n_pix())+b)   

    def bjd(self,spec_i):
        """ Return the barycentric julian date for a given spectra """
        return self.Bjd[spec_i].astype(np.float64)

    def berv(self,spec_i):
        """ Return the barycentric earth radial velocity for a given spectra """
        return self.Berv[spec_i].astype(np.float64)

    def berv_cor(self,spec_i):
        """ Return if the wavelength solution is in the rest (True) or barycenter (False) frames """
        return self.berv_Cor[spec_i]

    def vrad(self,spec_i):
        """ Return the "previous (from ccf or table)" radial velocity for a given spectra """
        return self.Rvc[spec_i].astype(np.float64)

    def vrad_avrg(self):
        """ Return the average radial velocity of the data set """
        return np.average(self.Rvc)

    def argbluer(self):
        """ Return the argument of the spectra which is bluer shifted """
        return ((c+self.Rvc)/(c+self.Berv)).argmax()

    def argredder(self):
        """ Return the argument of the spectra which is reder shifted """
        return ((c+self.Rvc)/(c+self.Berv)).argmin()

    def sigma_readout(self,spec_i):
        """ Return the read-out instrumental noise """
        return self.readout_noise[spec_i].astype(np.float64)

    def sn10(self,spec_i):
        """ Return the signal to noise for order 10 """
        return self.SN10[spec_i].astype(np.float64)

    def sn50(self,spec_i):
        """ Return the signal to noise for order 50 """
        return self.SN50[spec_i].astype(np.float64)

    def sn60(self,spec_i):
        """ Return the signal to noise for order 60 """
        return self.SN60[spec_i].astype(np.float64)

    def airmass(self,spec_i):
        """ Return the average airmass """
        return self.Airmass[spec_i].astype(np.float64)

    def airmassMed(self):
        """ Return the median of the average airmass """
        return np.median(self.Airmass)

    def seeing(self,spec_i):
        """ Return the seeing """
        return self.Seeing[spec_i].astype(np.float64)

    def n_spectra(self):
        """ Return the number of spectra in the data set """
        return self.N_spectra

    def n_order(self):
        """ Return the number of orders per spectra in the data set """
        return self.N_order

    def n_pix(self):
        """ Return the number of pixels for one order (CCD)"""
        return self.N_pix

    def resamp(self,wave1,flux1,wave2):
        """ Return the cubic spline of wave1 and flux1 evaluated on wave2 """
        try:
            spline = UnivariateSpline(wave1,flux1,k=3,s=0.0,ext='zeros')
            return spline(wave2)
        except:
            return np.zeros(len(wave2))

    def rootFile(self,spec_i):
        """ Return the root of .fits file """
        return self.root_file[spec_i]

    def e2ds_list(self):
        """ Return selected e2ds file list """
        return self.file_list


def create_master(data_path, spectra_list, tmpl_ans, skip_bjd = 0, save_name = ''):

    spec = Spec(data_path,spectra_list)
    spec.load_spec(tmpl_ans,rv_data_ans,special_bjd,special_bjd_interval,neglect_bjd, rv_avg, rv_avg_val)

    rv_ccf_avg = spec.vrad_avrg()
    print(' * Number of selected spectra:',spec.n_spectra())
    print(' * '+ target+' average RV_ccf =',round(rv_ccf_avg,5),'km/s')
    print(' * Creating Template...')

    if instrument == 'ESPRESSO':
        order_share_wave = 2
    else:
        order_share_wave = 1

    if flat_omc:
        n_iterations = 5
    else:
        n_iterations = 4

        
    trim = np.zeros(int(spec.n_order()/order_share_wave))

    med_spec_noBlaze = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))
    stellar_templ_noise = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))
    telluric_templ_noise = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))
    ref_wave = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))
    master_telluric = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))

    # Final atmospheric spectrum for every selected exposure.  The historical
    # implementation wrote each order to the same FITS path inside the order
    # loop, so only the last order survived.  Keep the complete echelle matrix
    # in memory and write it once after the final telluric iteration.
    telluric_matrices = np.full(
        (spec.n_spectra(), spec.n_order(), spec.n_pix()),
        np.nan,
        dtype=np.float32,
    )

    pos_masked_orders = []
    pos_masked_orders_binary = []
    offsets_to_flat = np.ones((spec.n_spectra(),spec.n_order(),spec.n_pix()))# multiplicative offsets computed from the deviation between spectra and the template
    #offsets_to_flat = np.zeros((spec.n_spectra(),spec.n_order(),spec.n_pix()))# additive offsets computed from the deviation between spectra and the template

    compute_template_std = False
    plot_tmpl = False

    for iteration in range(n_iterations):

        for order_i in range(spec.n_order()):

            if instrument == 'ESPRESSO' and (order_i % 2) != 0:
                continue

            resamp_spec_noBlaze = np.zeros((spec.n_spectra()*order_share_wave, spec.n_pix()))
            telluric_aux = np.zeros((spec.n_spectra()*order_share_wave, spec.n_pix()))
            order_tmpl_i = int(order_i/order_share_wave)

                
            if iteration == 0:

                ref_wave[order_tmpl_i] = spec.wave_restFrame(spec.argbluer(),order_i)
                # We trim the template spectrum edges for berv+vrad shifts
                trim_val = trim_rv(spec.wave_restFrame(spec.argredder(),order_i),spec.wave_restFrame(spec.argbluer(),order_i),spec.n_pix())
                trim[order_tmpl_i] = trim_val


            for spec_i in range(spec.n_spectra()):
                
                if iteration == 0 or iteration == 4 and round(spec.bjd(spec_i),6) != round(skip_bjd,6):

                    # We resamp spectra and we evaluate them in a reference spectrum

                    for order_ii in range(order_share_wave):

                        # the initial mask is adapted to the DRS that place 0. or np.nan values in the spectra
                        mask_nan_zero = np.logical_and(np.isfinite(spec.flux(spec_i,order_i+order_ii)), spec.flux(spec_i,order_i+order_ii) != 0)


                        # the mask is shifted to account for the wavelength shift in the wavelength solution
                        mask_shift = shift_mask(mask_nan_zero, spec.wave_restFrame(spec_i,order_i+order_ii), ref_wave[order_tmpl_i])

                        noBlaze_spec = spec.subBlaze(spec_i,order_i+order_ii, used_pix=mask_nan_zero)

                            
                        if iteration == 4:
                            #noBlaze_spec = noBlaze_spec #* (offsets_to_flat[spec_i, order_i+order_ii]+1.)
                            #noBlaze_spec = noBlaze_spec - offsets_to_flat[spec_i, order_i+order_ii]
                            #noBlaze_spec = spec.re_norm(spec.wave(spec_i,order_i+order_ii), noBlaze_spec, 200, trim=40)
                            noBlaze_spec = spec.re_norm(np.arange(spec.n_pix()),noBlaze_spec,512)
                            #noBlaze_spec = spec.re_norm2(spec_i,order_i+order_ii,noBlaze_spec)


                        # to deal with potential nan in the flux, we flag them as 0. to be able to resample
                        # without problems.
                        noBlaze_spec[~mask_nan_zero] = 0.

                        # allow orders with zeros
                        resampled = np.zeros(spec.n_pix())*np.nan
                        resampled[mask_shift] = spec.resamp(spec.wave_restFrame(spec_i,order_i+order_ii), noBlaze_spec, ref_wave[order_tmpl_i])[mask_shift]
                        resamp_spec_noBlaze[order_share_wave*spec_i+order_ii] = resampled

                            
                if iteration == 1 or iteration == 3 and round(spec.bjd(spec_i),6) != round(skip_bjd,6):

                    for order_ii in range(order_share_wave):

                        # This mask will be used to prevent the bad resampling
                        # on the edges but also near the detected tellurics (for
                        # iteration 3 mostly)
                        mask_nan_zero = np.logical_and(np.isfinite(med_spec_noBlaze[order_tmpl_i]), med_spec_noBlaze[order_tmpl_i] != 0)

                        # This mask is used in the addition of the blaze because
                        # the spec.flux and the med_spec_noBlaze may not have
                        # the same nan region so if we use the mask_nan_zero
                        # in the used_pix argument, we would uncorrectly remove
                        # values in the computation of the median of the flux
                        # and the blaze when computing the scale_spec coeff
                        mask_flux = np.logical_and(np.isfinite(spec.flux(spec_i,order_i+order_ii)), spec.flux(spec_i,order_i+order_ii) != 0)

                        mask_shift = shift_mask(mask_nan_zero, spec.wave_starFrame(ref_wave[order_tmpl_i],spec_i), spec.wave(spec_i,order_i+order_ii))


                        template_blaze = np.zeros(spec.n_pix())*np.nan
                        template_blaze[mask_shift] = spec.addBlaze(spec.resamp(spec.wave_starFrame(ref_wave[order_tmpl_i], spec_i), med_spec_noBlaze[order_tmpl_i], spec.wave(spec_i, order_i+order_ii)), spec_i, order_i+order_ii, used_pix=mask_flux)[mask_shift]

                        telluric_spectrum = spec.flux(spec_i,order_i+order_ii)/template_blaze

                        if order_i ==63:
                           pl.step(spec.wave(spec_i,order_i+order_ii),telluric_spectrum)
                            #pl.show()

                        if iteration == 1 or (iteration == 3 and not flat_omc):
                            telluric_aux[order_share_wave*spec_i+order_ii] = telluric_spectrum

                        elif iteration == 3 and flat_omc:
                            valid = np.ones(spec.n_pix(), dtype='bool')
                            valid[pos_masked_orders[order_tmpl_i]] = False
                            valid = np.logical_and(valid, np.isfinite(telluric_spectrum))

                            telluric_spectrum, offset = spec.flatfunc(spec_i, order_i+order_ii, telluric_spectrum-1., operator='-', used_pix=valid)
                            telluric_aux[order_share_wave*spec_i+order_ii] = telluric_spectrum+1.

                            offsets_to_flat[spec_i, order_i+order_ii] = offset

                        if iteration == 3:
                            telluric_matrices[spec_i, order_i+order_ii] = (
                                telluric_aux[order_share_wave*spec_i+order_ii]
                            )


                if iteration == 2 and round(spec.bjd(spec_i),6) != round(skip_bjd,6):
                    # We resamp spectra and we evaluate them in a reference spectrum neglecting telluric zones

                    for order_ii in range(order_share_wave):

                        mask_nan_zero = np.logical_and(np.isfinite(spec.flux(spec_i, order_i+order_ii)), spec.flux(spec_i, order_i+order_ii) != 0)

                        mask_shift = shift_mask(mask_nan_zero, spec.wave_restFrame(spec_i,order_i+order_ii), ref_wave[order_tmpl_i])

                        # This step allows to flag the tellurics found in iteration 1 to False in the mask
                        # However, for the moment, it is not possible to merge to mask_nan_zero with the mask of tellurics
                        # and then expand and shift, this enlarge unneccessarily the tellurics zones plus they would be
                        # badly placed with the shifting.
                        if opt_tmpl_ans:
                            pos_masked_shifted = np.searchsorted(ref_wave[order_tmpl_i], spec.wave_restFrame(spec_i, order_i+order_ii)[pos_masked_orders[order_tmpl_i]])
                            # The two wavelength solutions in the searchsorted function do not have the same spacing
                            # between two wavelengths. So to keep the order in the search, sometimes the position found
                            # may be shifted from one pixel and this will create an apparent moving pixel in the telluric
                            # windows, creating noise in the computation of the median close to the tellurics.
                            # To correct this, the function fill_nan_hole will "fill" this one pixel hole with a nan value.
                            pos_masked_shifted = pos_masked_shifted[np.where((pos_masked_shifted >= 0) & (pos_masked_shifted <= spec.n_pix()-1))]
                            tmp_mask = fill_nan_hole(pos_masked_shifted, spec.N_pix)
                            mask_shift_tell = np.logical_and(mask_shift, tmp_mask)

                        else:
                            mask_shift_tell = mask_shift


                        noBlaze_spec = spec.subBlaze(spec_i, order_i+order_ii, used_pix=mask_nan_zero)

                        # This operation flags to 0 the flux where the tellurics are located
                        # (potentially flagged to nan in the spectrum).
                        # This is necessary to ensure that the resampling will work fine.
                        # Without this flagging, the zones should be np.nan which prevent the resampling to work
                            
                        noBlaze_spec[~mask_nan_zero] = 0.

                        resampled = np.zeros(spec.n_pix())*np.nan
                        resampled[mask_shift_tell] = spec.resamp(spec.wave_restFrame(spec_i,order_i+order_ii), noBlaze_spec, ref_wave[order_tmpl_i])[mask_shift_tell]
                            
                        resamp_spec_noBlaze[order_share_wave*spec_i+order_ii] = resampled


            # Out of for spec_i in range(spec.n_spectra()):

            if iteration == 0 or iteration == 2 or iteration == 4:
                # We obtain the median template spectrum for each order
                median_values_pix = np.nanmedian(resamp_spec_noBlaze, axis=0)

                # We flag all the nan values to 0
                median_values_pix[np.isnan(median_values_pix)] = 0.

                med_spec_noBlaze[order_tmpl_i] = median_values_pix

                            
                if iteration == 2 and not flat_omc:
                    compute_template_std = True

                if iteration == 4:
                    compute_template_std = True
                
                if compute_template_std:
                    # We create count_nan to account for the number of finite
                    # values used in the computation of the std to be able
                    # to divide by the square root of this number
                    count_nan = np.sum(np.isfinite(resamp_spec_noBlaze), axis=0)
                    std_values_pix = np.nanstd(resamp_spec_noBlaze, axis=0)/np.sqrt(count_nan)
                    # We flag all the nan values to 0
                    std_values_pix[np.isnan(std_values_pix)] = 0.

                    stellar_templ_noise[order_tmpl_i] = std_values_pix


            if iteration == 1 or iteration == 3:

                telluric_med = np.nanmedian(telluric_aux,axis=0)
                telluric_std = np.nanstd(telluric_aux,axis=0)/np.sqrt(spec.n_spectra())

                if not spec.berv_cor(spec_i) and False:
                    zberv = (c+spec.berv(spec_i))/(c**2-spec.berv(spec_i)**2)**.5
                    telluric_med = spec.resamp(spec.wave(spec.argbluer(),order_i)/zberv, telluric_med, ref_wave[order_tmpl_i])
                    telluric_std = spec.resamp(spec.wave(spec.argbluer(),order_i)/zberv, telluric_std, ref_wave[order_tmpl_i])

                telluric_med[np.isnan(telluric_med)] = 0.
                telluric_std[np.isnan(telluric_std)] = 0.

                if iteration == 1:
                    # Tellurics detectet at 5sigma are neglected
                    med_filt_telluric = medfilt(telluric_med,49)
                    med_telluric_med = np.nanmedian(telluric_med)
                    #std_telluric_med = np.nanmedian(telluric_std)


                    # We reflag the values corresponding to 0 in the median of the tellurics
                    # which come from initial holes or 0 in the data.
                    # Otherwise, they are flagged as tellurics and we lose 10 pixels around
                    # the holes.
                    where_zero = np.where(telluric_med == 0)[0]
                    med_filt_telluric[where_zero] = 0

                    ### old telluric detector
                    pos_masked = np.where((telluric_med > med_filt_telluric+5.*telluric_std) | (telluric_med < med_filt_telluric-5.*telluric_std))[0]

                    ### new telluric detector
                    ### in saved template 8. below correspond to *_test35*
                    #telluric_std_filt = medfilt(telluric_std,127)
                    #med_filt_telluric_std = np.std(telluric_std_filt)
                    #pos_masked = np.where((telluric_med > ((1.+3.5*med_filt_telluric_std)+telluric_std_filt)) | (telluric_med < ((1.-3.5*med_filt_telluric_std)-telluric_std_filt)))[0]
                    
                    pos_masked = add_elements(pos_masked, 5, 0, spec.n_pix())
                    pos_masked_orders.append(pos_masked)

                    binary_mask = np.isfinite(telluric_med)
                    binary_mask[pos_masked] = False

                    # Binary mask with detected tellurics (to save)
                    pos_masked_orders_binary.append(binary_mask)

                    ###
                    #pl.figure(32)
                    #pl.step(ref_wave[order_tmpl_i],telluric_med,'k')
                    #pl.step(ref_wave[order_tmpl_i],(1.+8.*med_filt_telluric_std)+telluric_std_filt,'b')
                    #pl.step(ref_wave[order_tmpl_i],(1.-8.*med_filt_telluric_std)-telluric_std_filt,'b')
                    #pl.step(ref_wave[order_tmpl_i],(1.+6.*med_filt_telluric_std)+telluric_std_filt,'g')
                    #pl.step(ref_wave[order_tmpl_i],(1.-6.*med_filt_telluric_std)-telluric_std_filt,'g')
                    #pl.step(ref_wave[order_tmpl_i],(1.+4.*med_filt_telluric_std)+telluric_std_filt,'r')
                    #pl.step(ref_wave[order_tmpl_i],(1.-4.*med_filt_telluric_std)-telluric_std_filt,'r')     
                    
                if iteration == 3:
                    telluric_templ_noise[order_tmpl_i] = telluric_std
                    master_telluric[order_tmpl_i]  = telluric_med
                    
        # Out of for order_i in range(spec.n_order()):

        #pl.show()
           
        if iteration == 1:
            med_spec_noBlaze = np.zeros((int(spec.n_order()/order_share_wave), spec.n_pix()))


        if plot_ans_tmpl:
            if iteration == 3 and not flat_omc:
                plot_tmpl = True

            if iteration == 4:
                plot_tmpl = True

        if plot_tmpl:
            for spec_i in range(spec.n_spectra()):
                for order_i in range(spec.n_order()):
                    if instrument == 'ESPRESSO':
                        order_tmpl_i = int(order_i/2.)
                    else:
                        order_tmpl_i = order_i

                    template_blaze = spec.addBlaze(spec.resamp(spec.wave_starFrame(ref_wave[order_tmpl_i],spec_i),med_spec_noBlaze[order_tmpl_i],spec.wave(spec_i,order_i)),spec_i,order_i)
                    
                    noise_blaze = spec.addBlaze(stellar_templ_noise[order_tmpl_i],spec_i,order_i)

                    try:
                        template_process_plot(spec_i,order_i,spec.bjd(spec_i),trim_val,spec.wave(spec_i,order_i),spec.flux(spec_i,order_i),template_blaze,noise_blaze,spec.sigma_readout(spec_i),master_telluric[order_tmpl_i])

                    except:
                        print(' * Error ploting the template obtaining process for order',str(order_i))
                        continue


    # Write one matrix per exposure: rows are spectral orders and columns are
    # detector pixels.  The original science header is retained for
    # provenance, with explicit output-shape metadata added below.
    from astropy.io import fits
    import matplotlib.pyplot as plt

    output_dir = os.path.join(data_path, target, 'tell_spec')
    os.makedirs(output_dir, exist_ok=True)
    preview_order = min(63, spec.n_order() - 1)

    for spec_i, science_path in enumerate(spec.file_list):
        science_file = os.path.basename(science_path)
        base_name = os.path.splitext(science_file)[0]
        output_path = os.path.join(output_dir, f"{base_name}_telluric.fits")

        with openFits(science_path) as hdul:
            header = hdul[0].header.copy()
        header['TELLURIC'] = (True, 'Atmospheric transmission spectrum')
        header['TELLMAT'] = (True, 'Data stored as order x pixel matrix')
        header['NORDERS'] = (spec.n_order(), 'Number of spectral orders')
        header['NPIX'] = (spec.n_pix(), 'Pixels per spectral order')

        fits.PrimaryHDU(
            telluric_matrices[spec_i],
            header=header,
        ).writeto(output_path, overwrite=True)

        png_path = os.path.join(output_dir, f"{base_name}_telluric.png")
        plt.figure(figsize=(8, 3))
        plt.step(
            spec.wave(spec_i, preview_order),
            telluric_matrices[spec_i, preview_order],
        )
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.ylabel('Flux Tellur')
        plt.ylim(0.85, 1.05)
        plt.tight_layout()
        plt.savefig(png_path, dpi=150)
        plt.close()

    sys.stdout.write('\n')
    if save_ans or skip_bjd != 0:
        # Save template and the edge trim values
        if save_name == '':
            save_name = savename+'_template'


        np.savez(data_path+save_name, wave=ref_wave, template=med_spec_noBlaze, template_noise = stellar_templ_noise, telluric=master_telluric, telluric_noise = telluric_templ_noise, telluric_mask=pos_masked_orders_binary, trim=trim, medAirmass = spec.airmassMed())


    if plot_ans_tmpl:
        for order_i in range(spec.n_order()):
            pl.figure(0)
            pl.step(ref_wave[order_i][trim[order_i]:-trim[order_i]],med_spec_noBlaze[order_i][trim[order_i]:-trim[order_i]],'k',linewidth=1)

        pl.xlabel(r'Wavelength $[\AA]$')
        pl.ylabel('Template')
       # pl.show()




def add_elements(elem, n_elem, elem_a, elem_b):
    """ Add n_elem consecutive positionss to each side of elements in an array. Elements are in the range [elem_a, elem_b] """
    if len(elem) == 0:
        return elem

    for i in range(-n_elem,n_elem+1,1):
        aux1, aux2 = np.array([[]]), np.array([[]])
        if i < 0:
            aux1 = np.append(aux1,[elem[1:]],1)
            aux1 = np.append(aux1,[-elem[:-1]],0)
            elem = np.append(elem,(elem-1)[np.where(np.insert(np.sum(aux1,axis=0),0,0) != 1)])
            elem.sort()

        if i > 0:
            aux2 = np.append(aux2,[-elem[:-1]],1)
            aux2 = np.append(aux2,[elem[1:]],0)
            elem = np.append(elem,(elem+1)[np.where(np.append(np.sum(aux2,axis=0),0) != 1)])
            elem.sort()

    elem = elem[np.where((elem >= elem_a) & (elem <= elem_b-1))]

    return elem


def shift_mask(initial_mask, initial_wave, final_wave):
    """ Create a new mask that follows the wavelength shift that exists between
    initial_wave and final_wave.

    """

    final_mask = np.zeros_like(initial_mask, dtype=bool)

    # To account for the wavelength shift between the two solutions, we use
    # np.searchsorted that will give the indices where the new mask values
    # should be stored.
    where_mask = np.searchsorted(initial_wave, final_wave)

    # We retrieve the edges of the wavelength solutions that are disjointed
    # to flag them to False later
    left_edge_shift = np.where(where_mask == 0)[0]
    right_edge_shift = np.where(where_mask == final_wave.size)[0]

    # To avoid any trouble with index values, we retrieve only the 'inner'
    # region of the mask
    where_mask_trim_edge = where_mask[np.logical_and(where_mask > 0, where_mask < final_wave.size)]

    # We now retrieve the good mask values in the initial_mask to set them in
    # the final_mask. As the left edge might be removed in the previous step,
    # it has to be taken into account in the range of indices for the final_mask
    final_mask[np.arange(where_mask_trim_edge.size) + left_edge_shift.size] = initial_mask[where_mask_trim_edge]
    final_mask[left_edge_shift] = False
    final_mask[right_edge_shift] = False

    # To avoid any resampling problem on the edges of False regions, we enlarge
    # them by one pixel
    where_false = np.where(final_mask == False)[0]
    # Enlarge the False regions on the right
    final_mask[where_false[:-1] + 1] = False
    # Enlarge the False regions on the left
    final_mask[where_false[1:] - 1] = False

    return final_mask


def fill_nan_hole(pos_masked_orders_i, N_pix):
    """ Given the position of the masked tellurics, fill the gaps created due to the
    distorsion between two wavelegnths solutions.
    """

    mask = np.ones(N_pix, dtype=bool)
    mask[pos_masked_orders_i] = False
    where_hole = np.where(np.diff(pos_masked_orders_i) == 2)[0]
    mask[pos_masked_orders_i[where_hole] + 1] = False

    return mask


def trim_rv(wave1, wave2, N_pix):
    """ compute the number of elements before/after the beginning/end of two wavelength solutions  """
    if np.absolute(wave1[0]-wave2[0]) < np.absolute(wave1[-1]-wave2[-1]):
        if wave1[-1] < wave2[-1]:
            max_trim = N_pix - np.searchsorted(wave2,wave1[-1])
        else:
            max_trim = N_pix - np.searchsorted(wave1,wave2[-1])

    else:
        if wave1[0] < wave2[0]:
            max_trim = np.searchsorted(wave1,wave2[0])
        else:
            max_trim = np.searchsorted(wave2,wave1[0])

    return int(max_trim)


def template_process_plot(spec_i, order_i, bjd, trim, wave, flux, template, sigma_tmpl, sigma_ro, telluric):
    #pl.ion()
    pl.figure(0,figsize=(14,8))
    pl.clf()
    pl.suptitle('Spectrum: '+str(spec_i+1)+'; BJD-2400000: '+str(round(bjd-2400000.,5))+'; Order: '+str(order_i))
    gs = GridSpec(4,5,height_ratios=[1,1,1,1],hspace=0.001,wspace=0.001)

    plot1 = pl.subplot(gs[-1,-1])
    pl.axhline(0,color='0.5',ls='--')
    omc = (flux-template)/np.sqrt(template+sigma_ro**2)

    bins_width = 3.5*np.std(omc[trim:-trim])/(len(omc[trim:-trim])**(1./3))  #Scott's normal reference rule
    n, bins, patches = pl.hist(omc[trim:-trim],bins=np.arange(-10,10,bins_width),orientation='horizontal',histtype='step',align='left',edgecolor='black',linewidth=2)
    bins = bins[:-1]
    fitGauss_param = hist_param(omc,-10,10,bins_width)
    pl.plot(Gauss(bins,*fitGauss_param),bins,'r')
    xticklabels = plot1.get_xticklabels()
    pl.setp(xticklabels, visible=False)
    pl.ylim(-10,10)

    plot2 = pl.subplot(gs[0,:-1])
    pl.plot(wave[trim:-trim],template[trim:-trim],'k')
    #pl.errorbar(wave[trim:-trim],template[trim:-trim],yerr=sigma_tmpl[trim:-trim],color='k',marker='',linestyle='-')
    ylim_1, y_lim_2 = pl.yticks()[0][0],pl.yticks()[0][-1]
    pl.ylabel('Med Templ')
    pl.axvline(5895.92424,color='k',ls='--')########
    pl.axvline(5889.95095,color='k',ls='--')########

    plot3 = pl.subplot(gs[1,:-1])
    pl.plot(wave[trim:-trim],flux[trim:-trim],'k')
    pl.ylim(ylim_1,y_lim_2)
    pl.ylabel('Spec')
    pl.axvline(5895.92424,color='k',ls='--')########
    pl.axvline(5889.95095,color='k',ls='--')########

    plot4 = pl.subplot(gs[2,:-1])
    pl.axhline(1,color='0.5',ls='--')
    pl.plot(wave[trim:-trim],telluric[trim:-trim],'k')
    pl.yticks(np.arange(-0.4, 1.5, 0.2))
    pl.ylim(0.5,1.5)
    pl.ylabel('Tellur')
    pl.axvline(5895.92424,color='k',ls='--')########
    pl.axvline(5889.95095,color='k',ls='--')########

    plot5 = pl.subplot(gs[3,:-1])
    pl.axhline(0,color='0.5',ls='--')
    pl.plot(wave[trim:-trim],omc[trim:-trim],'k')
    pl.yticks(np.arange(-6, 6, 2))
    pl.ylim(-6,6)
    pl.ylabel(r'$(O-E)/ \sigma$')
    pl.xlabel(r'Wavelength $[\AA]$')
    pl.axvline(5895.92424,color='k',ls='--')########
    pl.axvline(5889.95095,color='k',ls='--')########

    xticklabels = plot2.get_xticklabels()+plot3.get_xticklabels()+plot3.get_xticklabels()+plot4.get_xticklabels()
    pl.setp(xticklabels, visible=False)

    #pl.draw()
    #pl.pause(0.001)
   # pl.show()



if tmpl_ans == 'y' and N_spectra > 0:
    create_master(data_path,spectra_list,tmpl_ans)
    if save_ans:
        print(' * Template created in '+data_path+savename+'_template.npz')

if N_spectra == 0:
    print('No spectra were selected')

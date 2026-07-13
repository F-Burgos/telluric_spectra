
import pandas as pd
import configparser
import os
import shutil
import subprocess
from astropy.io import fits
from tqdm import tqdm
import argparse


# Argument parsing for generalization
parser = argparse.ArgumentParser(description='Process object groups for telluric spectra.')
parser.add_argument('--parquet', type=str, default='/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/phase1/metadata_final.parquet', help='Path to metadata parquet file')
parser.add_argument('--template_config', type=str, default='/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/phase2/smoke_config.ini', help='Path to template config file')
parser.add_argument('--script', type=str, default='/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/phase2/telluric_spectra.py', help='Path to telluric spectra script')
parser.add_argument('--output', type=str, default='/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/phase2/object_groups', help='Output base directory for object groups')
parser.add_argument('--calib', type=str, default='/media/nicola/4000G/HARPS/calib', help='Calibration files base directory')
parser.add_argument('--tables_path', type=str, default='spectra_results', help='Path to spectra results tables')
args = parser.parse_args()

PARQUET_PATH = args.parquet
TEMPLATE_CONFIG = args.template_config
SCRIPT_PATH = args.script
OUTPUT_BASE = args.output
CALIB_BASE = args.calib
TABLES_PATH = args.tables_path


# Read metadata table
meta = pd.read_parquet(PARQUET_PATH)

    # Add more as needed from telluric_spectra.py

# Calibration header keys to check (expand as needed)
CALIB_KEYS = [
    'DRS CAL TH FILE', 'DRS BLAZE FILE', 'DRS CCD SIGDET', 'DRS CCD CONAD',
    'DRS BERV', 'DRS SPE EXT SN10', 'DRS SPE EXT SN50', 'DRS SPE EXT SN60',
    'TEL AIRM START', 'TEL AIRM END', 'TEL AMBI FWHM START', 'TEL AMBI FWHM END',
    'ORDER10 SNR', 'ORDER50 SNR', 'ORDER60 SNR', 'ORDER21 SNR', 'ORDER104 SNR', 'ORDER124 SNR',
    # Add more as needed from telluric_spectra.py
]

MISSING_LOG = 'missing_calibrations.txt'

# Helper: search for a file first in the night folder, then anywhere under CALIB_BASE
def find_calib_file(filename, night_folder=None):
    # 1. Try night folder in calib
    if night_folder:
        night_path = os.path.join(CALIB_BASE, night_folder)
        if os.path.isdir(night_path):
            for root, dirs, files in os.walk(night_path):
                if filename in files:
                    return os.path.join(root, filename)
    # 2. Try anywhere in calib
    for root, dirs, files in os.walk(CALIB_BASE):
        if filename in files:
            return os.path.join(root, filename)
    # 3. Not found
    return None

#object_names = meta['OBJECT'].unique()

object_names = meta[meta['OBJECT'] == 'HD10700']  # Se limita a un objeto
# Progress bar for OBJECT groups

for object_name in tqdm(object_names, desc='Processing OBJECT groups'):
    
    group = meta[meta['OBJECT'] == object_name]
    object_dir = os.path.join(OUTPUT_BASE, object_name)
    os.makedirs(object_dir, exist_ok=True)
    missing_files = []
    data_paths = []
    for idx, row in group.iterrows():
        src = row['PATH']
        dst = os.path.join(object_dir, os.path.basename(src))
        # Remove existing symlink if present
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(src), dst)
        data_paths.append(dst)
        night_folder = None
        try:
            night_folder = [p for p in src.split(os.sep) if p.startswith('20')][0]
        except IndexError:
            pass
        with fits.open(src) as hdul:
            hdr = hdul[0].header
            instrument = hdr.get('INSTRUME', 'HARPS')
            instrument_key = 'ESO' if instrument == 'HARPS' else 'ESO QC'
            for key, value in hdr.items():
                if isinstance(value, str) and (value.endswith('.fits') or value.endswith('.tbl')):
                    calib_file = value
                    found = find_calib_file(calib_file, night_folder)
                    calib_dst = os.path.join(object_dir, calib_file)
                    if found:
                        if os.path.islink(calib_dst) or os.path.exists(calib_dst):
                            os.remove(calib_dst)
                        os.symlink(os.path.abspath(found), calib_dst)
                    else:
                        warning = f'Missing calibration file: {calib_file} for science file: {src}'
                        missing_files.append(warning)
                    if calib_file.endswith('e2ds_A.fits'):
                        wave_file = calib_file[:-12] + '_wave_A.fits'
                        found_wave = find_calib_file(wave_file, night_folder)
                        wave_dst = os.path.join(object_dir, wave_file)
                        if found_wave:
                            if os.path.islink(wave_dst) or os.path.exists(wave_dst):
                                os.remove(wave_dst)
                            os.symlink(os.path.abspath(found_wave), wave_dst)
                        else:
                            warning = f'Missing calibration file: {wave_file} for science file: {src}'
                            missing_files.append(warning)
            for k in CALIB_KEYS:
                for prefix in ['', instrument_key + ' ']:
                    full_key = prefix + k
                    if full_key in hdr:
                        value = hdr[full_key]
                        if isinstance(value, str) and (value.endswith('.fits') or value.endswith('.tbl')):
                            calib_file = value
                            found = find_calib_file(calib_file, night_folder)
                            calib_dst = os.path.join(object_dir, calib_file)
                            if found:
                                if os.path.islink(calib_dst) or os.path.exists(calib_dst):
                                    os.remove(calib_dst)
                                os.symlink(os.path.abspath(found), calib_dst)
                            else:
                                warning = f'Missing calibration file: {calib_file} for science file: {src}'
                                missing_files.append(warning)
    total_calib_files = len(set([line.split(':')[1].split('for')[0].strip() for line in missing_files])) if missing_files else 0
    if missing_files:
        with open(os.path.join(object_dir, MISSING_LOG), 'w') as f:
            for line in missing_files:
                f.write(line + '\n')
    # Prepare config for this OBJECT
    data_path = os.path.abspath(OUTPUT_BASE) + '/'
    config = configparser.RawConfigParser(allow_no_value=True)
    config.read(TEMPLATE_CONFIG)
    config.set('data', 'target', object_name)
    config.set('data', 'data_path', data_path)
    config.set('data', 'tables_path', os.path.abspath(TABLES_PATH))
    config_out = os.path.join(object_dir, f'{object_name}.ini')
    with open(config_out, 'w') as f:
        config.write(f)
    telluric_cmd = ['python3', SCRIPT_PATH, config_out]
    print(f'Running: {" ".join(telluric_cmd)}')
    subprocess.run(telluric_cmd, check=True)
    print(f'Done with OBJECT {object_name}.')

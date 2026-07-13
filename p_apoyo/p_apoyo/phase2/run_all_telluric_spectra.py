"""Deprecated launcher.

Use ``run_telluric_pipeline.sh`` at the project root.  The maintained Stage 2
driver is ``smoke_object_group_telluric_v2.py``; this file is retained only as
historical reference.
"""

import os
import shutil
import configparser
import pandas as pd
from collections import defaultdict
import subprocess

# Paths
TABLE_PATH = '../telluric/tabla1_science_HARPS_50arcsec_groups.txt'
TEMPLATE_CONFIG = '../telluric/naira_config_HARPS_1.ini'
OUTPUT_DIR = './spectra_results'
CONFIGS_DIR = './configs'
SCRIPT_PATH = './telluric_spectra.py'

# Ensure output and configs directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CONFIGS_DIR, exist_ok=True)

# Read table
# The table is tab-separated, header in first row

df = pd.read_csv(TABLE_PATH, sep='\t')

# Group FITS files by OBJECT
object_groups = defaultdict(list)
for _, row in df.iterrows():
    object_groups[row['OBJECT']].append(row['FILENAME'])

for obj, fits_files in object_groups.items():
    # Prepare config for this object
    config = configparser.RawConfigParser(allow_no_value=True)
    config.read(TEMPLATE_CONFIG)
    config.set('data', 'target', obj)
    # Assume data_path is the directory containing the FITS files (take from first file)
    fits_dir = os.path.dirname(fits_files[0])
    config.set('data', 'data_path', fits_dir + '/')
    # Save config
    config_path = os.path.join(CONFIGS_DIR, f'{obj}_config.ini')
    with open(config_path, 'w') as f:
        config.write(f)
    # Run telluric_spectra.py
    print(f'Processing {obj}...')
    subprocess.run(['python3', SCRIPT_PATH, config_path], check=True)
    # Find output .npz file
    instrument = config.get('data', 'instrument')
    savename = config.get('data', 'savename')
    npz_name = f'{obj}_{instrument}_{savename}_template.npz'
    npz_path = os.path.join(fits_dir, npz_name)
    if os.path.exists(npz_path):
        shutil.move(npz_path, os.path.join(OUTPUT_DIR, npz_name))
    else:
        print(f'Warning: Output file {npz_path} not found for {obj}')

print('All done.')

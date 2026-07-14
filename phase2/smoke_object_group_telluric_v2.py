
import pandas as pd
import configparser
import os
import shutil
import subprocess
from astropy.io import fits
from tqdm import tqdm
import argparse
import json
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

# Argument parsing for generalization
parser = argparse.ArgumentParser(description='Process object groups for telluric spectra.')
parser.add_argument('--parquet', type=str, default=os.path.join(PROJECT_DIR, 'phase1', 'metadata_final.parquet'), help='Path to metadata parquet file')
parser.add_argument('--template_config', type=str, default=os.path.join(SCRIPT_DIR, 'smoke_config.ini'), help='Path to template config file')
parser.add_argument('--script', type=str, default=os.path.join(SCRIPT_DIR, 'telluric_spectra.py'), help='Path to telluric spectra script')
parser.add_argument('--output', type=str, default='spectra', help='Output base directory for object groups')
parser.add_argument('--calib', type=str, default='/mnt/disco_datos/data/HARPS/calib', help='Calibration files base directory')
parser.add_argument('--science_root', type=str, default=None, help='Current science root, used to resolve portable REL_PATH values')
parser.add_argument('--tables_path', type=str, default='spectra_results', help='Path to spectra results tables')
parser.add_argument('--star', type=str, default='all', help='Name of the star to process (default: all)')
parser.add_argument('--python', type=str, default='python3', help='Python executable used to run telluric_spectra.py')
parser.add_argument(
    '--plot_orders',
    type=str,
    default='',
    help='Optional telluric preview orders to plot, e.g. "63", "10,20", "50-55", or "all". Default: no plots.',
)
args = parser.parse_args()

PARQUET_PATH = args.parquet
TEMPLATE_CONFIG = args.template_config
SCRIPT_PATH = args.script
OUTPUT_BASE = args.output
CALIB_BASE = args.calib
TABLES_PATH = args.tables_path
PLOT_ORDERS = args.plot_orders.strip()


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
FAILED_LOG = 'telluric_run_failed.txt'


def resolve_science_path(row):
    """Resolve an observation from its stored path or the configured science root."""
    path = str(row['PATH'])
    candidates = [path]

    if args.science_root:
        rel_path = row.get('REL_PATH')
        if pd.notna(rel_path):
            candidates.append(os.path.join(args.science_root, str(rel_path)))
        candidates.append(
            os.path.join(args.science_root, str(row.get('NIGHT', '')), os.path.basename(path))
        )

    legacy_map = {
        '/media/nicola/T7ext4_data/HARPS/science': '/mnt/disco_datos/data/HARPS/science',
        '/media/nicola/4000G/HARPS/science': '/mnt/disco_datos/data/HARPS/science',
    }

    for old_prefix, new_prefix in legacy_map.items():
        if path.startswith(old_prefix):
            candidates.append(path.replace(old_prefix, new_prefix, 1))

    for candidate in candidates:
        if os.path.isfile(candidate) and os.access(candidate, os.R_OK):
            return candidate

    # Return the original path if no readable candidate is found.
    return path

# Index calibration files once. Re-walking a multi-year calibration tree for
# every science exposure made the old implementation scale as O(N_science *
# N_calibration).
calib_by_name = defaultdict(list)
calib_by_night_name = {}
print(f'Indexing calibration files under {CALIB_BASE}...')
for root, dirs, files in os.walk(CALIB_BASE):
    rel_root = os.path.relpath(root, CALIB_BASE)
    night = rel_root.split(os.sep, 1)[0]
    for filename in files:
        full_path = os.path.join(root, filename)
        calib_by_name[filename].append(full_path)
        calib_by_night_name.setdefault((night, filename), full_path)
print(f'Indexed {sum(len(paths) for paths in calib_by_name.values())} calibration files.')


# Helper: search for a file first in the night folder, then anywhere under CALIB_BASE
def find_calib_file(filename, night_folder=None):
    if night_folder:
        night_match = calib_by_night_name.get((night_folder, filename))
        if night_match:
            return night_match
    matches = calib_by_name.get(filename, [])
    return matches[0] if matches else None


def safe_name(value):
    """Return the same filesystem-safe object name used by telluric_spectra.py."""
    return ''.join(
        char if char.isalnum() or char in ('-', '_', '.') else '_'
        for char in str(value)
    )


# Progress bar for OBJECT groups
if args.star == 'all':
    object_names = meta['OBJECT'].unique()
else:
    object_names = [args.star]

    
for object_name in tqdm(object_names, desc='Processing OBJECT groups',
                                     unit="obj",
                                     total=len(object_names)):

    object_name = str(object_name)
    group = meta[meta['OBJECT'].astype(str) == object_name]

    object_dir = os.path.join(OUTPUT_BASE, object_name)
    os.makedirs(object_dir, exist_ok=True)
    missing_files = []
    data_paths = []
    staged_e2ds = []
    
    #for idx, row in group.iterrows():
    for idx, row in tqdm(group.iterrows(), desc='Processing files',
                                     unit="files",
                                     total=len(group)):
        
        src = resolve_science_path(row)
        dst = os.path.join(object_dir, os.path.basename(src))

        if not os.path.isfile(src) or not os.access(src, os.R_OK):
            warning = f'Unreadable science file (skipped): {src}'
            print(f'WARNING: {warning}')
            missing_files.append(warning)
            continue

        # Remove existing symlink if present
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)

        os.symlink(os.path.abspath(src), dst)
        if dst.endswith('_e2ds_A.fits'):
            staged_e2ds.append(dst)
            
        data_paths.append(dst)
        night_folder = str(row.get('NIGHT', ''))

        if not src.endswith('e2ds_A.fits'):
            continue
            
        try:
            with fits.open(src) as hdul:
                hdr = hdul[0].header
        except Exception as exc:
            warning = f'Cannot open FITS file (skipped): {src} | {exc}'
            print(f'WARNING: {warning}')
            missing_files.append(warning)
            continue

        instrument = hdr.get('INSTRUME', 'HARPS')
        instrument_key = 'ESO' if instrument == 'HARPS' else 'ESO QC'
        
        blaze_file = hdr.get(instrument_key + ' DRS BLAZE FILE')
        if not blaze_file:
            warning = f'Missing blaze key in header for science file: {src}'
            print(f'WARNING: {warning}')
            missing_files.append(warning)
            continue

        found_blaze = find_calib_file(blaze_file, night_folder)
        blaze_dst = os.path.join(object_dir, blaze_file)
        if found_blaze:
            if os.path.islink(blaze_dst) or os.path.exists(blaze_dst):
                os.remove(blaze_dst)
                
            os.symlink(os.path.abspath(found_blaze), blaze_dst)
        else:
            warning = f'Missing blaze file: {blaze_file} for science file: {src}'
            missing_files.append(warning)

        
        wave_file = hdr.get(instrument_key + ' DRS CAL TH FILE')
        if not wave_file:
            warning = f'Missing wavelength key in header for science file: {src}'
            print(f'WARNING: {warning}')
            missing_files.append(warning)
            continue

        if wave_file.find('e2ds_A.fits') != -1:
            wave_file = wave_file.replace('e2ds','wave')

        found_wave = find_calib_file(wave_file, night_folder)
        wave_dst = os.path.join(object_dir, wave_file)
        
        if found_wave:
            if os.path.islink(wave_dst) or os.path.exists(wave_dst):
                os.remove(wave_dst)
                
            os.symlink(os.path.abspath(found_wave), wave_dst)
        else:
            warning = f'Missing wave file: {wave_file} for science file: {src}'
            missing_files.append(warning)
                

        # Por ahora no se hace nada con los valores buscados en header
        #with fits.open(src) as hdul:
        #    hdr = hdul[0].header
        #    instrument = hdr.get('INSTRUME', 'HARPS')
        #    instrument_key = 'ESO' if instrument == 'HARPS' else 'ESO QC'

        #    for k in CALIB_KEYS:
        #        for prefix in ['', instrument_key + ' ']:
        #            full_key = prefix + k
        #            if full_key in hdr:
        #                value = hdr[full_key]


    
    total_calib_files = len(set([line.split(':')[1].split('for')[0].strip() for line in missing_files])) if missing_files else 0
    if missing_files:
        with open(os.path.join(object_dir, MISSING_LOG), 'w') as f:
            for line in missing_files:
                f.write(line + '\n')
    elif os.path.exists(os.path.join(object_dir, MISSING_LOG)):
        os.remove(os.path.join(object_dir, MISSING_LOG))
    # Prepare config for this OBJECT
    data_path = os.path.abspath(OUTPUT_BASE) + '/'
    config = configparser.RawConfigParser(allow_no_value=True)
    config.read(TEMPLATE_CONFIG)
    config.set('data', 'target', object_name)
    config.set('data', 'data_path', data_path)
    config.set('data', 'tables_path', os.path.abspath(TABLES_PATH))
    if not config.has_section('output'):
        config.add_section('output')
    config.set('output', 'plot_telluric', 'yes' if PLOT_ORDERS else 'no')
    config.set('output', 'plot_orders', PLOT_ORDERS)
    config_out = os.path.join(object_dir, f'{object_name}.ini')
    with open(config_out, 'w') as f:
        config.write(f)
    telluric_cmd = [args.python, SCRIPT_PATH, config_out]
    print(f'Running: {" ".join(telluric_cmd)}')
    try:
        subprocess.run(telluric_cmd, check=True)
    except subprocess.CalledProcessError as exc:
        warning = (
            f'telluric_spectra.py failed for OBJECT {object_name} '
            f'with exit code {exc.returncode}'
        )
        print(f'WARNING: {warning}')
        with open(os.path.join(object_dir, FAILED_LOG), 'w') as f:
            f.write(warning + '\n')
        continue
    if os.path.exists(os.path.join(object_dir, FAILED_LOG)):
        os.remove(os.path.join(object_dir, FAILED_LOG))

    tell_spec_dir = os.path.join(object_dir, object_name, 'tell_spec')
    cube_path = os.path.join(tell_spec_dir, f'{safe_name(object_name)}_telluric_cube.fits')
    missing_outputs = []
    invalid_outputs = []
    expected_cube_shape = None
    cube_shape = None
    if not staged_e2ds:
        missing_outputs.append('No staged e2ds_A spectra for object')
    elif not os.path.isfile(cube_path):
        missing_outputs.append(cube_path)
    else:
        try:
            with fits.open(staged_e2ds[0], memmap=True) as science_hdul:
                expected_cube_shape = (len(staged_e2ds),) + tuple(science_hdul[0].data.shape)
            with fits.open(cube_path, memmap=True) as output_hdul:
                cube_shape = tuple(output_hdul[0].data.shape)
                header = output_hdul[0].header
                is_cube = bool(header.get('TELLCUBE', False))
            if cube_shape != expected_cube_shape or not is_cube:
                invalid_outputs.append(
                    {
                        'path': cube_path,
                        'expected_shape': expected_cube_shape,
                        'shape': cube_shape,
                        'tellcube_header': is_cube,
                    }
                )
        except Exception as exc:
            invalid_outputs.append({'path': cube_path, 'error': str(exc)})

    summary = {
        'object': object_name,
        'input_e2ds': len(staged_e2ds),
        'telluric_cube': cube_path,
        'telluric_outputs': 1 if staged_e2ds and not missing_outputs and not invalid_outputs else 0,
        'expected_cube_shape': expected_cube_shape,
        'cube_shape': cube_shape,
        'plot_orders': PLOT_ORDERS,
        'missing_outputs': missing_outputs,
        'invalid_outputs': invalid_outputs,
        'missing_calibration_messages': len(missing_files),
    }
    with open(os.path.join(object_dir, 'processing_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    if missing_outputs or invalid_outputs:
        print(
            f'WARNING: {object_name} did not produce a valid all-order telluric cube '
            f'for {summary["input_e2ds"]} input spectra.'
        )
    print(f'Done with OBJECT {object_name}.')

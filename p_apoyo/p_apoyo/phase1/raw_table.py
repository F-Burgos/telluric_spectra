# -*- coding: utf-8 -*-
# =============================================
# Script: raw_table.py
# Description: Extracts and compiles raw metadata from FITS files in science directories.
# Inputs: FITS files from science directories (ROOT_DIR)
# Outputs: metadata_raw.parquet (compiled raw metadata table)
# =============================================
"""
Created on Thu Dec 11 14:27:19 2025

@author: pipeb
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from glob import glob

import numpy as np
import pandas as pd
from astropy.io import fits
from tqdm import tqdm

# =========================

# =========================
# 1. Configuración de rutas y parámetros por argumentos
# =========================
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract and compile raw metadata from FITS files.")
    parser.add_argument('--root_dir', type=str, default="/media/nicola/4000G/HARPS/science/", help='Root directory containing science FITS files')
    parser.add_argument('--output_dir', type=str, default="/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/", help='Directory to save the output table')
    parser.add_argument('--output_name', type=str, default="metadata_raw.parquet", help='Name of the output parquet file')
    parser.add_argument(
        '--products',
        type=str,
        default='e2ds_A,ccf_A,s1d_A,bis_A',
        help=(
            'Comma-separated HARPS products to inspect. Supported values: '
            'e2ds_A, ccf_A, s1d_A, bis_A, all. Filtering occurs before FITS files are opened.'
        ),
    )
    return parser.parse_args()

args = parse_args()
ROOT_DIR = args.root_dir
OUTPUT_DIR = args.output_dir
OUTPUT_NAME = args.output_name

# =========================
# 2. Keywords a extraer del header
# =========================
# Puedes ajustar esta lista según lo que te interese
KEYWORDS = [
    "RA",
    "DEC",
    "OBJECT",
    "MJD-OBS",
    "ESO DRS BJD",
    "ESO DRS BERV",
    "ESO DRS SPE EXT SN50",
    "EXPTIME",
    "ESO INS DET1 TMMEAN",
    "ESO TEL AMBI FWHM START",
    "ESO TEL AMBI FWHM END",
    "ESO TEL AIRM START",
    "ESO TEL AIRM END",
    "ESO TEL AMBI PRES START",
    "ESO TEL AMBI PRES END",
    "ESO TEL AMBI RHUM",
    "ESO TEL AMBI TEMP",
    "ESO TEL AMBI WINDDIR",
    "ESO TEL AMBI WINDSP",
    "ESO TEL MOON DEC",
    "ESO TEL MOON RA",
    "ESO DRS VERSION",
    "ESO DPR TYPE",
]

# Normalizamos los nombres de columnas para pandas
BASE_COLUMNS = [
    k.replace(" ", "_").replace("-", "_").replace(".", "_") for k in KEYWORDS
]

# Añadimos columnas adicionales que queremos siempre
EXTRA_COLUMNS = ["NIGHT", "FILENAME", "PATH", "REL_PATH", "PRODUCT_TYPE", "OBJ_ID"]

ALL_COLUMNS = BASE_COLUMNS + EXTRA_COLUMNS

# =========================
# 3. Localizar todos los archivos FITS
# =========================
# Patrón: subcarpeta (noche) / archivo .fits
# Si quieres restringir a cierto tipo de producto, por ejemplo "*e2ds_A.fits", cambia el patrón.
def product_type(filename):
    """Return the science-product role encoded in a HARPS filename."""
    if filename.endswith('_e2ds_A.fits'):
        return 'e2ds_A'
    if filename.endswith('_s1d_A.fits'):
        return 's1d_A'
    if '_ccf_' in filename and filename.endswith('_A.fits'):
        return 'ccf_A'
    if '_bis_' in filename and filename.endswith('_A.fits'):
        return 'bis_A'
    return 'other'


requested_products = {item.strip() for item in args.products.split(',') if item.strip()}
pattern = os.path.join(ROOT_DIR, "*", "*.fits")
candidate_files = sorted(glob(pattern))
if 'all' in requested_products:
    fits_files = candidate_files
else:
    fits_files = [
        path for path in candidate_files
        if product_type(os.path.basename(path)) in requested_products
    ]

print(
    f"Encontrados {len(candidate_files)} FITS bajo {ROOT_DIR}; "
    f"se leerán {len(fits_files)} productos seleccionados"
)

# =========================
# 4. Función para obtener un ID de objeto
# =========================
def get_obj_id(header, object_str):
    """
Asigna OBJ_ID, este sera el nombre original existente en el header.
    """
    return str(object_str).strip() if object_str is not None else np.nan

# =========================
# 5. Recorrer archivos y extraer metadatos
# =========================
rows = []

for path in tqdm(fits_files, desc="Leyendo FITS", unit="archivo"):
    try:
        with fits.open(path) as hdul:
            header = hdul[0].header

            # OBJECT más seguro
            obj = (header.get("OBJECT", "") or "").strip()
            row = {}

            # 5.1 Extraemos todos los keywords
            for key, col_name in zip(KEYWORDS, BASE_COLUMNS):
                val = header.get(key, np.nan)
                # Para OBJECT quitamos espacios internos
                if key == "OBJECT" and val is not None:
                    val = str(val).replace(" ", "")
                row[col_name] = val

            # 5.2 Información de la estructura de directorios
            night = os.path.basename(os.path.dirname(path))
            filename = os.path.basename(path)
            abs_path = os.path.abspath(path)

            row["NIGHT"] = night
            row["FILENAME"] = filename
            row["PATH"] = abs_path
            row["REL_PATH"] = os.path.relpath(abs_path, os.path.abspath(ROOT_DIR))
            row["PRODUCT_TYPE"] = product_type(filename)

            # 5.3 ID de objeto observado
            row["OBJ_ID"] = get_obj_id(header, obj)

            rows.append(row)
            # Fin de lectura FITS AGREGAR CLOSE

    except Exception as e:
        # Si algo falla con un FITS, lo reportamos y seguimos
        print(f"Error leyendo {path}: {e}")

# =========================
# 6. Construir DataFrame metadata_raw
# =========================
metadata_raw = pd.DataFrame(rows, columns=ALL_COLUMNS)

print(f"metadata_raw construido con {len(metadata_raw)} filas.")

# Opcional: conversión de RA/DEC a numérico
for col in ["RA", "DEC"]:
    if col in metadata_raw.columns:
        metadata_raw[col] = pd.to_numeric(metadata_raw[col], errors="coerce")

# =========================
# 7. Guardar a disco
# =========================
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)

metadata_raw.to_parquet(output_path, index=False)
print(f"Tabla metadata_raw guardada en: {output_path}")

# Si prefieres CSV:
# metadata_raw.to_csv(os.path.join(OUTPUT_DIR, "metadata_raw.csv"), index=False)

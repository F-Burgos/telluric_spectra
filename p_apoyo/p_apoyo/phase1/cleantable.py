# -*- coding: utf-8 -*-
# =============================================
# Script: cleantable.py
# Description: Cleans and clusters astronomical metadata.
# Inputs: metadata_raw.parquet (raw metadata table)
# Outputs: metadata.parquet (cleaned table), metadata_final.parquet (final clustered table), metadata_remaining.parquet (optional, ungrouped rows)
# =============================================
"""
Created on Thu Dec 11 17:06:27 2025

@author: pipeb
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from collections import deque

# =========================

# =========================
# Configuración de rutas y parámetros por argumentos
# =========================
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Clean and cluster astronomical metadata.")
    parser.add_argument('--base_dir', type=str, default="/media/nicola/T7ext4/fondecyt24/telluric/p-apoyo/phase1/", help='Base directory for input/output files')
    parser.add_argument('--raw_path', type=str, default=None, help='Path to raw metadata parquet file')
    parser.add_argument('--clean_path', type=str, default=None, help='Path to cleaned metadata parquet file')
    parser.add_argument('--final_path', type=str, default=None, help='Path to final clustered metadata parquet file')
    parser.add_argument('--remaining_path', type=str, default=None, help='Path to remaining (ungrouped) metadata parquet file')
    parser.add_argument('--radius_near', type=float, default=50.0, help='Radius (arcsec) to consider objects as neighbors')
    return parser.parse_args()

args = parse_args()
BASE_DIR = args.base_dir
RAW_PATH = args.raw_path if args.raw_path else os.path.join(BASE_DIR, "metadata_raw.parquet")
CLEAN_PATH = args.clean_path if args.clean_path else os.path.join(BASE_DIR, "metadata.parquet")
FINAL_PATH = args.final_path if args.final_path else os.path.join(BASE_DIR, "metadata_final.parquet")
REMAINING_PATH = args.remaining_path if args.remaining_path else os.path.join(BASE_DIR, "metadata_remaining.parquet")
RADIUS_NEAR = args.radius_near


# =========================
# Paso 0: limpiar si es necesario
# =========================
def ensure_clean_metadata():
    """
    Si 'metadata.parquet' no existe:
      - lee metadata_raw.parquet
      - limpia RA/DEC nulos
      - calcula OBS_COUNT (radio definido en --radius_near [arcsec])
      - guarda metadata.parquet

    Si ya existe, no hace nada.
    """
    if os.path.exists(CLEAN_PATH):
        print(f"Encontrado archivo limpio: {CLEAN_PATH}. Se omite la fase de limpieza.")
        return

    if not os.path.exists(RAW_PATH):
        raise FileNotFoundError(
            f"No se encontró {CLEAN_PATH} ni {RAW_PATH}. "
            f"Primero debe ejecutarse el script que genera metadata_raw.parquet."
        )

    print(f"No existe {CLEAN_PATH}. Leyendo tabla cruda {RAW_PATH} para limpiar...")

    # Read in chunks if possible (for very large files)
    # For Parquet, chunked reading is not natively supported, but for CSV use pd.read_csv(..., chunksize=...)
    metadata = pd.read_parquet(RAW_PATH)

    # Optimize dtypes for memory
    dtype_map = {
        "RA": np.float32,
        "DEC": np.float32,
        "OBJ_ID": "category",
        "OBJECT": "category"
    }
    for col, dtype in dtype_map.items():
        if col in metadata.columns:
            metadata[col] = metadata[col].astype(dtype)


    # solo selecciona la fibra A
    metadata = metadata[metadata['PATH'].str.endswith('A.fits')]
    
    n_original = len(metadata)
    print(f"Filas originales en metadata_raw: {n_original}")

    # Asegurar RA/DEC numéricos
    for col in ["RA", "DEC"]:
        if col not in metadata.columns:
            raise ValueError(f"La columna '{col}' no existe en metadata_raw.")
        metadata[col] = pd.to_numeric(metadata[col], errors="coerce")

    # Eliminar filas con RA o DEC nulos
    mask_valid = metadata["RA"].notna() & metadata["DEC"].notna()
    n_valid = mask_valid.sum()
    n_eliminadas = n_original - n_valid

    print(f"Filas con RA/DEC no válidos eliminadas: {n_eliminadas}")
    print(f"Filas restantes tras limpieza: {n_valid}")

    metadata = metadata.loc[mask_valid].reset_index(drop=True)

    # Opcional: ordenar por RA
    metadata = metadata.sort_values(by="RA").reset_index(drop=True)

    # Calcular OBS_COUNT en un radio de XX arcsec (como en el diseño original)
    #print("Calculando OBS_COUNT (objetos dentro de XX arcsec)...")
    # For very large datasets, consider using BallTree from sklearn for neighbor search
    # from sklearn.neighbors import BallTree
    # coords_rad = np.deg2rad(np.c_[metadata["RA"], metadata["DEC"]])
    # tree = BallTree(coords_rad, metric='haversine')
    # radius = 50.0 / 3600.0 * np.pi / 180.0  # arcsec to radians
    # ind = tree.query_radius(coords_rad, r=radius)
    # obs_count = np.array([len(i) for i in ind])
    # metadata["OBS_COUNT"] = obs_count

    coords = SkyCoord(
    ra=metadata["RA"].values * u.deg,
    dec=metadata["DEC"].values * u.deg,
    frame="fk5",
    )

    coords_roll = SkyCoord(
    ra=np.roll(metadata["RA"].values,shift=1) * u.deg,
    dec=np.roll(metadata["DEC"].values,shift=1) * u.deg,
    frame="fk5",
    )

    obs_count = np.zeros(len(metadata), dtype=int)
    sep = coords.separation(coords_roll)

    other_star = np.asarray(sep.to(u.arcsec)) < RADIUS_NEAR


    i = 0
    while i < len(other_star):
        # Inicio de un nuevo bloque: contar cuántas entradas consecutivas son True
        j = i
        sameroot = 0
        while j < len(other_star) and other_star[j]:
            
            if not metadata.iloc[j]["PATH"].endswith('_e2ds_A.fits'):
                sameroot += 1
                
            j += 1
            
        block_size = j - i - sameroot  # número de entradas True en este bloque

        # Asignar el total del bloque a cada posición correspondiente
        obs_count[i-1:j] = block_size

        # Avanzar al siguiente elemento (el False que cierra el bloque, si existe)
        i = j + 1  # saltar el False (cambio de estrella)


    metadata["OBS_COUNT"] = obs_count 

    print("OBS_COUNT calculado.")

    # Save intermediate result to disk to free memory
    metadata.to_parquet(CLEAN_PATH, index=False)
    print(f"Tabla limpia guardada en: {CLEAN_PATH}")


# =========================
# Paso 2: agrupar por OBJECT + vecinos cercanos
# =========================
def build_final_metadata():
    print(f"Leyendo tabla limpia desde: {CLEAN_PATH}")
    metadata = pd.read_parquet(CLEAN_PATH).reset_index(drop=True)

    for col in ["RA", "DEC", "OBJECT", "OBJ_ID"]:
        if col not in metadata.columns:
            raise ValueError(f"La tabla limpia debe contener la columna '{col}'.")

    n_total = len(metadata)
    print(f"Filas en metadata (limpia): {n_total}")

    # Coordenadas
    coords = SkyCoord(
        ra=metadata["RA"].values * u.deg,
        dec=metadata["DEC"].values * u.deg,
        frame="fk5",
    )

    print(f"Buscando vecinos dentro de {RADIUS_NEAR} arcsec usando separación directa...")

    # Coordenadas desplazadas (equivalente al Paso 1)
    coords_roll = SkyCoord(
        ra=np.roll(metadata["RA"].values, shift=1) * u.deg,
        dec=np.roll(metadata["DEC"].values, shift=1) * u.deg,
        frame="fk5",
    )

    # Separación angular entre elementos consecutivos
    sep = coords.separation(coords_roll)

    # Condición de vecindad (en grados)
    is_neighbor = sep < (RADIUS_NEAR / 3600.0) * u.deg

    # Construcción de lista de vecinos
    n = len(metadata)
    neighbors = [set() for _ in range(n)]

    print("Construyendo listas de vecinos (adyacencia secuencial)...")

    for i in tqdm(range(n), desc="Vecinos", unit="obj"):
        if is_neighbor[i]:
            j = i - 1  # vecino anterior (por construcción del roll)
            if j >= 0:
                neighbors[i].add(j)
                neighbors[j].add(i)


    print("Construyendo clusters por conectividad geométrica (sin sesgo por OBJECT)...")

    assigned = np.zeros(n, dtype=bool)
    metadata_final_list = []

    for i in tqdm(range(n), desc="Clusters", unit="obj"):
        if assigned[i]:
            continue

        # BFS: construir componente conexa completa
        queue = deque([i])
        cluster_indices = set([i])
        assigned[i] = True

        while queue:
            current = queue.popleft()

            for neighbor in neighbors[current]:
                if not assigned[neighbor]:
                    assigned[neighbor] = True
                    cluster_indices.add(neighbor)
                    queue.append(neighbor)

        cluster_indices = sorted(cluster_indices)
        cluster_df = metadata.iloc[cluster_indices].copy()

        # -------------------------
        # ASIGNACIÓN CONSISTENTE DE OBJECT
        # -------------------------
        value_counts = cluster_df["OBJECT"].value_counts(dropna=True)

        if not value_counts.empty:
            majority_object = value_counts.idxmax()
            cluster_df["OBJECT"] = majority_object

            majority_count = value_counts.max()
            total_count = value_counts.sum()
            ambiguous = majority_count <= 0.5 * total_count
            cluster_df["ambiguous"] = ambiguous
        else:
            cluster_df["ambiguous"] = True

        metadata_final_list.append(cluster_df)

    
    # Concatenar clusters
    if metadata_final_list:
        metadata_final = pd.concat(metadata_final_list, ignore_index=True)
        # Remove duplicates based on PATH column
        if "PATH" in metadata_final.columns:
            before = len(metadata_final)
            metadata_final = metadata_final.drop_duplicates(subset=["PATH"])
            after = len(metadata_final)
            print(f"Filas eliminadas por duplicados en PATH: {before - after}")
    else:
        metadata_final = pd.DataFrame(columns=metadata.columns)

    # Print the number of clusters formed
    print(f"Número de clusters formados: {len(metadata_final_list)}")

    n_final = len(metadata_final)
    print(f"Filas en metadata_final: {n_final}")

    # Filas no agrupadas
    remaining = metadata.loc[~assigned].reset_index(drop=True)
    n_remaining = len(remaining)
    print(f"Filas sin agrupar: {n_remaining}")

    # Guardar
    metadata_final.to_parquet(FINAL_PATH, index=False)
    metadata.to_parquet(CLEAN_PATH, index=False)
    remaining.to_parquet(REMAINING_PATH, index=False)

    print(f"metadata_final guardado en: {FINAL_PATH}")
    print(f"metadata.parquet (actualizado) guardado en: {CLEAN_PATH}")
    print(f"metadata_remaining guardado en: {REMAINING_PATH}")
    
    
# =========================
# Main
# =========================
if __name__ == "__main__":
    # Paso 0: asegurar que metadata.parquet exista y esté limpio
    ensure_clean_metadata()

    # Paso 2: construir metadata_final a partir de metadata.parquet
    build_final_metadata()


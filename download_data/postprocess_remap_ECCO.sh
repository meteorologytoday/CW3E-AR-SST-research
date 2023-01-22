#!/bin/bash

remap_filename="remap_ECCO_to_ERA5_bilinear.nc"
raw_dir="data/ECCO/processed_0p50deg"
remapped_dir="data/ECCO/processed_0p25deg"

if [ ! -f "$remap_filename" ]; then
    echo "Remapping file $remap_filename does not exist. Generate it now..."
    python generate_remapping_file_ECCO_to_ERA5.py
fi

echo "Start remapping"
ncremap -m $remap_filename -I $raw_dir -O $remapped_dir

#!/bin/bash

remap_filename="NEMO_to_ECMWF_bilinear.nc"
raw_dir="data/ORA5/processed"
remapped_dir="data/ORA5/processed_remapped"

if [ ! -f "$remap_filename" ]; then
    echo "Remapping file $remap_filename does not exist. Generate it now..."
    python generate_remapping_file_ORA5_to_ERA5.py
fi

echo "Start remapping"
ncremap -m $remap_filename -I $raw_dir -O $remapped_dir

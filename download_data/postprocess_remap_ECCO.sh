#!/bin/bash

remap_filename="remap_ECCO_to_ERA5_bilinear.nc"

for MLD_algo in "ECCO-MLD" "RHO-03-MLD"; do

    raw_dir="data/ECCO/processed_0p50deg_${MLD_algo}"
    remapped_dir="data/ECCO/processed_0p25deg_${MLD_algo}"

    if [ ! -f "$remap_filename" ]; then
        echo "Remapping file $remap_filename does not exist. Generate it now..."
        python generate_remapping_file_ECCO_to_ERA5.py
    fi

    if [ ! -d $remapped_dir ] ; then

        echo "Start remapping"
        ncremap -m $remap_filename -I $raw_dir -O $remapped_dir

    else
        
        echo "Directory $remapped_dir already exists. Skip it."

    fi
done

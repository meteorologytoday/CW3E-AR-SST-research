#!/bin/bash

for method in "ANOM_LEN" "TOTIVT250" ; do 
    python3 make_ERA5_AR_objects.py  --method $method --AR-clim-dir data/ERA5/AR_processed_clim --nproc 10 &
done


wait

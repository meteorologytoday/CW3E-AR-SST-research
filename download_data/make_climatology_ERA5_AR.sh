#!/bin/bash

beg_year=2015
end_year=2017

python3 make_climatology.py \
    --beg-year $beg_year \
    --end-year $end_year \
    --input-dir data/ERA5/AR_processed \
    --output-dir data/ERA5/AR_processed_clim \
    --filename-prefix ERA5_AR_ \
    --nproc 3

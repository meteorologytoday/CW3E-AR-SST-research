#!/bin/bash


output_dir=figures
mkdir -p $output_dir

python3 plot_AR_freq.py \
    --input-dir output_ECCO/1993-2017_10N-60N-n25_120E-120W-n60 \
    --output $output_dir/AR_freq1.png \
    --beg-year 1993 \
    --end-year 2017  \
    --freq-max 0.25 \
    --ndays 5 \
    --threshold-days 4 &

python3 plot_AR_freq.py \
    --input-dir output_ECCO/1993-2017_10N-60N-n25_120E-120W-n60 \
    --output $output_dir/AR_freq2.png \
    --beg-year 1993 \
    --end-year 2017 \
    --freq-max 0.15  \
    --ndays 10 \
    --threshold-days 8 &

wait


#python3 plot_G_terms_map.py --input-dir output_ECCO/1993-2017_10N-60N-n25_120E-120W-n60 --output $output_dir/AR_forcing_partition1.png --no-display

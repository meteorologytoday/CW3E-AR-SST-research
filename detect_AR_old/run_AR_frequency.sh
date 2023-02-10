#!/bin/bash

. pretty_latlon.sh


#output_AR=AR_frequency_statistics.nc
#output_nonAR=ARnon_frequency_statistics.nc



beg_year=1998
end_year=2017

output_dir="output_AR_stat"
output_dir_AR_vs_nonAR=$output_dir/AR_vs_nonAR
output_dir_AR_vs_clim=$output_dir/AR_vs_clim


mkdir -p $output_dir
mkdir -p $output_dir_AR_vs_clim
mkdir -p $output_dir_AR_vs_nonAR

spatial_rngs=(
    30 65 -160 -110
)

yrng_text=$( printf "%04d-%04d" $beg_year $end_year )
output_AR=$output_dir/AR_frequency_statistics_${yrng_text}_AR.nc
output_nonAR=$output_dir/AR_frequency_statistics_${yrng_text}_nonAR.nc
output_clim=$output_dir/AR_frequency_statistics_${yrng_text}_clim.nc



mask_file="ERA5_mask.nc"
if [ ! -f "$mask_file" ]; then
    echo "Mask file $mask_file does not exist. Generating now..."
    python3 make_mask_ERA5.py
fi


if [ ! -f "$output_AR" ]; then
    python3 AR_frequency_by_map.py --output $output_AR --mask ERA5_mask.nc \
        --beg-year $beg_year \
        --end-year $end_year \
        --lat-rng ${spatial_rngs[0]} ${spatial_rngs[1]} \
        --lon-rng ${spatial_rngs[2]} ${spatial_rngs[3]} \
        --IVT-rng 250.0 1e9
fi


if [ ! -f "$output_nonAR" ]; then
    python3 AR_frequency_by_map.py --output $output_nonAR --mask ERA5_mask.nc \
        --beg-year $beg_year \
        --end-year $end_year \
        --lat-rng ${spatial_rngs[0]} ${spatial_rngs[1]} \
        --lon-rng ${spatial_rngs[2]} ${spatial_rngs[3]} \
        --IVT-rng 0.0 250.0
fi

if [ ! -f "$output_clim" ]; then
    python3 AR_frequency_by_map.py --output $output_clim --mask ERA5_mask.nc \
        --beg-year $beg_year \
        --end-year $end_year \
        --lat-rng ${spatial_rngs[0]} ${spatial_rngs[1]} \
        --lon-rng ${spatial_rngs[2]} ${spatial_rngs[3]} \
        --IVT-rng 0.0 1e9
fi


python3 plot_AR_frequency_map.py \
    --lat-rng 30 65 \
    --lon-rng -160 -110 \
    --input  $output_AR $output_clim \
    --output-dir $output_dir_AR_vs_clim \
    --no-display

python3 plot_AR_frequency_map.py \
    --lat-rng 30 65 \
    --lon-rng -160 -110 \
    --input  $output_AR $output_nonAR \
    --output-dir $output_dir_AR_vs_nonAR \
    --no-display


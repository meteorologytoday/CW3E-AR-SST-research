#!/bin/bash

. pretty_latlon.sh


#output_AR=AR_frequency_statistics.nc
#output_nonAR=ARnon_frequency_statistics.nc


beg_year=1993
end_year=2022

beg_year=1997
end_year=2018


spatial_rngs=(
    30 65 -160 -110
)

yrng_text=$( printf "%04d-%04d" $beg_year $end_year )
output_AR=AR_frequency_statistics_${yrng_text}_AR.nc
output_nonAR=AR_frequency_statistics_${yrng_text}_nonAR.nc



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

python3 plot_AR_frequency_map.py \
    --lat-rng 30 65 \
    --lon-rng -160 -110 \
    --input  $output_AR $output_nonAR


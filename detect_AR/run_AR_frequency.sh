#!/bin/bash

. pretty_latlon.sh


output_1=AR_frequency_statistics.nc
output_2=ARnon_frequency_statistics.nc

beg_year=1993
end_year=2022

spatial_rngs=(
    30 65 -160 -110
)

mask_file="ERA5_mask.nc"
if [ ! -f "$mask_file" ]; then
    echo "Mask file $mask_file does not exist. Generating now..."
    python3 make_mask_ERA5.py
fi


python3 AR_frequency_by_map.py --output $output_2 --mask ERA5_mask.nc \
    --beg-year $beg_year \
    --end-year $end_year \
    --lat-rng ${spatial_rngs[0]} ${spatial_rngs[1]} \
    --lon-rng ${spatial_rngs[2]} ${spatial_rngs[3]} \
    --IVT-rng 0.0 250.0

python3 AR_frequency_by_map.py --output $output_1 --mask ERA5_mask.nc \
    --beg-year $beg_year \
    --end-year $end_year \
    --lat-rng ${spatial_rngs[0]} ${spatial_rngs[1]} \
    --lon-rng ${spatial_rngs[2]} ${spatial_rngs[3]} \
    --IVT-rng 250.0 1e9

ncdiff -O $output_1 $output_2 ARdiff.nc

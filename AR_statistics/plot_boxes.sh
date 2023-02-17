#!/bin/bash

. pretty_latlon.sh


# format: lat_m lat_M lon_m lon_M lat_nbox lon_nbox
spatial_rngs=(
    10 60 -160 -120  10 8
    48 52 -147 -143  1  1
)

spatial_rngs=(
    10 60 120 -120  10 24
)



nparms=6


for i in $( seq 1 $(( "${#spatial_rngs[@]}" / $nparms )) ); do

    lat_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 0 ))]}
    lat_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 1 ))]}
    lon_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 2 ))]}
    lon_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 3 ))]}

    lat_nbox=${spatial_rngs[$(( ( i - 1 ) * $nparms + 4 ))]}
    lon_nbox=${spatial_rngs[$(( ( i - 1 ) * $nparms + 5 ))]}


    python3 plot_boxes.py \
        --lat-rng $lat_min $lat_max \
        --lon-rng $lon_min $lon_max \
        --lat-nbox $lat_nbox \
        --lon-nbox $lon_nbox \
        --plot-lat-rng 0 70         \
        --plot-lon-rng 100 270     \
        --output fig_boxmap_$i.png 

done

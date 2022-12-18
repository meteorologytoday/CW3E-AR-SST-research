#!/bin/bash

. pretty_latlon.sh

beg_year=2001
end_year=2022



spatial_rng=(
    30 50 -160 -130
    30 40 -160 -145
    30 40 -145 -130
    40 50 -160 -145
    40 50 -145 -130
)

if [ ] ; then
python3 count_days_map.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --output $output_dir/AR_days.nc    &
fi

for i in $( seq 1 $(( "${#spatial_rng[@]}" / 4 )) ); do

    lat_min=${spatial_rng[$(( ( i - 1 ) * 4 + 0 ))]}
    lat_max=${spatial_rng[$(( ( i - 1 ) * 4 + 1 ))]}
    lon_min=${spatial_rng[$(( ( i - 1 ) * 4 + 2 ))]}
    lon_max=${spatial_rng[$(( ( i - 1 ) * 4 + 3 ))]}

    time_str=$( printf "%04d-%04d" $beg_year $end_year )
    spatial_str=$( printf "%s-%s_%s-%s" $( pretty_lat $lat_min ) $( pretty_lat $lat_max ) $( pretty_lon $lon_min ) $( pretty_lon $lon_max )  )
    output_dir=output/${time_str}_${spatial_str}

    echo "time_str    : $time_str"    
    echo "spatial_str : $spatial_str"
    echo "output_dir  : $output_dir"

    mkdir -p $output_dir


    python3 plot_rectangular.py \
        --lat-rng $lat_min $lat_max \
        --lon-rng $lon_min $lon_max \
        --plot-lat-rng 0 70         \
        --plot-lon-rng 180 270      \
        --output $output_dir/AR_region.png  &


    python3 construct_timeseries.py \
        --beg-year=$beg_year \
        --end-year=$end_year \
        --lat-rng $lat_min $lat_max \
        --lon-rng $lon_min $lon_max \
        --output-dir $output_dir

done

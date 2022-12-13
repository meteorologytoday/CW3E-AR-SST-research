#!/bin/bash

beg_date=2015-05-01
end_date=2020-05-01
end_date=2015-12-31
#
lat_min=35
lat_max=40
#lon_min=220
lon_min=$(( 360 - 140 ))
lon_max=$(( 360 - 135 ))

output_dir=output

mkdir -p $output_dir

if [ ] ; then
python3 plot_rectangular.py \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --plot-lat-rng 0 70         \
    --plot-lon-rng 180 270      \
    --output AR_region.png  &


python3 count_days_map.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --output $output_dir/AR_days.nc

fi
python3 construct_timeseries.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --output $output_dir/AR.nc

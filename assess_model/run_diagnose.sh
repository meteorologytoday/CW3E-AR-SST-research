#!/bin/bash

beg_date=2017-11-01
end_date=2018-03-01
#end_date=2017-12-01

lat_min=20
lat_max=50
#lon_min=220
lon_min=$(( 360 - 170 ))
lon_max=$(( 360 - 130 ))

AR_npz=AR.npz
AR_glue_threshold=24
ts_npz_f120="output_ts/${beg_date}_${end_date}_f120.npz"
ts_npz_f240="output_ts/${beg_date}_${end_date}_f240.npz"

#if [ ] ; then
python3 plot_rectangular.py \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --plot-lat-rng 0 70         \
    --plot-lon-rng 180 270      \
    --output AR_region.png  &
#fi



#if [ ] ; then
python3 detect_AR.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --save-npz $AR_npz \

#if [ ] ; then

python3 compute_timeseries.py   \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --fcst=120           \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --save

python3 compute_timeseries.py   \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --fcst=240           \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --save


#fi

python3 plot_skills_and_AR.py \
    --AR-npz $AR_npz \
    --ts-npz "5d" $ts_npz_f120 "10d" $ts_npz_f240 \
    --lines  "blue" "solid" "red" "solid"         \
    --AR-glue-threshold $AR_glue_threshold \
    --output AR_ts_${beg_date}_${end_date}.png   
#fi

wait 

#!/bin/bash

fcst=240
method_skill=integrated
#method_skill=pt-by-pt

beg_date=2018-11-01
end_date=2019-03-01

lat_min=30
lat_max=50
#lon_min=220
lon_min=210
lon_max=235

AR_npz=AR.npz
skill_npz=skill

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

python3 plot_skills.py   \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --fcst=$fcst         \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --save-npz $skill_npz \
    --no-display

#fi

python3 plot_skills_and_AR.py \
    --AR-npz=$AR_npz \
    --skill-npz=$skill_npz.$method_skill.npz  \
    --output=AR_skill.png   
#fi

wait 

#!/bin/bash

AR_npz=AR.npz

skill_npz=skill


method_skill=integrated
#method_skill=pt-by-pt

beg_date=2017-01-01
end_date=2017-04-01

lat_min=30
lat_max=40
lon_min=220
lon_max=235
fcst=240

#if [ ] ; then
python3 detect_AR.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --lat-rng $lat_min $lat_max \
    --lon-rng $lon_min $lon_max \
    --save-npz $AR_npz \



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
    --skill-npz=$skill_npz.$method_skill.npz 
    

#!/bin/bash

download_dir=./data/ECCO


start_date="1996-09-30T00:00:00Z"
end_date="2018-01-01T00:00:00Z"
#end_date="1997-09-30T00:00:00Z"


datasets=(
    ECCO_L4_MIXED_LAYER_DEPTH_05DEG_DAILY_V4R4
    ECCO_L4_OCEAN_VEL_05DEG_DAILY_V4R4
    ECCO_L4_SSH_05DEG_DAILY_V4R4B 
    ECCO_L4_MIXED_LAYER_DEPTH_05DEG_DAILY_V4R4
    ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4
)

datasets=(
    ECCO_L4_MIXED_LAYER_DEPTH_05DEG_DAILY_V4R4
)

for dataset in "${datasets[@]}"; do

    echo "# Download dataset: $dataset"
    echo "# Start date: $start_date"
    echo "# End   date: $end_date"

    podaac-data-subscriber          \
        -c $dataset                 \
        -d $download_dir/$dataset   \
        --start-date $start_date    \
        --end-date $end_date


done

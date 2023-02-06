#!/bin/bash

download_dir=./data/ECCO_LLC


start_date="1996-09-30T00:00:00Z"
end_date="2007-01-01T00:00:00Z"


datasets_LLC=(
    ECCO_L4_SSH_LLC0090GRID_SNAPSHOT_V4R4
    ECCO_L4_GEOMETRY_LLC0090GRID_V4R4
    ECCO_L4_MIXED_LAYER_DEPTH_LLC0090GRID_DAILY_V4R4
    ECCO_L4_HEAT_FLUX_LLC0090GRID_DAILY_V4R4
    ECCO_L4_FRESH_FLUX_LLC0090GRID_DAILY_V4R4
    ECCO_L4_SSH_LLC0090GRID_DAILY_V4R4
    ECCO_L4_TEMP_SALINITY_LLC0090GRID_SNAPSHOT_V4R4
    ECCO_L4_STRESS_LLC0090GRID_DAILY_V4R4
    ECCO_L4_OCEAN_3D_TEMPERATURE_FLUX_LLC0090GRID_DAILY_V4R4
    ECCO_L4_OCEAN_3D_SALINITY_FLUX_LLC0090GRID_DAILY_V4R4
    ECCO_L4_OCEAN_VEL_LLC0090GRID_DAILY_V4R4
    ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_DAILY_V4R4
)

for dataset in "${datasets_LLC[@]}"; do

    echo "# Download dataset: $dataset"
    echo "# Start date: $start_date"
    echo "# End   date: $end_date"

    podaac-data-downloader          \
        -c $dataset                 \
        -d $download_dir/$dataset   \
        --start-date $start_date    \
        --end-date $end_date


done

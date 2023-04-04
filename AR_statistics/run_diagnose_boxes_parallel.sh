#!/bin/bash

. pretty_latlon.sh

ncpu=26
ocn_dataset="ECCO"

        
output_root=output_${ocn_dataset}

beg_year=1993
end_year=2017



# format: lat_m lat_M lon_m lon_M lat_nbox lon_nbox
spatial_rngs=(
    10 60 100 -100  25 80
)

# Test box
spatial_rngs=(
    30 50 170 -170 1 1
)


nparms=6


mask_ERA5="mask_ERA5.nc"
mask_ECCO="mask_ECCO.nc"

if [ ! -f "$mask_ECCO" ]; then
    echo "Mask file $mask_ECCO does not exist. Generating now..."
    #python3 make_mask_ECCO.py
fi

if [ ! -f "$mask_ERA5" ]; then
    echo "Mask file $mask_ERA5 does not exist. Generating now..."
    #python3 make_mask_ERA5.py
fi

#for fixed_500m_opt in "" "--fixed-500m" ; do
for fixed_500m_opt in "--fixed-500m" ; do

for i in $( seq 1 $(( ${#spatial_rngs[@]} / $nparms )) ); do

    lat_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 0 ))]}
    lat_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 1 ))]}
    lon_min=${spatial_rngs[$(( ( i - 1 ) * $nparms + 2 ))]}
    lon_max=${spatial_rngs[$(( ( i - 1 ) * $nparms + 3 ))]}

    lat_nbox=${spatial_rngs[$(( ( i - 1 ) * $nparms + 4 ))]}
    lon_nbox=${spatial_rngs[$(( ( i - 1 ) * $nparms + 5 ))]}


    time_str=$( printf "%04d-%04d" $beg_year $end_year )
    spatial_str=$( printf "%s-%s-n%d_%s-%s-n%d" $( pretty_lat $lat_min ) $( pretty_lat $lat_max ) $lat_nbox $( pretty_lon $lon_min ) $( pretty_lon $lon_max ) $lon_nbox )


    if [ "$fixed_500m_opt" = "--fixed-500m" ]; then
    
        extra_suffix="_500m"

    else
        extra_suffix=""
    fi

    output_dir=$output_root/${time_str}_${spatial_str}${extra_suffix}

    echo "time_str    : $time_str"    
    echo "spatial_str : $spatial_str"
    echo "extra_suffix: $extra_suffix"
    echo "output_dir  : $output_dir"

    mkdir -p $output_dir

    eval "python3 construct_timeseries_by_boxes_parallel.py \\
        --beg-year=$beg_year \\
        --end-year=$end_year \\
        --lat-rng $lat_min $lat_max \\
        --lon-rng $lon_min $lon_max \\
        --lat-nbox $lat_nbox \\
        --lon-nbox $lon_nbox \\
        --mask-ERA5 $mask_ERA5 \\
        --mask-ECCO $mask_ECCO \\
        --ignore-empty-box \\
        --output-dir $output_dir \\
        --ncpu $ncpu  \\
        $fixed_500m_opt
    "
    wait

done
done

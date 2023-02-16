#!/bin/bash

. pretty_latlon.sh


#ocn_dataset="ORA5-clim"
ocn_dataset="ECCO"

        
output_root=output_${ocn_dataset}

beg_year=1998
end_year=2017

spatial_rngs=(
    48 52 -147 -143 
    30 50 -160 -130
    40 50 -150 -130
    40 50 -145 -130
    30 40 -160 -145
    30 40 -145 -130
    40 50 -160 -145
    31 43 230 244
)



mask_ERA5="mask_ERA5.nc"
mask_ECCO="mask_ECCO.nc"

if [ ! -f "$mask_ECCO" ]; then
    echo "Mask file $mask_ECCO does not exist. Generating now..."
    python3 make_mask_ECCO.py
fi

if [ ! -f "$mask_ERA5" ]; then
    echo "Mask file $mask_ERA5 does not exist. Generating now..."
    python3 make_mask_ERA5.py
fi

if [ ] ; then
python3 count_days_map.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --output $output_dir/AR_days.nc    &
fi

for i in $( seq 1 $(( "${#spatial_rngs[@]}" / 4 )) ); do

    lat_min=${spatial_rngs[$(( ( i - 1 ) * 4 + 0 ))]}
    lat_max=${spatial_rngs[$(( ( i - 1 ) * 4 + 1 ))]}
    lon_min=${spatial_rngs[$(( ( i - 1 ) * 4 + 2 ))]}
    lon_max=${spatial_rngs[$(( ( i - 1 ) * 4 + 3 ))]}

    time_str=$( printf "%04d-%04d" $beg_year $end_year )
    spatial_str=$( printf "%s-%s_%s-%s" $( pretty_lat $lat_min ) $( pretty_lat $lat_max ) $( pretty_lon $lon_min ) $( pretty_lon $lon_max )  )
    output_dir=$output_root/${time_str}_${spatial_str}

    echo "time_str    : $time_str"    
    echo "spatial_str : $spatial_str"
    echo "output_dir  : $output_dir"

    mkdir -p $output_dir

    output_AR_file=$output_dir/AR_timeseries.nc
    if [ ! -f "$output_AR_file" ] ; then
        python3 construct_timeseries_point_by_point.py \
            --beg-year=$beg_year \
            --end-year=$end_year \
            --lat-rng $lat_min $lat_max \
            --lon-rng $lon_min $lon_max \
            --mask-ERA5 $mask_ERA5 \
            --mask-ECCO $mask_ECCO \
            --output-dir $output_dir
    fi

    output_img=$output_dir/fig_dTdt_pdf_AR.png
    python3 plot_dTdt_pdf.py --input $output_AR_file \
        --IVT-rng 250 1e5 \
        --output $output_img \
        --no-display
 
    output_img=$output_dir/fig_dTdt_pdf_AR_free.png
    python3 plot_dTdt_pdf.py --input $output_AR_file \
        --IVT-rng 0 250 \
        --output $output_img \
        --no-display
 
    output_img=$output_dir/fig_dTdt_scatter_AR_123.png
    python3 plot_dTdt_scatter_by_ARday.py \
        --input $output_AR_file \
        --IVT-rng 250 1e5 \
        --output $output_img \
        --watermonths 1 2 3 \
        --no-display

    output_img=$output_dir/fig_dTdt_scatter_AR_456.png
    python3 plot_dTdt_scatter_by_ARday.py \
        --input $output_AR_file \
        --IVT-rng 250 1e5 \
        --output $output_img \
        --watermonths 4 5 6 \
        --no-display


    output_img=$output_dir/fig_dTdt_scatter_AR_free_123.png
    python3 plot_dTdt_scatter_by_ARday.py \
        --input $output_AR_file \
        --IVT-rng 0 250 \
        --output $output_img \
        --watermonths 1 2 3 \
        --no-display

    output_img=$output_dir/fig_dTdt_scatter_AR_free_456.png
    python3 plot_dTdt_scatter_by_ARday.py \
        --input $output_AR_file \
        --IVT-rng 0 250 \
        --output $output_img \
        --watermonths 4 5 6 \
        --no-display



done

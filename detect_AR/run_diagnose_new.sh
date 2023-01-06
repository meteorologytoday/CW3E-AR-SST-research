#!/bin/bash

. pretty_latlon.sh

beg_year=2001
end_year=2014



spatial_rngs=(
    30 50 -160 -130
    30 40 -160 -145
    30 40 -145 -130
    40 50 -160 -145
    40 50 -145 -130
)

AR_dt_rngs=(
    0 50
)

# Testing

#if [ ] ; then
beg_year=2001
end_year=2014

spatial_rngs=(
    31 43 230 244
)

AR_dt_rngs=(
   0 50
)
#fi

mask_file="ERA5_mask.nc"
if [ ! -f "$mask_file" ]; then
    echo "Mask file $mask_file does not exist. Generating now..."
    python3 make_mask_ERA5.py
fi

if [ ] ; then
python3 count_days_map.py \
    --beg-date=$beg_date \
    --end-date=$end_date \
    --output $output_dir/AR_days.nc    &
fi

for mld in somxl030  ; do

    for i in $( seq 1 $(( "${#spatial_rngs[@]}" / 4 )) ); do

        lat_min=${spatial_rngs[$(( ( i - 1 ) * 4 + 0 ))]}
        lat_max=${spatial_rngs[$(( ( i - 1 ) * 4 + 1 ))]}
        lon_min=${spatial_rngs[$(( ( i - 1 ) * 4 + 2 ))]}
        lon_max=${spatial_rngs[$(( ( i - 1 ) * 4 + 3 ))]}

        time_str=$( printf "%04d-%04d" $beg_year $end_year )
        spatial_str=$( printf "%s-%s_%s-%s" $( pretty_lat $lat_min ) $( pretty_lat $lat_max ) $( pretty_lon $lon_min ) $( pretty_lon $lon_max )  )
        mld_str="mld${mld}"
        output_dir=output/${time_str}_${spatial_str}_${mld_str}

        echo "time_str    : $time_str"    
        echo "spatial_str : $spatial_str"
        echo "mld_str     : $mld_str"
        echo "output_dir  : $output_dir"

        mkdir -p $output_dir


        output_AR_region=$output_dir/AR_region.png
        output_AR_file=$output_dir/AR_timeseries_climanom.nc


        if [ ! -f "$output_AR_region" ] ; then
            python3 plot_rectangular.py \
                --lat-rng $lat_min $lat_max \
                --lon-rng $lon_min $lon_max \
                --plot-lat-rng 0 70         \
                --plot-lon-rng 180 270      \
                --output $output_dir/AR_region.png  &
        fi

        if [ ! -f "$output_AR_file" ] ; then
            python3 construct_timeseries_point_by_point.py \
                --beg-year=$beg_year \
                --end-year=$end_year \
                --lat-rng $lat_min $lat_max \
                --lon-rng $lon_min $lon_max \
                --mld $mld \
                --mask $mask_file \
                --output-dir $output_dir
        fi
       
        #if [ ]; then 
        #if [ ] ; then
        output_AR_evts_database=$output_dir/AR_evts.csv
        echo "Generating analysis: $output_AR_evts_database"
        python3 dTdt_decomposition_and_output.py --input $output_AR_file --output "$output_AR_evts_database"
        
        #fi

        for j in $( seq 1 $(( "${#AR_dt_rngs[@]}" / 2 )) ); do

            AR_dt_min=${AR_dt_rngs[$(( ( j - 1 ) * 2 + 0 ))]}
            AR_dt_max=${AR_dt_rngs[$(( ( j - 1 ) * 2 + 1 ))]}
            
            output_AR_analysis_fig=$output_dir/analysis_${AR_dt_min}-${AR_dt_max}.png

            echo "Generating analysis: $output_AR_analysis_fig"
            python3 dTdt_analysis.py --input $output_AR_file --AR-dt-rng $AR_dt_min $AR_dt_max --output "$output_AR_analysis_fig" --no-display

        done

        #fi
    done

done

#!/bin/bash

input_dir="data/ORA5/processed_remapped"


avg_yr_beg=2001
avg_yr_end=2014
output_dirname=$( printf "data/ORA5/processed_remapped_clim_%04d-%04d" $avg_yr_beg $avg_yr_end )

prefixes=(
    ORA5_NillerKrausMixedLayerDynamics_somxl010_
    ORA5_NillerKrausMixedLayerDynamics_somxl030_
)


mkdir -p $output_dirname


echo "Output dirname : $output_dirname"

for prefix in "${prefixes[@]}"; do

    echo "Doing prefix: $prefix"

    for mon in $( seq 1 12 ) ; do



        mon_str=$( printf "%02d" $mon )
        yrng_str=$( printf "{%04d..%04d}" $avg_yr_beg $avg_yr_end )
        
        output_filename="${prefix}${mon_str}.nc"
        cmd="ncra -O data/ORA5/processed_remapped/${prefixes}${yrng_str}-${mon_str}.nc $output_dirname/$output_filename"
        echo ">> $cmd"
        eval "$cmd"
    done
done



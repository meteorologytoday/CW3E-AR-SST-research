#!/bin/bash

#iter_beg=133632
iter_beg=142272
diter=576
N=5
output_dir=figures


mkdir -p $output_dir


for i in $( seq 1 $N ); do

    iter=$(( $iter_beg + ( $i - 1 ) * $diter ))

    input_filename=$( printf "output_daily_5days/diag_%010d.nc" $iter )
    output_filename=$( printf "figures/diag_%010d.png" $iter )

    echo "[$i] $input_filename => $output_filename"

    python3 plot.py --input-file $input_filename --no-display --output $output_filename --scale-days 5



done

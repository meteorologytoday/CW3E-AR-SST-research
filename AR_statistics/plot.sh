#!/bin/bash


output_dir=figures
mkdir -p $output_dir

data_dir=output_ECCO_old/1993-2017_10N-60N-n25_120E-120W-n60

echo "Fig S1"
python3 plot_AR_basic_diagnostics.py --input output_ECCO_old/1993-2017_10N-60N-n25_120E-120W-n60/AR_simple_statistics_1993-2017.nc --output $output_dir/AR_simple_stat.png &


echo "Fig 1"
#python3 plot_G_terms_map.py --input-dir $data_dir --output AR_freq.png --no-display


echo "Fig 2"
if [ ] ; then
python3 plot_AR_freq.py \
    --input-dir $data_dir \
    --output $output_dir/AR_freq1.png \
    --beg-year 1993 \
    --end-year 2017  \
    --freq-max 0.50 \
    --ndays 1       \
    --threshold-days 1 &


python3 plot_AR_freq.py \
    --input-dir $data_dir \
    --output $output_dir/AR_freq2.png \
    --beg-year 1993 \
    --end-year 2017  \
    --freq-max 0.05 \
    --ndays 7 \
    --threshold-days 7 &

python3 plot_AR_freq.py \
    --input-dir $data_dir \
    --output $output_dir/AR_freq3.png \
    --beg-year 1993 \
    --end-year 2017 \
    --freq-max 0.05  \
    --ndays 10 \
    --threshold-days 8 \
    --markers &


fi

echo "Fig 3"
#python3 plot_G_terms_map.py --input-dir $data_dir --output $output_dir/AR_forcing_partition1.png --no-display

echo "Fig 4"
#python3 plot_G_terms_map.py --input-dir $data_dir --output $output_dir/AR_forcing_partition1.png --no-display

echo "Fig 5"

. pretty_latlon.sh
spatial_rngs=(
    33 183
    33 220
    45 183    
)

nparms=2
if [ ] ; then
for i in $( seq 1 $(( "${#spatial_rngs[@]}" / $nparms )) ); do

    lat=${spatial_rngs[$(( ( i - 1 ) * $nparms + 0 ))]}
    lon=${spatial_rngs[$(( ( i - 1 ) * $nparms + 1 ))]}
    
    for bkdn in atmocn atm ocn ; do

        output_filename="$output_dir/Gterms_$( pretty_lat $lat )_$( pretty_lon $lon )_${bkdn}.png"
        python3 plot_G_terms_point.py \
            --input-dir $data_dir     \
            --lat $lat                \
            --lon $lon                \
            --breakdown $bkdn         \
            --title-style latlon      \
            --output $output_filename &
        
        
    done
done
fi

wait

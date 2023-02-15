#!/bin/bash

output_dir=figures
#input_file_AR_stat=output_ECCO/1998-2017_30N-50N_160W-130W/AR_timeseries.nc
input_file_AR_stat=output_ECCO/1998-2017_48N-52N_147W-143W/AR_timeseries.nc

rm -rf $output_dir
mkdir -p $output_dir


python3 plot_dTdt_scatter_by_ARday_frc_nonfrc.py --input $input_file_AR_stat --output $output_dir/Fig2.png & 
python3 plot_dTdt_scatter_by_ARday_frc_nonfrc.py --watermonths 1 2 3 --input $input_file_AR_stat --output $output_dir/Fig2_123.png & 
python3 plot_dTdt_scatter_by_ARday_frc_nonfrc.py --watermonths 4 5 6 --input $input_file_AR_stat --output $output_dir/Fig2_456.png & 



python3 plot_G_terms_monthly.py --input $input_file_AR_stat --breakdown atmocn --output $output_dir/Fig1_atmocn.png &
python3 plot_G_terms_monthly.py --input $input_file_AR_stat --breakdown atm    --output $output_dir/Fig1_atm.png &
python3 plot_G_terms_monthly.py --input $input_file_AR_stat --breakdown ocn    --output $output_dir/Fig1_ocn.png &


python3 plot_dTdt_scatter_by_ARday_dTdzh_vdiff.py --input $input_file_AR_stat --output $output_dir/Fig3.png &

wait

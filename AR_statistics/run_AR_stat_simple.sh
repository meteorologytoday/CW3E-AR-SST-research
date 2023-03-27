#!/bin/bash

ncpu=50
beg_year=1993
end_year=2017

input_dir=output_ECCO/1993-2017_10N-60N-n25_120E-120W-n60
yrng_str="${beg_year}-${end_year}"



#python3 mk_IVT_IWV_stat.py --input $input_dir --beg-year $beg_year --end-year $end_year

python3 mk_AR_stat.py --input-dir $input_dir --beg-year $beg_year --end-year $end_year --ncpu $ncpu --overwrite

#python3 mk_AR_interannual_stat.py --input $input_dir --beg-year $beg_year --end-year $end_year
#python3 mk_AR_EOF.py --input $input_dir/AR_interannual_statistics_${yrng_str}.nc --output-dir $input_dir


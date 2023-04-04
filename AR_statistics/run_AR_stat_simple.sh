#!/bin/bash

ncpu=50
beg_year=1993
end_year=2017

input_dir=output_ECCO/1993-2017_10N-60N-n25_100E-100W-n80
yrng_str="${beg_year}-${end_year}"

AR_algo="ANOM_LEN"

#for suffix in "" "_500m"; do
for suffix in "_500m"; do
    _input_dir="${input_dir}${suffix}"
    python3 mk_AR_stat.py --input-dir $_input_dir --beg-year $beg_year --end-year $end_year --ncpu $ncpu --AR-algo $AR_algo
    python3 mk_AR_interannual_stat.py --input $_input_dir --beg-year $beg_year --end-year $end_year --AR-algo $AR_algo
    python3 mk_AR_EOF.py --input $_input_dir/AR_interannual_statistics_${AR_algo}_${yrng_str}.nc --output-dir $input_dir
done
#python3 mk_AR_stat.py --input-dir $input_dir --beg-year $beg_year --end-year $end_year --ncpu $ncpu --overwrite --AR-algo $AR_algo --fixed-500m

#python3 mk_IVT_IWV_stat.py --input $input_dir --beg-year $beg_year --end-year $end_year







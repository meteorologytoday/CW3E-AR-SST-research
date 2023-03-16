#!/bin/bash

beg_year=1993
end_year=2017

input_dir=output_ECCO/1993-2017_10N-60N-n25_120E-120W-n60

python3 mk_IVT_IWV_stat.py --input $input_dir --beg-year $beg_year --end-year $end_year

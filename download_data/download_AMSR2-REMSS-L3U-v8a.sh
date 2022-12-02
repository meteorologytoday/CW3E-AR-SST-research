#!/bin/bash

cp shared_header_login_PODAAC_netrc.txt ~/.netrc

output_dir=data/AMSR2/SST
beg_datetime=2017-01-01T00:00:00Z
end_datetime=2019-01-01T00:00:00Z


mkdir -p $output_dir

podaac-data-downloader \
    -c AMSR2-REMSS-L3U-v8a \
    -d ./data/AMSR2-REMSS-L3U-v8a \
    --start-date "$beg_datetime" \
    --end-date   "$end_datetime" \
    -d $output_dir

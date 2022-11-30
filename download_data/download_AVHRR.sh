#!/bin/bash

cp shared_header_login_PODAAC_netrc.txt ~/.netrc

beg_datetime=2017-01-01T00:00:00Z
end_datetime=2019-01-01T00:00:00Z

podaac-data-downloader \
    -c AVHRR_OI-NCEI-L4-GLOB-v2.0 \
    -d ./data/AVHRR_OI-NCEI-L4-GLOB-v2.0 \
    --start-date "$beg_datetime" \
    --end-date   "$end_datetime" \
    -d data_AVHRR

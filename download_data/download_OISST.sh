#!/bin/bash


download_folder="data/OISST"

mkdir -p $download_folder

for y in $( seq 2016 2018 ) ; do
   
    filename=$( printf "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.%04d.nc" $y ) 

    echo "Download file: $filename"
    wget $filename -P $download_folder

done

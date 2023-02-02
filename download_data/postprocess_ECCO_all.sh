#!/bin/bash


for MLD_method in RHO ECCO ; do
    python3 postprocess_ECCO.py --MLD-method $MLD_method
done
    

./postprocess_remap_ECCO.sh

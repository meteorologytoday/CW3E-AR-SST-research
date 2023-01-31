

# ECCOv4 data

- Download data: run `./download_ECCO.sh` which uses `podaac-data-subscriber` to download data.
- Generate mixed layer temperature, salinity, and depth: Run `python3 postprocess_ECCO.py`.
- To generate ERA5 compatible grid resolution (i.e. 0.5 deg regrid to 0.25 deg), run `postprocess_remap_ECCO.sh`.
- Lazy option: `postprocess_ECCO_all.sh` do all MLD methods (`RHO` and `ECCO`) and remapping automatically.

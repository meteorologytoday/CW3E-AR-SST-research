using PyCall
using JLD2

include("MITgcmTools.jl")
include("Operators.jl")
include("Operators_ML.jl")

function findArgRange(
    arr :: AbstractArray{T, 1},
    lb :: T,
    ub :: T,
) where T

    if lb > ub
        throw(ErrorException("Lower bound should be no larger than upper bound"))
    end


    if any( (arr[2:end] - arr[1:end-1]) .<= 0 )
        throw(ErrorException("input array should be monotonically increasing"))
    end

    idx = lb .<= arr .<= ub
    
    idx_low = findfirst(idx)
    idx_max = findlast(idx)

    return idx_low, idx_max
end


println("*** Test SurfaceTendency package ***")

mitgcm = pyimport("MITgcmutils")

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"

lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0]

coo_tmp = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true)

lat_idx = findArgRange(coo_tmp.gd.y_T[1, :, 1], lat_rng[1], lat_rng[2])
lon_idx = findArgRange(coo_tmp.gd.x_T[:, 1, 1], lon_rng[1], lon_rng[2])

mitgcm_region = (lon_idx[1]-1, lon_idx[2], lat_idx[1]-1, lat_idx[2])

coo = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true, lev=lev, region=mitgcm_region)

#iters = 386400
iters = 388800
data_3D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iters, region=mitgcm_region, lev=collect(lev), returnmeta=true))
data_2D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iters, region=mitgcm_region, returnmeta=true))


println("Varnames of data_3D: ", keys(data_3D))
println("Varnames of data_2D: ", keys(data_2D))






println("Loading data...")

TEMP   = data_3D["THETA"]
SALT   = data_3D["SALT"]
KPPhbl = data_2D["KPPhbl"]

MLD, bundle = Operators_ML.detectMLD(TEMP, SALT, coo)
MLD_compromised = maxarr(MLD, KPPhbl)

using NCDatasets
output_file = "output_MLD.nc"
println("Writing output: $output_file")
Dataset(output_file, "c") do ds

    defDim(ds, "lon", coo.gd.Nx)
    defDim(ds, "lat", coo.gd.Ny)
    defDim(ds, "z",   coo.gd.Nz)

    # Define the variables temperature with the attribute units

    for (varname, vardata, datatype, dimnames) in [
        ("z",      coo.gd.z_T[:], eltype(coo.gd.z_T), ("z",)),
        ("rho",    bundle["ρ"], eltype(bundle["ρ"]), ("lon", "lat", "z",)),
        ("MLD",    MLD,         eltype(MLD),         ("lon", "lat",)),
        ("KPPhbl", KPPhbl,      eltype(KPPhbl),      ("lon", "lat",)),
        ("MLD_compromised", MLD_compromised, eltype(MLD_compromised),  ("lon", "lat",)),
    ]

        _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
        _var[:] = vardata
    end

end

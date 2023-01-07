using PyCall
using JLD2

ρ   = 1027.5
c_p = 3994.0
fill_value = 1e20

include("MITgcmTools.jl")
include("Operators.jl")
include("Operators_ML.jl")

function py2jl(x::AbstractArray)

    dim_len = length(size(x))
    permute = [dim_len - (i-1) for i=1:dim_len]

    return permutedims(x, permute)
end

function maxarr(
    a :: AbstractArray{T},
    b :: AbstractArray{T},
) where T

    m = zeros(T, size(a)...)

    for i=1:length(a)
        m[i] = max(a[i], b[i])
    end

    return m
end

println("*** Testing mitgcm ***")

mitgcm = pyimport("MITgcmutils")

#data_dir = "/data/SO2/SWOT/MARA/RUN4_LY_NoRainOct11to18/DIAGS"
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"


coo = readMITgcmGrid_MM(grid_dir, verbose=true)

lev = 1:50
#iters = 386400
iters = 388800
data_3D = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iters, lev=collect(lev), returnmeta=True))
data_2D = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iters,    lev=collect(lev), returnmeta=True))


println("Size of data_3D: ", size(data_3D))
println("Size of data_2D: ", size(data_2D))






println("Loading data...")

TEMP   = data_3D["THETA"]
SALT   = data_3D["SALIN"]
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

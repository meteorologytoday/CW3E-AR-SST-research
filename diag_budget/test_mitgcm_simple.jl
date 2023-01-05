using PyCall
using JLD2

fill_value = 1e20

include("read_mitgcm_grid.jl")
include("Operators.jl")

function py2jl(x::AbstractArray)

    dim_len = length(size(x))
    permute = [dim_len - (i-1) for i=1:dim_len]

    return permutedims(x, permute)
end


println("*** Testing mitgcm ***")

mitgcm = pyimport("MITgcmutils")

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"


coo = readMITgcmGrid_MM(grid_dir, verbose=true)

data = mitgcm.mds.rdmds("$data_dir/diag_Tbdgt", 512064) |> py2jl

println("Size of data: ", size(data))

TOTTTEND = data[:, :, :, 5]
dTOTTTENDdx = Operators.T_mask_T(Operators.T_interp_U( Operators.U_∂x_T(TOTTTEND, coo), coo), coo; fill_value=fill_value)

println(size(dTOTTTENDdx))

using NCDatasets

Dataset("test.nc", "c") do ds

    defDim(ds, "lon", coo.gd.Nx)
    defDim(ds, "lat", coo.gd.Ny)
    defDim(ds, "z",   coo.gd.Nz)

    # Define the variables temperature with the attribute units

    for (varname, vardata, datatype, dimnames) in [
        ("TOTTTEND", TOTTTEND, eltype(TOTTTEND), ("lon", "lat", "z")),
        ("dTOTTTENDdx", dTOTTTENDdx, eltype(dTOTTTENDdx), ("lon", "lat", "z")),
    ]

        _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
        _var[:] = vardata
    end

end

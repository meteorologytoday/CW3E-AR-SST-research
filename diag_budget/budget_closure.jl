using PyCall
using JLD2

ρ   = 1029.0
c_p = 3996.0
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

#data_dir = "/data/SO2/SWOT/MARA/RUN4_LY_NoRainOct11to18/DIAGS"
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"


coo = readMITgcmGrid_MM(grid_dir, verbose=true)

#iters = 386400
iters = 388800
data_3D = mitgcm.mds.rdmds("$data_dir/diag_Tbdgt", iters) |> py2jl
data_2D = mitgcm.mds.rdmds("$data_dir/diag_2D", iters) |> py2jl

println("Size of data_3D: ", size(data_3D))
println("Size of data_2D: ", size(data_2D))

function nest(a :: AbstractArray{T, 3}, grid::Symbol) where T

    Nx, Ny, Nz = size(a)

    newa = nothing

    if grid == :U
        newa = zeros(T, Nx+1, Ny, Nz)
        newa[1:Nx, :, :] = a
    elseif grid == :V
        newa = zeros(T, Nx, Ny+1, Nz)
        newa[:, 1:Ny, :] = a
    elseif grid == :W
        newa = zeros(T, Nx, Ny, Nz+1)
        newa[:, :, 1:Nz] = a

    else
        throw(ErrorException("Unknown grid: $grid"))
    end

    return newa
end

println("Loading data...")

ADVr_TH  = nest(data_3D[:, :, :, 1], :W)
ADVx_TH  = nest(data_3D[:, :, :, 2], :U)
ADVy_TH  = nest(data_3D[:, :, :, 3], :V)

DFrI_TH  = nest(data_3D[:, :, :, 4], :W)

TOTTTEND = data_3D[:, :, :, 5] / 86400.0

WTHMASS  = nest(data_3D[:, :, :, 9], :W)
KPPg_TH  = nest(data_3D[:, :, :, 10], :W)

oceQnet = data_2D[:, :,  9]
oceQnet  = reshape(oceQnet, size(oceQnet)..., 1)

oceQsw  = data_2D[:, :, 10]
oceQsw  = reshape(oceQsw, size(oceQsw)..., 1)



#=
ADVr_TH  = nest(data_3D[:, :, :, 1], :W)
ADVx_TH  = nest(data_3D[:, :, :, 2], :U)
ADVy_TH  = nest(data_3D[:, :, :, 3], :V)

DFxE_TH  = nest(data_3D[:, :, :, 4], :U)
DFyE_TH  = nest(data_3D[:, :, :, 5], :V)
DFrI_TH  = nest(data_3D[:, :, :, 6], :W)

TOTTTEND = data_3D[:, :, :, 7] / 86400.0

#UVELTH   = data_3D[:, :, :, 6]
#VVELTH   = data_3D[:, :, :, 7]
#WVELTH   = data_4D[:, :, :, 8]

WTHMASS  = nest(data_3D[:, :, :, 11], :W)
KPPg_TH  = nest(data_3D[:, :, :, 12], :W)

oceQsw  = data_2D[:, :, 10]
oceQsw = reshape(oceQsw, size(oceQsw)..., 1)
=#

println("Compute Qsw flux...")

Qsw_shape(z) = 0.62 * exp.(z/0.6) + (1 - 0.62) * exp.(z/20.0)
SWFLX = - oceQsw .* Qsw_shape(coo.gd.z_W)

SFCFLX_shape(z) = z .>= 0
SFCFLX = - (oceQnet - oceQsw) .* SFCFLX_shape(coo.gd.z_W)


println("Compute tendency...")

TEND_ADV = - (
        Operators.T_DIVx_U(ADVx_TH, coo)
    +   Operators.T_DIVy_V(ADVy_TH, coo)
    +   Operators.T_DIVz_W(ADVr_TH, coo)
)

TEND_DIFF = - (
       Operators.T_DIVz_W(DFrI_TH, coo)
)

#=
TEND_DIFF = - (
        Operators.T_DIVx_U(DFxE_TH, coo)
    +   Operators.T_DIVy_V(DFyE_TH, coo)
    +   Operators.T_DIVz_W(DFrI_TH, coo)
)
=#

TEND_KPP = - Operators.T_DIVz_W(KPPg_TH, coo)

TEND_SWFLX = - Operators.T_DIVz_W(SWFLX, coo; already_weighted=false) / (ρ*c_p)

TEND_SFCFLX = - Operators.T_DIVz_W(SFCFLX, coo; already_weighted=false) / (ρ*c_p)

TEND_SUM = TEND_ADV + TEND_DIFF + TEND_KPP + TEND_SWFLX + TEND_SFCFLX

TEND_RES = TEND_SUM - TOTTTEND

using NCDatasets
output_file = "check_budget.nc"
println("Writing output: $output_file")
Dataset(output_file, "c") do ds

    defDim(ds, "lon", coo.gd.Nx)
    defDim(ds, "lat", coo.gd.Ny)
    defDim(ds, "z",   coo.gd.Nz)

    # Define the variables temperature with the attribute units

    for (varname, vardata, datatype, dimnames) in [
        ("TOTTTEND",   TOTTTEND,   eltype(TOTTTEND),   ("lon", "lat", "z")),
        ("TEND_ADV",   TEND_ADV,   eltype(TEND_ADV),   ("lon", "lat", "z")),
        ("TEND_DIFF",  TEND_DIFF,  eltype(TEND_DIFF),  ("lon", "lat", "z")),
        ("TEND_KPP",   TEND_KPP,   eltype(TEND_KPP),   ("lon", "lat", "z")),
        ("TEND_SWFLX", TEND_SWFLX, eltype(TEND_SWFLX), ("lon", "lat", "z")),
        ("TEND_SFCFLX",TEND_SFCFLX,eltype(TEND_SFCFLX),("lon", "lat", "z")),
        ("TEND_SUM",   TEND_SUM,   eltype(TEND_SUM),   ("lon", "lat", "z")),
        ("TEND_RES",   TEND_RES,   eltype(TEND_RES),   ("lon", "lat", "z")),
    ]

        _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
        _var[:] = vardata
    end

end

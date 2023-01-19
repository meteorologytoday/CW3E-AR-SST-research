include("MITgcmTools.jl")
include("Operators.jl")

using Statistics
using PyCall
using JLD2

mitgcm = pyimport("MITgcmutils")

ρ   = 1027.5
c_p = 3994.0
fill_value = 1e20

println("*** Testing mitgcm ***")

#data_dir = "/data/SO2/SWOT/MARA/RUN4_LY_NoRainOct11to18/DIAGS"
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
grid_dir = "/data/SO2/SWOT/GRID/BIN"



lat_rng = [31.0, 43.0]
lon_rng = [230.0, 244.0] #360 .- [130.0, 116.0]

#lat_rng = [0.0, 90.0]
#lat_rng = [35.0, 37.0]
#lon_rng = [238.0, 240.0]
lev = 1:40
mitgcm_lev = collect(lev) .- 1

coo_tmp = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true)

lat_idx_rng = MITgcmTools.findArgRange(coo_tmp.gd.ϕ_T[1, :, 1], lat_rng[1], lat_rng[2])
lon_idx_rng = MITgcmTools.findArgRange(coo_tmp.gd.λ_T[:, 1, 1], lon_rng[1], lon_rng[2])

region        = (lon_idx_rng..., lat_idx_rng...)
mitgcm_region = (lon_idx_rng[1]-1, lon_idx_rng[2], lat_idx_rng[1]-1, lat_idx_rng[2])


println("# region : ", region)
println("# lev    : ", lev)

iters = 142272
#388800

coo = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true, lev=lev, region=region)
data_3D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iters, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
data_Tbdgt, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_Tbdgt", iters, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
data_2D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iters, region=mitgcm_region, returnmeta=true))

#println("Size of data_3D: ", size(data_3D))
#println("Size of data_2D: ", size(data_2D))



println("z_W: ", coo.gd.z_W[:])
println("Δz_W: ", coo.gsp.Δz_T)

println("Loading data...")

d = Dict()

mapping_grid3D = Dict(
    "ADVr_TH"  => :W,
    "ADVx_TH"  => :U,
    "ADVy_TH"  => :V,
    "DFrI_TH"  => :W,
    "TOTTTEND" => :T,
    "WTHMASS"  => :W,
    "KPPg_TH"  => :W,
)

for (varname, grid) in mapping_grid3D
    
    _target_data = nothing

    if haskey(data_3D, varname)
        _target_data = data_3D        
    elseif haskey(data_Tbdgt, varname)
        _target_data = data_Tbdgt
    else
        throw(ErrorException("Unknown varname: $varname"))
    end
    d[varname] = MITgcmTools.nest3D(_target_data[varname], grid) 
end

for varname in ["oceQnet", "oceQsw"]
    d[varname] = MITgcmTools.nestSlab(data_2D[varname])
end


d["TOTTTEND"] ./= 86400.0

println("Compute Qsw flux...")

Qsw_shape(z) = ( 0.62 * exp.(z/0.6) + (1 - 0.62) * exp.(z/20.0) ) .* (z .>= -200.0)
d["SWFLX"] = - d["oceQsw"] .* Qsw_shape(coo.gd.z_W)

SFCFLX_shape(z) = z .>= 0
d["SFCFLX"] = - (d["oceQnet"] - d["oceQsw"]) .* SFCFLX_shape(coo.gd.z_W)

println("Compute tendency...")

d["TEND_ADV_X"] = - (
        Operators.T_DIVx_U(d["ADVx_TH"], coo)
)

d["TEND_ADV_Y"] = - (
       Operators.T_DIVy_V(d["ADVy_TH"], coo)
)

d["TEND_ADV_Z"] = - (
       Operators.T_DIVz_W(d["ADVr_TH"], coo)
)

d["TEND_ADV"] = d["TEND_ADV_X"] + d["TEND_ADV_Y"] + d["TEND_ADV_Z"]


d["TEND_DIFF"] = - (
       Operators.T_DIVz_W(d["DFrI_TH"], coo)
)


global_area = 1.022854491990098E+12

d["WTHMASS_masked"] = d["WTHMASS"] .* SFCFLX_shape(coo.gd.z_W)
d["TEND_SFC_WTHMASS"] = - Operators.T_DIVz_W(d["WTHMASS_masked"], coo; already_weighted=false)

TEMP_SURF_CORR_MEAN = sum(d["WTHMASS"][:, :, 1:1] .* coo.gsp.Δa_T) / global_area
println("TEMP_SURF_CORR_MEAN: ", TEMP_SURF_CORR_MEAN)

TEMP_SURF_CORR_MEAN *= 0

d["TEND_SFC_WTHMASS"] = - Operators.T_DIVz_W( d["WTHMASS_masked"] .- TEMP_SURF_CORR_MEAN, coo; already_weighted=false)

println("Global area from summing: ", sum(coo.gsp.Δa_T[:, :, 1]))
println("Global area from mitGCM STDOUT: ", 1.022854491990098E+12)


d["TEND_KPP"]    = - Operators.T_DIVz_W(d["KPPg_TH"], coo)

d["TEND_SWFLX"]  = - Operators.T_DIVz_W(d["SWFLX"], coo; already_weighted=false) / (ρ*c_p)

d["TEND_SFCFLX"] = - Operators.T_DIVz_W(d["SFCFLX"], coo; already_weighted=false) / (ρ*c_p)

d["TEND_SUM"] = d["TEND_ADV_X"] + d["TEND_ADV_Y"] + d["TEND_ADV_Z"] + d["TEND_DIFF"] + d["TEND_SWFLX"] + d["TEND_SFCFLX"] + d["TEND_SFC_WTHMASS"]

d["TEND_RES"] = d["TEND_SUM"] - d["TOTTTEND"]

elm_type = eltype(d["TOTTTEND"])

using NCDatasets
output_file = "check_budget.nc"
println("Writing output: $output_file")
Dataset(output_file, "c") do ds

    defDim(ds, "lon", coo.gd.Nx)
    defDim(ds, "lat", coo.gd.Ny)
    defDim(ds, "z",   coo.gd.Nz)
    defDim(ds, "zp1",   coo.gd.Nz+1)

    # Define the variables temperature with the attribute units

    for (varname, vardata, datatype, dimnames) in [
        ("TOTTTEND",   d["TOTTTEND"],   elm_type, ("lon", "lat", "z")),
        ("TEND_ADV_X", d["TEND_ADV_X"], elm_type, ("lon", "lat", "z")),
        ("TEND_ADV_Y", d["TEND_ADV_Y"], elm_type, ("lon", "lat", "z")),
        ("TEND_ADV_Z", d["TEND_ADV_Z"], elm_type, ("lon", "lat", "z")),
        ("TEND_ADV",   d["TEND_ADV"],   elm_type, ("lon", "lat", "z")),
        ("TEND_DIFF",  d["TEND_DIFF"],  elm_type, ("lon", "lat", "z")),
        ("TEND_KPP",   d["TEND_KPP"],   elm_type, ("lon", "lat", "z")),
        ("TEND_SWFLX", d["TEND_SWFLX"], elm_type, ("lon", "lat", "z")),
        ("TEND_SFCFLX",d["TEND_SFCFLX"],elm_type, ("lon", "lat", "z")),
        ("TEND_SFC_WTHMASS",d["TEND_SFC_WTHMASS"],elm_type, ("lon", "lat", "z")),
        ("TEND_SUM",   d["TEND_SUM"],   elm_type, ("lon", "lat", "z")),
        ("TEND_RES",   d["TEND_RES"],   elm_type, ("lon", "lat", "z")),
        ("VDIFF",      d["DFrI_TH"],   elm_type, ("lon", "lat", "zp1")),
        ("WTHMASS",    d["WTHMASS"],   elm_type, ("lon", "lat", "zp1")),
        ("WTHMASS_masked", d["WTHMASS_masked"],   elm_type, ("lon", "lat", "zp1")),
        ("z",          coo.gd.z_W[1, 1, :], elm_type, ("zp1",)),
    ]

        _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
        _var[:] = vardata
    end

end

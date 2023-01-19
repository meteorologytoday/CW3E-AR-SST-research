include("MITgcmTools.jl")
include("SurfaceTendency.jl")

using Formatting
using NCDatasets
using PyCall
mitgcm = pyimport("MITgcmutils")

fill_value = 1e20

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


println("*** Testing SurfaceTendency.jl ***")



grid_dir = "/data/SO2/SWOT/GRID/BIN"

model_dt = 150.0

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
N = 120
diter = 576
beg_iter = 133632

#=
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY_NoRainOct11to18/DIAGS"
N = 744
diter = 24
beg_iter = 368688
=#

lat_rng = [35.0, 37.0]
lat_rng = [0.0, 90.0]
lon_rng = [230.0, 240.0]
lev = 1:40
mitgcm_lev = collect(lev) .- 1

coo_tmp = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true)

lat_idx_rng = MITgcmTools.findArgRange(coo_tmp.gd.ϕ_T[1, :, 1], lat_rng[1], lat_rng[2])
lon_idx_rng = MITgcmTools.findArgRange(coo_tmp.gd.λ_T[:, 1, 1], lon_rng[1], lon_rng[2])

region        = (lon_idx_rng..., lat_idx_rng...)
mitgcm_region = (lon_idx_rng[1]-1, lon_idx_rng[2], lat_idx_rng[1]-1, lat_idx_rng[2])


println("# region : ", region)
println("# lev    : ", lev)



process_iters = [ beg_iter + (i-1) * diter for i=1:N ]

coo = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true, lev=lev, region=region)
coo_s = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true, lev=1:1, region=region)

mapping_grid3D = Dict(
    "TOTTTEND" => :T,
    "DFrI_TH"  => :W,
    "DFxE_TH"  => :U,
    "DFyE_TH"  => :V,
    "THETA"    => :T,
    "SALT"     => :T,
    "UVEL"     => :U,
    "VVEL"     => :V,
    "WVEL"     => :W,
)


for (i, iter_now) in enumerate(process_iters)
   
    println("Iteration $i / $N : $iter_now") 
    println("Loading data...")


    data_3D_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    data_Tbdgt_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_Tbdgt", iter_now, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))


    data_2D_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now, region=mitgcm_region, returnmeta=true))


    println("Loading extra timestep")    
    
    data_3D_l, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now - diter, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    data_3D_r, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now + diter, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    
    data_2D_l, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now - diter, region=mitgcm_region, returnmeta=true))
    data_2D_r, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now + diter, region=mitgcm_region, returnmeta=true))


    d_c = Dict()
    d_l = Dict()
    d_r = Dict()

    for varname in ["TOTTTEND", "DFrI_TH"]#, "DFxE_TH", "DFyE_TH"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_Tbdgt_c[varname], grid) 
    end
    d_c["TOTTTEND"] ./= 86400.0
    #d_c["DFrI_TH"]  ./= coo.gsp.Δa_T
    #d_c["DFxE_TH"]  ./= coo.gsp.Δz_T
    #d_c["DFyE_TH"]  ./= coo.gsp.Δz_T
    
    #println("Sum of DFxE_TH: ", sum(d_c["DFxE_TH"]))

    for varname in ["UVEL", "VVEL", "WVEL"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_3D_c[varname], grid) 
    end

    for varname in ["oceQnet", "oceQsw",]
        d_c[varname] = data_2D_c[varname]
    end

    for varname in ["THETA", "SALT"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_3D_c[varname], grid) 
        d_l[varname] = MITgcmTools.nest3D(data_3D_l[varname], grid) 
        d_r[varname] = MITgcmTools.nest3D(data_3D_r[varname], grid) 
    end

    for varname in ["KPPhbl",]
        d_c[varname] = data_2D_c[varname]
        d_l[varname] = data_2D_l[varname]
        d_r[varname] = data_2D_r[varname]
    end


    # determin mixed-layer depth
    h_l_algo, _ = Operators_ML.detectMLD(d_l["THETA"], d_l["SALT"], coo)
    h_c_algo, _ = Operators_ML.detectMLD(d_c["THETA"], d_c["SALT"], coo)
    h_r_algo, _ = Operators_ML.detectMLD(d_r["THETA"], d_r["SALT"], coo)
    
    h_l = maxarr(h_l_algo, data_2D_l["KPPhbl"])
    h_c = maxarr(h_c_algo, data_2D_c["KPPhbl"])
    h_r = maxarr(h_r_algo, data_2D_r["KPPhbl"])



    #h_l .= 400.0
    #h_c .= 400.0
    #h_r .= 400.0

    TOTTTEND_mean = Operators_ML.computeMLMean(
        d_c["TOTTTEND"],
        coo,
        h = h_c,
        do_avg = true,
    )

    println("Compute SurfaceTendency terms...")
    terms = SurfaceTendency.computeSurfaceTendencyTerms(;
        X      = d_c["THETA"],
        h_l    = h_l,
        h_c    = h_c,
        h_r    = h_r,
        UVEL   = d_c["UVEL"],
        VVEL   = d_c["VVEL"],
        WVEL   = d_c["WVEL"],
#        XDIFFFLX = d_c["DFxE_TH"],
#        YDIFFFLX = d_c["DFyE_TH"],
        XDIFFFLX = d_c["UVEL"] * 0.0,
        YDIFFFLX = d_c["VVEL"] * 0.0,
        ZDIFFFLX = d_c["DFrI_TH"],
        Fsol   = d_c["oceQsw"],
        Fnet   = d_c["oceQnet"],
        Δt     = diter * model_dt,
        coo    =  coo,
        coo_s  = coo_s,
        tracer = "TEMP",
    )

        
    terms["TEND_SUM_NOVDIFF"] = terms["TEND_SFCFLX"] + terms["TEND_ENT_dhdt"] + terms["TEND_ENT_wb"] + terms["TEND_ENT_hadv"] + terms["TEND_BT_hadv"] + terms["TEND_EDDY"] + terms["TEND_HDIFF"] #+ terms["TEND_VDIFF"]
    terms["TEND_SUM"] = terms["TEND_SFCFLX"] + terms["TEND_ENT_dhdt"] + terms["TEND_ENT_wb"] + terms["TEND_ENT_hadv"] + terms["TEND_BT_hadv"] + terms["TEND_EDDY"] + terms["TEND_HDIFF"] + terms["TEND_VDIFF"]

    terms["TEND_RES_NOVDIFF"] = terms["TEND_SUM_NOVDIFF"] - TOTTTEND_mean
    terms["TEND_RES"] = terms["TEND_SUM"] - TOTTTEND_mean

    elm_type = eltype(d_c["THETA"])


    DFrI_TH_bot = Operators_ML.evalAtMLD_W(
        d_c["DFrI_TH"],
        coo,
        h = h_c,
    )


    output_file = format("output/diag_{:010d}.nc", iter_now)
    println("Writing output: $output_file")
    Dataset(output_file, "c") do ds

        defDim(ds, "lon", coo.gd.Nx)
        defDim(ds, "lat", coo.gd.Ny)
        defDim(ds, "z",   coo.gd.Nz)
        defDim(ds, "z_w", coo.gd.Nz+1)
        defDim(ds, "time", Inf)

        # Define the variables temperature with the attribute units

        for varname in keys(terms)
            _var = defVar(ds, varname, elm_type, ("lon", "lat", "time"), fillvalue=fill_value)
            _var[:, :, 1] = terms[varname]
        end

        println(coo.gd.z_W[1, 1, :])
        for (varname, vardata, datatype, dimnames) in [
            ("z_w",            -coo.gd.z_W[1, 1, :], eltype(coo.gd.z_W), ("z_w", )),
            ("z",              coo.gd.z_T[1, 1, :], elm_type, ("z", )),
            ("TOTTTEND_mean",  TOTTTEND_mean,   elm_type, ("lon", "lat",)),
            ("TOTTTEND",       d_c["TOTTTEND"], elm_type, ("lon", "lat", "z",)),
            ("THETA",          d_c["THETA"], elm_type, ("lon", "lat", "z",)),
            ("SALT",           d_c["SALT"], elm_type, ("lon", "lat", "z",)),
            ("h_c_algo",       h_c_algo, elm_type, ("lon", "lat", )),
            ("KPPhbl",         d_c["KPPhbl"], elm_type, ("lon", "lat", )),
#            ("DFrI_TH",        d_c["DFrI_TH"], elm_type, ("lon", "lat", "z_w",)),
#            ("DFrI_TH_bot",    DFrI_TH_bot, elm_type, ("lon", "lat",)),
            ("h_c",            h_c, elm_type, ("lon", "lat",)),
#            ("da",             coo.gsp.Δa_T[:, :, 1], elm_type, ("lon", "lat",)),
        ]

            colons = [Colon() for i=1:length(dimnames)]

            dimnames = [dimnames..., "time"]

            _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
            _var[colons..., 1] = vardata
        end

    end

end


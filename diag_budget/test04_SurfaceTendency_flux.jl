include("MITgcmTools.jl")
include("SurfaceTendency_flux.jl")

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

#=
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/DIAGS_DLY"
N = 360
diter = 576
#beg_iter = 133632
beg_iter = 142272
=#

data_dir = "/data/SO2/SWOT/MARA/RUN4_LY/TEST_TFLUX"
snapshot_data_dir = "/data/SO2/SWOT/MARA/RUN4_LY"
N = 1
diter = 576
beg_iter = 150336



output_dir = "output_daily_$(N)days_flux"
mkpath(output_dir)

#=
data_dir = "/data/SO2/SWOT/MARA/RUN4_LY_NoRainOct11to18/DIAGS"
N = 744
diter = 24
beg_iter = 368688
output_dir = "output_mara_hourly_$(N)"
mkpath(output_dir)
=#

lat_rng = [35.0, 37.0]
lon_rng = [235.0, 238.0]
#lat_rng = [0.0, 90.0]
#lon_rng = [230.0, 240.0]

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
    "ADVx_TH"  => :U,
    "ADVy_TH"  => :V,
    "ADVr_TH"  => :W,
    "WTHMASS"  => :W,
)


for (i, iter_now) in enumerate(process_iters)
   
    println("Iteration $i / $N : $iter_now") 

    output_file = format("$output_dir/diag_{:010d}.nc", iter_now)
    println("Target output: $output_file")

    if isfile(output_file)
        println("File $output_file already exists. Move on to the next..")
        continue
    end


    println("Loading data...")


    snapshot_3D_l, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$snapshot_data_dir/pickup", iter_now - diter, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    snapshot_3D_r, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$snapshot_data_dir/pickup", iter_now, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))


    data_3D_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    data_Tbdgt_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_Tbdgt", iter_now, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))


    data_2D_c, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now, region=mitgcm_region, returnmeta=true))


    println("Loading extra timestep")    
    
    data_3D_lc, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now - diter, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    data_3D_cr, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iter_now + diter, region=mitgcm_region, lev=mitgcm_lev, returnmeta=true))
    
    data_2D_l, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now - diter, region=mitgcm_region, returnmeta=true))
    data_2D_r, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iter_now + diter, region=mitgcm_region, returnmeta=true))


    d_c = Dict()
    d_l = Dict()
    d_r = Dict()
    d_lc = Dict()
    d_cr = Dict()

    for varname in ["TOTTTEND", "DFrI_TH", "WTHMASS"]#, "DFxE_TH", "DFyE_TH"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_Tbdgt_c[varname], grid) 
    end
    d_c["TOTTTEND"] ./= 86400.0
    
    #println("Sum of DFxE_TH: ", sum(d_c["DFxE_TH"]))

    for varname in ["ADVx_TH", "ADVy_TH", "ADVr_TH"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_Tbdgt_c[varname], grid) 
    end

    for varname in ["oceQnet", "oceQsw", "oceFWflx", "EXFuwind", "EXFvwind", "TFLUX"]
        d_c[varname] = data_2D_c[varname]
    end

    for varname in ["THETA", "SALT"]
        grid = mapping_grid3D[varname]
        d_c[varname] = MITgcmTools.nest3D(data_3D_c[varname], grid) 
        d_l[varname] = MITgcmTools.nest3D(data_3D_l[varname], grid) 
        d_r[varname] = MITgcmTools.nest3D(data_3D_r[varname], grid) 
    end

    for varname in ["Theta", "Salt"]
        grid = mapping_grid3D[varname]
        d_lc[varname] = MITgcmTools.nest3D(data_3D_lc[varname], grid) 
        d_cr[varname] = MITgcmTools.nest3D(data_3D_cr[varname], grid) 
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


    #=
    cutoff_h = - coo.gd.z_T[:, :, 30]
    h_l .= cutoff_h
    h_c .= cutoff_h
    h_r .= cutoff_h
    =#

    #h_l .= - coo.gd.z_T[:, :, 30]
    #h_c .= - coo.gd.z_T[:, :, 30]
    #h_r .= - coo.gd.z_T[:, :, 30]

    TOTTTEND_mean = Operators_ML.computeMLMean(
        d_c["TOTTTEND"],
        coo,
        h = h_c,
        do_avg = true,
    )

    println("Compute SurfaceTendency terms...")
    terms, bundle = SurfaceTendency.computeSurfaceTendencyTerms(;
        X_lc    = d_lc["Theta"],
#        X_c    = d_c["THETA"],
        h_l    = h_l,
        h_c    = h_c,
        h_r    = h_r,
        XUVEL   = d_c["ADVx_TH"],
        XVVEL   = d_c["ADVy_TH"],
        XWVEL   = d_c["ADVr_TH"],
        XDIFFFLX = d_c["ADVx_TH"] * 0.0,
        YDIFFFLX = d_c["ADVy_TH"] * 0.0,
        ZDIFFFLX = d_c["DFrI_TH"],
        Fsol   = d_c["oceQsw"],
        Fnet   = d_c["TFLUX"],
        CORRECTION_SFCFLX = d_c["WTHMASS"][:, :, 1],
        Δt     = diter * model_dt,
        coo    =  coo,
        coo_s  = coo_s,
        tracer = "TEMP",
    )

        
    terms["TEND_SUM"] = (
          terms["TEND_SW"]
        + terms["TEND_SFCFLX"]
        + terms["TEND_CORRFLX"]
        + terms["TEND_ENT_dhdt"] 
        + terms["TEND_TRANSPORT"] 
        + terms["TEND_HDIFF"] 
        + terms["TEND_VDIFF"]
    )



    #terms["TEND_RES_NOVDIFF"] = terms["TEND_SUM_NOVDIFF"] - TOTTTEND_mean
    terms["TEND_RES"] = terms["TEND_SUM"] - TOTTTEND_mean

    bundle["TEND_SUM"] = (
          bundle["TEND_SW"]
        + bundle["TEND_SFCFLX"]
        + bundle["TEND_CORRFLX"]
        + bundle["TEND_TRANSPORT"]
        + bundle["TEND_HDIFF"]
        + bundle["TEND_VDIFF"]
    )

    bundle_sum_mean = Operators_ML.computeMLMean(
        bundle["TEND_SUM"],
        coo,
        h = h_c,
        do_avg = true,
    )

    bundle["TEND_RES"] = bundle["TEND_SUM"] - d_c["TOTTTEND"]

    elm_type = eltype(d_c["THETA"])

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

        for varname in keys(bundle)
            varname_new = "$(varname)_3D"
            _var = defVar(ds, varname_new, elm_type, ("lon", "lat", "z", "time"), fillvalue=fill_value)
            _var[:, :, :, 1] = bundle[varname]
        end


        println(coo.gd.z_W[1, 1, :])
        for (varname, vardata, datatype, dimnames) in [
            ("TOTTTEND_mean",  TOTTTEND_mean,   elm_type, ("lon", "lat",)),
            ("THETA",          d_c["THETA"], elm_type, ("lon", "lat", "z",)),
            ("SALT",           d_c["SALT"], elm_type, ("lon", "lat", "z",)),
            ("h_c_algo",       h_c_algo, elm_type, ("lon", "lat", )),
            ("bundle_sum_mean",bundle_sum_mean, elm_type, ("lon", "lat", )),
            ("KPPhbl",         d_c["KPPhbl"], elm_type, ("lon", "lat", )),
            ("h_c",            h_c, elm_type, ("lon", "lat",)),
            ("oceFWflx",       d_c["oceFWflx"], elm_type, ("lon", "lat",)),
            ("EXFuwind",       d_c["EXFuwind"], elm_type, ("lon", "lat",)),
            ("EXFvwind",       d_c["EXFvwind"], elm_type, ("lon", "lat",)),
        ]

            colons = [Colon() for i=1:length(dimnames)]

            dimnames = [dimnames..., "time"]

            _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
            _var[colons..., 1] = vardata
        end

        for (varname, vardata, datatype, dimnames) in [
            ("z_w",            coo.gd.z_W[1, 1, :], elm_type, ("z_w", )),
            ("z",              coo.gd.z_T[1, 1, :], elm_type, ("z", )),
            ("lon",            coo.gd.λ_T[:, 1, 1], elm_type, ("lon", )),
            ("lat",            coo.gd.ϕ_T[1, :, 1], elm_type, ("lat", )),
        ]
            _var = defVar(ds, varname, datatype, dimnames; fillvalue=fill_value)
            _var[:] = vardata
        end

    end

end


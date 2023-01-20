if !(:SurfaceFluxes in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "SurfaceFluxes.jl")))
end

if !(:Operators in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "Operators.jl")))
end

if !(:Operators_ML in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "Operators_ML.jl")))
end

if !(:OceanConstants in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "OceanConstants.jl")))
end

if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end



module SurfaceTendency


    using ..Operators 
    using ..Operators_ML
    using ..SurfaceFluxes
    using ..OceanConstants
    using ..CoordinateModule


    function computeSurfaceTendencyTerms(;
        X_lc   :: AbstractArray{T, 3},  # T-grid X is any tracer
        h_l    :: AbstractArray{T, 2},  # T-grid  l = left   = past
        h_c    :: AbstractArray{T, 2},  # T-grid  c = center = now
        h_r    :: AbstractArray{T, 2},  # T-grid  r = right  = future
        XUVEL   :: AbstractArray{T, 3},  # U-grid
        XVVEL   :: AbstractArray{T, 3},  # V-grid
        XWVEL   :: AbstractArray{T, 3},  # W-grid
        XDIFFFLX :: AbstractArray{T, 3},  # U-grid
        YDIFFFLX :: AbstractArray{T, 3},  # V-grid
        ZDIFFFLX :: AbstractArray{T, 3},  # W-grid
        CORRECTION_SFCFLX :: AbstractArray{T, 2},
        Fsol   :: Union{AbstractArray{T, 2}, Nothing} = nothing,  # W-grid (Solar radiation input)
        Fnet   :: Union{AbstractArray{T, 2}, Nothing} = nothing,  # W-grid (top-grid heat fluxes)
        Δt     :: T,
        coo    :: CoordinateModule.Coordinate,
        coo_s  :: CoordinateModule.Coordinate,  # slab coo
        tracer :: String = "TEMP"
    ) where T

        gsp = coo.gsp
        gd = coo.gd

        Δz_T = view(gsp.Δz_T, 1, 1, :)
        z_W = view(coo.gd.z_W, 1, 1, :)
        mask_T = coo.gd.mask_T

        h_lc = (h_l + h_c) / 2
        h_cr = (h_c + h_r) / 2

        h_lc .= - coo.gd.z_W[:, :, 2]
        h_cr .= - coo.gd.z_W[:, :, 30]


        #X_lc = (X_l + X_c) / 2

        # clear them in case reused
        h_l = nothing
        h_c = nothing
        h_r = nothing


        result = Dict()

        bundle = Dict()
        SWFLX = Fsol
        SWFLX = - reshape(SWFLX, size(SWFLX)..., 1) .* SurfaceFluxes.RADFLX_shape.(coo.gd.z_W)
        bundle["TEND_SW"] = - Operators.T_DIVz_W(
            SWFLX,
            coo,
            already_weighted = false,
        ) / (OceanConstants.ρ * OceanConstants.c_p)


        SFCFLX = Fnet - Fsol
        SFCFLX = - reshape(SFCFLX, size(SFCFLX)..., 1) .* SurfaceFluxes.SFCFLX_shape.(coo.gd.z_W)
        bundle["TEND_SFCFLX"] = - Operators.T_DIVz_W(
            SFCFLX,
            coo,
            already_weighted = false,
        ) / (OceanConstants.ρ * OceanConstants.c_p) 

        CORRFLX = reshape(CORRECTION_SFCFLX, size(CORRECTION_SFCFLX)..., 1) .* SurfaceFluxes.SFCFLX_shape.(coo.gd.z_W)
        bundle["TEND_CORRFLX"] = - Operators.T_DIVz_W(
            CORRFLX,
            coo,
            already_weighted = false,
        )
        
        bundle["TEND_HDIFF"] = - (
              Operators.T_DIVx_U(XDIFFFLX, coo; already_weighted=true)
            + Operators.T_DIVy_V(YDIFFFLX, coo; already_weighted=true)
        )

        bundle["TEND_VDIFF"] = - Operators.T_DIVz_W(
            ZDIFFFLX,
            coo,
            already_weighted = true,
        )

        bundle["TEND_TRANSPORT"] = - (
              Operators.T_DIVx_U(XUVEL, coo; already_weighted=true)
            + Operators.T_DIVy_V(XVVEL, coo; already_weighted=true)
            + Operators.T_DIVz_W(XWVEL, coo; already_weighted=true)
        )


        # Entrainment dhdt is special

        X_mean_lc = Operators_ML.computeMLMean(
            X_lc,
            coo,
            h = h_lc,
            do_avg = true,
        )

        X_mean_cr = Operators_ML.computeMLMean(
            X_lc,
            coo,
            h = h_cr,
            do_avg = true,
        )

        result["TEND_ENT_dhdt"] = ( X_mean_cr - X_mean_lc ) / Δt

        result["TEND_SW"] = Operators_ML.computeMLMean(
            bundle["TEND_SW"],
            coo,
            h = h_cr,
            do_avg = true,
        )

        result["TEND_SFCFLX"] = Operators_ML.computeMLMean(
            bundle["TEND_SFCFLX"],
            coo,
            h = h_cr,
            do_avg = true,
        )

        result["TEND_TRANSPORT"] = Operators_ML.computeMLMean(
            bundle["TEND_TRANSPORT"],
            coo,
            h = h_cr,
            do_avg = true,
        )
       
        result["TEND_HDIFF"] = Operators_ML.computeMLMean(
            bundle["TEND_HDIFF"],
            coo,
            h = h_cr,
            do_avg = true,
        )
 
        result["TEND_VDIFF"] = Operators_ML.computeMLMean(
            bundle["TEND_VDIFF"],
            coo,
            h = h_cr,
            do_avg = true,
        )

        result["TEND_CORRFLX"] = Operators_ML.computeMLMean(
            bundle["TEND_CORRFLX"],
            coo,
            h = h_cr,
            do_avg = true,
        )



        return result, bundle
    end   





end

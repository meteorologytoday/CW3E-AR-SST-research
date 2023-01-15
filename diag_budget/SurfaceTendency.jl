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
        X      :: AbstractArray{T, 3},  # T-grid X is any tracer
        h_l    :: AbstractArray{T, 2},  # T-grid  l = left   = past
        h_c    :: AbstractArray{T, 2},  # T-grid  c = center = now
        h_r    :: AbstractArray{T, 2},  # T-grid  r = right  = future
        UVEL   :: AbstractArray{T, 3},  # U-grid
        VVEL   :: AbstractArray{T, 3},  # V-grid
        WVEL   :: AbstractArray{T, 3},  # W-grid
        VDIFFFLX :: AbstractArray{T, 3},  # W-grid
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

        dhdt = (h_r - h_l) / (2*Δt) 

        result = Dict()

        # First term : Surface flux
       
        if tracer == "TEMP"
            F0 = Fnet
            Fb = Fsol .* SurfaceFluxes.RADFLX_shape.(- h_c)
            
            result["TEND_SFCFLX"] = - (F0 - Fb) ./ (OceanConstants.ρ * OceanConstants.c_p * h_c)
        else
            throw(ErrorException("Unknown tracer: $tracer"))
        end    



        # Preparation for entrainment 
        X_mean = Operators_ML.sT_ML∫dz_T(
            X,
            coo,
            h = h_c,
            do_avg = true,
        )

        X_b = Operators_ML.evalAtMLD_T(
            X,
            coo,
            h = h_c,
        )

        ΔX = X_mean - X_b

        # Preparation for velocity

        UVEL = Operators.T_interp_U( UVEL, coo)
        VVEL = Operators.T_interp_V( VVEL, coo)


        UVEL_mean = Operators_ML.sT_ML∫dz_T(
            UVEL,
            coo,
            h = h_c,
            do_avg = true,
        )
        
        VVEL_mean = Operators_ML.sT_ML∫dz_T(
            VVEL,
            coo,
            h = h_c,
            do_avg = true,
        )



        # Second term : Entrainment dhdt
        result["TEND_ENT_dhdt"] = - dhdt ./ h_c .* ΔX


        # Thrid term : Entrainment wb
        wb = Operators_ML.evalAtMLD_W(
            WVEL,
            coo,
            h = h_c,
        )
        result["TEND_ENT_wb"] = - wb ./ h_c .* ΔX
        
        # Fourth term : Entrainment due to horizontal adv
      
        reshaped_h_c = reshape(h_c, size(h_c)..., 1) 

        Udhdx = UVEL_mean .* Operators.T_interp_U( Operators.U_∂x_T(reshaped_h_c, coo_s), coo_s)[:, :, 1]
        Vdhdy = VVEL_mean .* Operators.T_interp_V( Operators.V_∂y_T(reshaped_h_c, coo_s), coo_s)[:, :, 1]

        hadv = - (Udhdx + Vdhdy) 
        
        result["TEND_ENT_hadv"] = hadv ./ h_c .* ΔX


        # Fifth term : barotropic advection 

        dX_meandx = UVEL_mean .* Operators.T_interp_U( Operators.U_∂x_T(X, coo_s), coo_s)[:, :, 1]
        dX_meandy = VVEL_mean .* Operators.T_interp_V( Operators.V_∂y_T(X, coo_s), coo_s)[:, :, 1]
        result["TEND_BT_hadv"] = - (dX_meandx + dX_meandy)
 
        # Eighth term : Entrainment wb
        result["TEND_VDIFF"] = Operators_ML.evalAtMLD_W(
            VDIFFFLX,
            coo,
            h = h_c,
        ) ./ h_c
         
        return result
    end   





end

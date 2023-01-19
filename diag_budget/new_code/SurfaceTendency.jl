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
        X_p    :: AbstractArray{T, 3},  # T-grid X of the past
        X_n    :: AbstractArray{T, 3},  # T-grid X of now
        h_p    :: AbstractArray{T, 2},  # T-grid  p = past
        h_n    :: AbstractArray{T, 2},  # T-grid  n = now
        U      :: AbstractArray{T, 3},  # U-grid
        V      :: AbstractArray{T, 3},  # V-grid
        W      :: AbstractArray{T, 3},  # W-grid
        UX     :: AbstractArray{T, 3},  # U-grid
        VX     :: AbstractArray{T, 3},  # V-grid
        WX     :: AbstractArray{T, 3},  # W-grid
        XDIFFFLX :: AbstractArray{T, 3},  # U-grid
        YDIFFFLX :: AbstractArray{T, 3},  # V-grid
        ZDIFFFLX :: AbstractArray{T, 3},  # W-grid
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

        dhdt = (h_n - h_p) / Δt 

        result = Dict()

        # First term : Surface flux
       
        if tracer == "TEMP"
            F0 = Fnet
            Fb = Fsol .* SurfaceFluxes.RADFLX_shape.(- h_n)
            
            result["TEND_SFCFLX"] = (F0 - Fb) ./ (OceanConstants.ρ * OceanConstants.c_p * h_n)
        else
            throw(ErrorException("Unknown tracer: $tracer"))
        end    


        # Preparation for entrainment 
        X_mean_n = Operators_ML.computeMLMean(
            X,
            coo,
            h = h_n,
            do_avg = true,
            keep_dim = true,
        )
        X_prime_n = X .- X_mean

        X_b = Operators_ML.evalAtMLD_T(
            X,
            coo,
            h = h_c,
        )

        ΔX = X_mean - X_b

        # Preparation for velocity

        h_c_V = Operators.V_interp_T(reshape(h_c, size(h_c)..., 1), coo_s)[:, :, 1]
        h_c_U = Operators.U_interp_T(reshape(h_c, size(h_c)..., 1), coo_s)[:, :, 1]

        UVEL_mean = Operators_ML.computeMLMean(
            UVEL,
            coo,
            h = h_c_U,
            do_avg = true,
            keep_dim = true,
        )
        
        VVEL_mean = Operators_ML.computeMLMean(
            VVEL,
            coo,
            h = h_c_V,
            do_avg = true,
            keep_dim = true,
        )

        
        UVEL_prime = UVEL .- UVEL_mean
        VVEL_prime = VVEL .- VVEL_mean


        
        # Second term : Entrainment dhdt
        X_mean_p = Operators_ML.computeMLMean(
            X_p,
            coo,
            h = h_p,
            do_avg = true,
        )

        X_mean_n = Operators_ML.computeMLMean(
            X_p,
            coo,
            h = h_n,
            do_avg = true,
        )

        result["TEND_ENT_dhdt"] = ( X_mean_n - X_mean_p ) / Δt
        
        
        # Thrid term : Entrainment wb
        wb = Operators_ML.evalAtMLD_W(
            WVEL,
            coo,
            h = h_c,
        )
        result["TEND_ENT_wb"] = - wb ./ h_c .* ΔX
        
        # Fourth term : Entrainment due to horizontal adv
      
        reshaped_h_c = reshape(h_c, size(h_c)..., 1) 

        #Udhdx = UVEL_mean .* Operators.T_interp_U( Operators.U_∂x_T(reshaped_h_c, coo_s), coo_s)[:, :, 1]
        #Vdhdy = VVEL_mean .* Operators.T_interp_V( Operators.V_∂y_T(reshaped_h_c, coo_s), coo_s)[:, :, 1]

        Udhdx = view(Operators.T_interp_U( UVEL_mean .* Operators.U_∂x_T(reshaped_h_c, coo_s), coo_s), :, :, 1)
        Vdhdy = view(Operators.T_interp_V( VVEL_mean .* Operators.V_∂y_T(reshaped_h_c, coo_s), coo_s), :, :, 1)

        hadv = - (Udhdx + Vdhdy) 
        
        result["TEND_ENT_hadv"] = hadv ./ h_c .* ΔX


        # Fifth term : barotropic advection 
        dX_meandx = view(Operators.T_interp_U( UVEL_mean .* Operators.U_∂x_T(X, coo_s), coo_s), :, :, 1)
        dX_meandy = view(Operators.T_interp_V( VVEL_mean .* Operators.V_∂y_T(X, coo_s), coo_s), :, :, 1)
        result["TEND_BT_hadv"] = - (dX_meandx + dX_meandy)

        # Sixth term: eddy

        XU_eddy = Operators.T_DIVx_U( Operators.U_interp_T(X_prime, coo) .* UVEL_prime , coo ; already_weighted=false)
        XV_eddy = Operators.T_DIVy_V( Operators.V_interp_T(X_prime, coo) .* VVEL_prime , coo ; already_weighted=false)
        X_eddy = - ( XU_eddy + XV_eddy )       
        
        result["TEND_EDDY"] = Operators_ML.computeMLMean(
            X_eddy,
            coo,
            h = h_c,
            do_avg = true,
        )

        # Seventh term : Horizontal diffusion

        HDIV = (
              Operators.T_DIVx_U(XDIFFFLX, coo; already_weighted=true)
            + Operators.T_DIVy_V(YDIFFFLX, coo; already_weighted=true)
        )

        result["TEND_HDIFF"] = Operators_ML.computeMLMean(
            HDIV,
            coo,
            h = h_c,
            do_avg = true,
        )
 
        # Eighth term : Entrainment wb
        result["TEND_VDIFF"] = Operators_ML.evalAtMLD_W(
            ZDIFFFLX,
            coo,
            h = h_c,
        ) ./ coo_s.gsp.Δa_W[:, :, 1] ./ h_c
         
        return result
    end   





end

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


module SurfaceTendency


    using ..Operators 
    using ..Operators_ML
    using ..SurfaceFluxes
    using ..OceanConstants


    function computeSurfaceTendencyTerms(;
        X      :: AbstractArray{T, 3},  # T-grid X is any tracer
        TEMP_l :: AbstractArray{T, 3},  # T-grid  l = left   = past
        SALT_l :: AbstractArray{T, 3},  # T-grid  
        TEMP   :: AbstractArray{T, 3},  # T-grid  c = center = now
        SALT   :: AbstractArray{T, 3},  # T-grid  
        TEMP_r :: AbstractArray{T, 3},  # T-grid  r = right  = future
        SALT_r :: AbstractArray{T, 3},  # T-grid
        VVEL   :: AbstractArray{T, 3},  # V-grid
        UVEL   :: AbstractArray{T, 3},  # U-grid
        WVEL   :: AbstractArray{T, 3},  # W-grid
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


        # Determine MLD
        h_l = detectMLD(TEMP_l, SALT_l, coo)
        h_c = detectMLD(TEMP_c, SALT_c, coo)
        h_r = detectMLD(TEMP_r, SALT_r, coo)

        dhdt = (h_r - h_l) / (2*Δt) 

        # First term : Surface flux
       
        if tracer == "TEMP"
            F0 = Fnet
            Fb = Fsol .* RADFLX_shape.(- h_c)
            
            TEND_SFCFLX = - (F0 - Fb) ./ (OceanConstants.ρ * OceanConstants.c_p * h_c)
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


        # Second term : Entrainment dhdt
        TEND_ENT_dhdt = - dhdt ./ h_c .* ΔX


        # Thrid term : Entrainment wb
        wb = Operators_ML.evalAtMLD_W(
            WVEL,
            coo,
            h = h_c,
        )
        TEND_ENT_wb = - wb ./ h_c .* ΔX
        
        # Fourth term : Entrainment due to horizontal adv
        
        dhdx = Operators.T_interp_U(Operators.U_∂x_T(h, coo_s), coo_s)
        dhdy = Operators.T_interp_V(Operators.V_∂y_T(h, coo_s), coo_s)
        hadv = - (UVEL_mean .* dhdx + VVEL_mean .* dhdy)
        
        TEND_ENT_hadv = hadv ./ h_c .* ΔX


        # 
         
        
        
        
    end   





end

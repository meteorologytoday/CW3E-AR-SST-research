if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end

if !(:Operators in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "Operators.jl")))
end


if !(:BuoyancyNonlinear in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "BuoyancyNonlinear.jl")))
end



module Operators_ML

    using ..CoordinateModule
    using ..Operators
    using ..BuoyancyNonlinear


    fill_value = 1e20

    function sT_ML∫dz_T(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        h    :: AbstractArray{T, 2},
        Nz_h :: Union{AbstractArray{T, 2}, Nothing} = nothing,
        do_avg :: Bool = false,
    ) where T

        gsp = coo.gsp
        gd = coo.gd
        Δz_T = view(gsp.Δz_T, 1, 1, :)
        z_W = view(coo.gd.z_W, 1, 1, :)
        mask_T = coo.gd.mask_T

        fo = zeros(T, gsp.Nx, gsp.Ny)
        
        if Nz_h == nothing
            Nz_h = detectMLNz(h, coo)
        end


        for j=1:gsp.Ny, i=1:gsp.Nx
            
            _tmp = 0.0
            
            if mask_T[i, j, 1] == 0
                fo[i, j] = fill_value
            else
                _Nz = Nz_h[i, j]
                for k=1:_Nz-1
                    _tmp += Δz_T[k] * fi[i, j, k]
                end

                _tmp += (z_W[_Nz] + h[i, j]) * fi[i, j, _Nz]
                
            end

            fo[i, j] = _tmp
        end

        if do_avg
            for j=1:gsp.Ny, i=1:gsp.Nx
                if mask_T[i, j, 1] != 0
                    fo[i, j] /= h[i, j]
                end
            end
        end

        
        return fo
    end

    function detectMLD(
        TEMP :: AbstractArray{T, 3},
        SALT :: AbstractArray{T, 3},
        coo  :: CoordinateModule.Coordinate;
        algo :: String = "ρ_threshold",
        algo_params :: Dict = Dict("ρ_threshold" => 0.03),
    ) where T

        gd = coo.gd
        mask = view(gd.mask_T, :, :, 1)

        bundle = Dict()
        MLD    = zeros(T,     gd.Nx, gd.Ny)
        z_T = view(gd.z_T, 1, 1, :)
       

        MLD[mask .== 0] .= NaN

        if algo == "ρ_threshold"

            ρ_threshold = algo_params["ρ_threshold"]
            ρ = zeros(T, gd.Nx, gd.Ny, gd.Nz)
            BuoyancyNonlinear.TS2ρ!(TEMP, SALT, ρ)
            
            bundle["ρ"] = ρ

            for j=1:gd.Ny, i=1:gd.Nx
                
                if isnan(MLD[i, j])
                    continue
                end
                 
                _ρ = view(ρ, i, j, :)
    
                δρ = _ρ .- _ρ[1]
                
                for k=2:gd.Nz
                    
                    if δρ[k] > ρ_threshold # since δρ[1] = 0 by def, no need to test δρ[k-1]
                        MLD[i, j] = - z_T[k-1] + (z_T[k-1] - z_T[k]) * (ρ_threshold - δρ[k-1])/(δρ[k] - δρ[k-1]) 
                        #MLD[i, j] = - z_T[k]
                        break
                    end

                    if k == gd.Nz  # cannot find the mixed layer depth
                        MLD[i, j] = - z_T[k]
                    end
                    
                end
                
            end
        
        else
            throw(ErrorException("Unknown algorithm: $algo"))
        end

        # Sanity check
        if any(MLD[isfinite.(MLD)] .<= 0)
            throw(ErrorException("Some MLD is negative."))
        end

        return MLD, bundle
        
    end


    function detectMLNz(
        h   :: AbstractArray{T, 2},
        coo :: CoordinateModule.Coordinate,
    ) where T

        gd = coo.gd

        MLNz = zeros(Int64, gd.Nx, gd.Ny)
        z_W = view(gd.z_W, 1, 1, :)
       
        for j=1:gd.Ny, i=1:gd.Nx
            
            z = - h[i, j]
            for k=1:gd.Nz
                if z_W[k] >= z >= z_W[k+1]   # Using 2 equalities so that the last grid-box will include z = z_bottom
                    MLNz[i, j] = k
                    break
                end
            end
            
        end
        
        return MLNz
        
    end
    function evalAtMLD_W(
        fi  :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        h    :: AbstractArray{T, 2},
        Nz_h :: Union{AbstractArray{T, 2}, Nothing} = nothing,
        mask_T :: Union{AbstractArray, Nothing} = nothing,
    ) where T

        gd = coo.gd
        z_W = view(gd.z_W, 1, 1, :)
        Δz_T = view(coo.gsp.Δz_T, 1, 1, :)
        mask_T = gd.mask_T

        fo = zeros(T, gd.Nx, gd.Ny)

        if Nz_h == nothing
            Nz_h = detectMLNz(h, coo)
        end

        for j=1:gd.Ny, i=1:gd.Nx

            if mask_T[i, j] == 0
                fo[i, j] = fill_value
            else
                _Nz = Nz_h[i, j]
                _h  = h[i, j]
                fo[i, j] = fi[i, j, _Nz] + (fi[i, j, _Nz+1] - fi[i, j, _Nz]) * (z_W[_Nz] + _h) / Δz_T[_Nz]
            end        
        end
        
        return fo
        
    end


    function evalAtMLD_T(
        fi  :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        h    :: AbstractArray{T, 2},
        Nz_h :: Union{AbstractArray{T, 2}, Nothing} = nothing,
    ) where T

        gd = coo.gd
        z_W = view(gd.z_W, 1, 1, :)
        Δz_T = view(coo.gsp.Δz_T, 1, 1, :)
        mask_T = gd.mask_T

        fo = zeros(T, gd.Nx, gd.Ny)

        if Nz_h == nothing
            Nz_h = detectMLNz(h, coo)
        end

        for j=1:gd.Ny, i=1:gd.Nx

            if mask_T[i, j] == 0
                fo[i, j] = fill_value
            else
                _Nz = Nz_h[i, j]
                fo[i, j] = fi[i, j, _Nz]
            end        
        end
        
        return fo
        
    end


end

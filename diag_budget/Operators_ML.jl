

module Operators_ML

    using ..CoordinateModule

    function sT_ML∫dz_T(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        h    :: AbstractArray{T, 2},
        Nz_h :: Union{AbstractArray{T, 2}, Nothing} = nothing,
    ) where T

        gsp = coo.gsp

        fo = zeros(T, gsp.Nx, gsp.Ny, 1)
        
        if Nz_h == nothing
            Nz_h = detectMLNz(h, coo)
        end


        for k=1:gsp.Nz, j=1:gsp.Ny, i=1:gsp.Nx
            if mask_T[i, j, 1] == 0
                fo[i, j, k] = fill_value
            end
        end

        return fo
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


end

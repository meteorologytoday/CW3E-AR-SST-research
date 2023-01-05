

module Operators

    using ..CoordinateModule

    function T_mask_T(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        fill_value :: T = NaN,
    ) where T

        gsp = coo.gsp
        mask_T = coo.gd.mask_T

        fo = zeros(T, gsp.Nx, gsp.Ny, gsp.Nz)
        fo .= fi
        for k=1:gsp.Nz, j=1:gsp.Ny, i=1:gsp.Nx
            if mask_T[i, j, 1] == 0
                fo[i, j, k] = fill_value
            end
        end

        return fo
    end


    function T_interp_U(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
    ) where T

        gsp = coo.gsp

        fo = zeros(T, gsp.Nx, gsp.Ny, gsp.Nz)
        for k=1:gsp.Nz, j=1:gsp.Ny, i=1:gsp.Nx
            fo[i, j, k] = (fi[i, j, k] + fi[i+1, j, k]) / 2
        end

        return fo
    end


    function U_∂x_T(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate,
    ) where T

        gsp = coo.gsp

        fo = zeros(T, gsp.Nx+1, gsp.Ny, gsp.Nz)
        for k=1:gsp.Nz, j=1:gsp.Ny, i=2:gsp.Nx
            fo[i, j, k] = (fi[i, j, k] - fi[i-1, j, k]) / gsp.Δx_U[i, j, 1]
        end

        return fo
    end


end

if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end

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

    function T_interp_V(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
    ) where T

        gsp = coo.gsp

        fo = zeros(T, gsp.Nx, gsp.Ny, gsp.Nz)
        for k=1:gsp.Nz, j=1:gsp.Ny, i=1:gsp.Nx
            fo[i, j, k] = (fi[i, j, k] + fi[i, j+1, k]) / 2
        end

        return fo
    end


    function V_∂y_T(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate,
    ) where T

        gsp = coo.gsp

        fo = zeros(T, gsp.Nx, gsp.Ny+1, gsp.Nz)
        for k=1:gsp.Nz, j=2:gsp.Ny, i=1:gsp.Nx
            fo[i, j, k] = (fi[i, j, k] - fi[i, j-1, k]) / gsp.Δy_V[i, j, 1]
        end

        return fo
    end

    function T_DIVy_V(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        already_weighted :: Bool = true
    ) where T

        gsp = coo.gsp
        Δx_V = gsp.Δx_V
        Δa_T = gsp.Δa_T
        Δz_T = gsp.Δz_T
        Nx, Ny, Nz = gsp.Nx, gsp.Ny, gsp.Nz


        fo = zeros(T, Nx, Ny, Nz)

        if already_weighted
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i, j+1, k] - fi[i, j, k]) / (Δa_T[i, j, 1] * Δz_T[1, 1, k])
            end
        else
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i, j+1, k] * Δx_V[i, j+1, 1] - fi[i, j, k] * Δx_V[i, j, 1]) / Δa_T[i, j, 1]
            end
        end

        return fo
    end

    function T_DIVx_U(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        already_weighted :: Bool = true,
    ) where T

        gsp = coo.gsp
        Δy_U = gsp.Δy_U
        Δa_T = gsp.Δa_T
        Δz_T = gsp.Δz_T
        Nx, Ny, Nz = gsp.Nx, gsp.Ny, gsp.Nz


        fo = zeros(T, Nx, Ny, Nz)

        if already_weighted
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i+1, j, k] - fi[i, j, k]) / (Δa_T[i, j, 1] * Δz_T[1, 1, k])
            end
        else
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i+1, j, k] * Δy_U[i+1, j, 1] - fi[i, j, k] * Δy_U[i, j, 1]) / Δa_T[i, j, 1]
            end
        end


        return fo
    end

    function T_DIVz_W(
        fi :: AbstractArray{T, 3},
        coo :: CoordinateModule.Coordinate;
        already_weighted :: Bool = true,
    ) where T

        gsp = coo.gsp
        Δa_T = gsp.Δa_T
        Δz_T = gsp.Δz_T
        Nx, Ny, Nz = gsp.Nx, gsp.Ny, gsp.Nz


        fo = zeros(T, Nx, Ny, Nz)

        if already_weighted
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i, j, k] - fi[i, j, k+1]) / (Δa_T[i, j, 1] * Δz_T[1, 1, k])
            end
        else
            for k=1:Nz, j=1:Ny, i=1:Nx
                fo[i, j, k] = (fi[i, j, k] - fi[i, j, k+1]) / Δz_T[1, 1, k]
            end
        end


        return fo
    end




end

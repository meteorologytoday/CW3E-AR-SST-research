if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end

module TestTools

    using Formatting
    using ..CoordinateModule

    function genGrid(
        raw_λ_U :: AbstractArray{T, 1},
        raw_ϕ_V :: AbstractArray{T, 1},
        raw_z_W :: AbstractArray{T, 1};
        R :: Float64 = 6.371e6,
    ) where T

        Nx = length(raw_λ_U) - 1
        Ny = length(raw_ϕ_V) - 1
        Nz = length(raw_z_W) - 1

        println(format("(Nx, Ny, Nz) = ({:d}, {:d}, {:d})", Nx, Ny, Nz))

        # ===== [ λ ] ===== 

        λ_U = repeat(reshape(raw_λ_U, :, 1, 1), outer=(1, Ny, 1))
        λ_T  = (λ_U[1:end-1, :, :] + λ_U[2:end, :, :]) / 2
        λ_V  = repeat(λ_T[:, 1:1, :], outer=(1, Ny+1, 1))
        λ_UV = repeat(reshape(raw_λ_U, 1, :, 1), outer=(1, Ny+1, 1))

        λ_W = copy(λ_T) 

        # ===== [ ϕ ] ===== 

        ϕ_V = repeat(reshape(raw_ϕ_V, 1, :, 1), outer=(Nx, 1, 1))
        ϕ_T  = (ϕ_V[:, 1:end-1, :] + ϕ_V[:, 2:end, :]) / 2
        ϕ_U  = repeat(ϕ_T[1:1, :, :], outer=(Nx+1, 1, 1))
        ϕ_UV = repeat(reshape(raw_ϕ_V, 1, :, 1), outer=(Nx+1, 1, 1))

        ϕ_W = copy(ϕ_T)
     
        # ===== [ z ] ===== 

        z_W = reshape(copy(raw_z_W), 1, 1, :)

        z_T = ( z_W[:, :, 1:end-1] + z_W[:, :, 2:end] ) / 2
        z_U  = copy(z_T)
        z_V  = copy(z_T)
        z_UV = copy(z_T)

        # ===== [ mask ] ===== 

        mask_T = ones(T, Nx, Ny, 1)

        # ===== [ Δx ] ===== 

        raw_Δλ_T = raw_λ_U[2:end] - raw_λ_U[1:end-1]
        raw_Δλ_U = zeros(T, Nx+1)
        raw_Δλ_U[2:end-1] = (raw_Δλ_T[1:end-1] + raw_Δλ_T[2:end]) / 2
        raw_Δλ_U[1]   = raw_Δλ_U[2]
        raw_Δλ_U[end] = raw_Δλ_U[end-1]
        
        Δλ_T   = repeat( reshape(raw_Δλ_T, :, 1, 1), outer=(1, Ny,   1) )
        Δλ_U   = repeat( reshape(raw_Δλ_U, :, 1, 1), outer=(1, Ny,   1) )
        Δλ_V   = repeat( reshape(raw_Δλ_T, :, 1, 1), outer=(1, Ny+1, 1) )
        Δλ_UV  = repeat( reshape(raw_Δλ_U, :, 1, 1), outer=(1, Ny+1, 1) )

        Δx_T   = R * cos.(deg2rad.(ϕ_T))  .* deg2rad.(Δλ_T)
        Δx_U   = R * cos.(deg2rad.(ϕ_U))  .* deg2rad.(Δλ_U)
        Δx_V   = R * cos.(deg2rad.(ϕ_V))  .* deg2rad.(Δλ_V)
        Δx_UV  = R * cos.(deg2rad.(ϕ_UV)) .* deg2rad.(Δλ_UV)
        Δx_W   = copy(Δx_T)
        
        # ===== [ Δy ] ===== 

        raw_Δϕ_T = raw_ϕ_V[2:end] - raw_ϕ_V[1:end-1]
        raw_Δϕ_V = zeros(T, Ny+1)
        raw_Δϕ_V[2:end-1] = (raw_Δϕ_T[1:end-1] + raw_Δϕ_T[2:end]) / 2
        raw_Δϕ_V[1]   = raw_Δϕ_V[2]
        raw_Δϕ_V[end] = raw_Δϕ_V[end-1]
 
        Δϕ_T   = repeat( reshape(raw_Δϕ_T, 1, :, 1), outer=(Nx,   1,   1) )
        Δϕ_U   = repeat( reshape(raw_Δϕ_T, 1, :, 1), outer=(Nx+1, 1,   1) )
        Δϕ_V   = repeat( reshape(raw_Δϕ_V, 1, :, 1), outer=(Nx,   1,   1) )
        Δϕ_UV  = repeat( reshape(raw_Δϕ_V, 1, :, 1), outer=(Nx+1, 1,   1) )

        Δy_T  = R * deg2rad.(Δϕ_T) 
        Δy_U  = R * deg2rad.(Δϕ_T) 
        Δy_V  = R * deg2rad.(Δϕ_T) 
        Δy_UV = R * deg2rad.(Δϕ_T) 
        Δy_W  = copy(Δy_T)

        # ===== [ Δz ] ===== 

        Δz_T  = z_W[:, :, 1:end-1] - z_W[:, :, 2:end]
        Δz_U  = copy(Δz_T)
        Δz_V  = copy(Δz_T)
        Δz_UV = copy(Δz_T)
        Δz_W  = zeros(T, 1, 1, Nz+1)

        Δz_W[:, :, 2:end-1] = (Δz_T[:, :, 1:end-1] + Δz_T[:, :, 2:end]) / 2
        Δz_W[:, :, 1]   = Δz_T[:, :, 1]
        Δz_W[:, :, end] = Δz_T[:, :, end]
        
        # ===== [ Δa ] =====
        Δa_T  = Δx_T .* Δy_T

        gd = CoordinateModule.Grid(

            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
        
            R = 6.371e6,
            
            λ_T  = λ_T,
            λ_U  = λ_U,
            λ_V  = λ_V,
            λ_W  = λ_W,
            λ_UV = λ_UV,

            ϕ_T  = ϕ_T,
            ϕ_U  = ϕ_U,
            ϕ_V  = ϕ_V,
            ϕ_W  = ϕ_W,
            ϕ_UV = ϕ_UV,

            z_T  = z_T,
            z_U  = z_U,
            z_V  = z_V,
            z_W  = z_W,
            z_UV = z_UV,

            mask_T = mask_T,
        )

        gsp = CoordinateModule.GridSpacing(

            Nx = Nx,
            Ny = Ny,
            Nz = Nz,

            Δx_T  = Δx_T,
            Δx_U  = Δx_U,
            Δx_V  = Δx_V,
            Δx_W  = Δx_W,
            Δx_UV = Δx_UV,

            Δy_T  = Δy_T,
            Δy_U  = Δy_U,
            Δy_V  = Δy_V,
            Δy_W  = Δy_W,
            Δy_UV = Δy_UV,

            Δz_T  = Δz_T,
            Δz_U  = Δz_U,
            Δz_V  = Δz_V,
            Δz_W  = Δz_W,
            Δz_UV = Δz_UV,

            Δa_T  = Δa_T,
        )

        return CoordinateModule.Coordinate(gd, gsp)
        

    end


    #=
            




    function loadMITGCM(
        input_dir :: String = ".";
        XG :: String = "XG",
        YG :: String = "YG",
        RF :: String = "RF",
    )

        m = pyimport("MITgcmutils")


        c = Dict()

        for varname, loaded_varname in Dict(
            :XG => XG,
            :YG => YG,
            :RF => RF,
        )
            c[varname] = m.mds.rdmds(joinpath(input_dir, varname))[:]
        end

        z_W = reshape(c[:RF], 1, 1, :)

        λ_U = zeros(eltype(c[:XG]), length(c[:XG])+1, 1, 1)
        λ_U[1:end-1, 1, 1] = c[:XG]
        λ_U[end, 1, 1] = c[:XG][end] + (c[:XG][end] - c[:XG][end-1])
     
        ϕ_V = zeros(eltype(c[:YG]), 1, length(c[:YG])+1, 1)
        ϕ_V[1, 1:end-1, 1] = c[:YG]
        ϕ_V[1,     end, 1] = c[:YG][end] + (c[:YG][end] - c[:YG][end-1])


    end

=#

end

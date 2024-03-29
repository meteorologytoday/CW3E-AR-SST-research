

Base.@kwdef struct GridSpacing

        Nx    :: Integer
        Ny    :: Integer
        Nz    :: Integer

        Δx_T  :: AbstractArray{Float64, 3}
        Δx_U  :: AbstractArray{Float64, 3}
        Δx_V  :: AbstractArray{Float64, 3}
        Δx_W  :: AbstractArray{Float64, 3}
        Δx_UV  :: AbstractArray{Float64, 3}

        Δy_T  :: AbstractArray{Float64, 3}
        Δy_U  :: AbstractArray{Float64, 3}
        Δy_V  :: AbstractArray{Float64, 3}
        Δy_W  :: AbstractArray{Float64, 3}
        Δy_UV  :: AbstractArray{Float64, 3}

        Δz_T  :: AbstractArray{Float64, 3}
        Δz_U  :: AbstractArray{Float64, 3}
        Δz_V  :: AbstractArray{Float64, 3}
        Δz_W  :: AbstractArray{Float64, 3}
        Δz_UV  :: AbstractArray{Float64, 3}

        Δa_T  :: AbstractArray{Float64, 3}
        Δa_W  :: AbstractArray{Float64, 3}
        
        Δvol_T  :: AbstractArray{Float64, 3}

end

#=


function genGrid_simple(;
    λ_U :: AbstractArray{Float64, 3},
    ϕ_V :: AbstractArray{Float64, 3},
    z_W :: AbstractArray{Float64, 3};
)

    Δλ_T = λ_U[2:end] - λ_U[1:end-1] 
    if sum(Δλ_T < 0) > 1
        throw(ErrorException("Only one jump in longitude is allowed."))
    end

    Δϕ_T = ϕ_V[2:end] - ϕ_V[1:end-1] 
    if any(Δϕ_T < 0)
        throw(ErrorException("Latitude should be monotonically increasing."))
    end

    Δz_T = z_W[1:end-1] - z_W[2:end] 
    if any(Δz_T < 0)
        throw(ErrorException("Z should be monotonically decreasing."))
    end

    jumped = false
    for i = 1:length(Δλ)
        if Δλ[i] < 0  # Cross the line
            λ_U[i+1:end] += 360.0
            break
        end
    end
   
    
end
=#
#=
        function


            vs_λ = copy(vs_λ)
            vs_ϕ = copy(vs_ϕ)
           
            if angle_unit == :deg

                λ_t .*= d2r
                ϕ_t .*= d2r

                vs_λ .*= d2r 
                vs_ϕ .*= d2r
     
            elseif angle_unit == :rad
                
                # do nothing

            else
                throw(ErrorException("Unknown `angle_unit`: " * angle_unit))
            end

            #println("Compute vectors and sides")
            # First, compute the needed vectors and sides

            permuted_vs_λ = permutedims(vs_λ, [2,3,1])
            permuted_vs_ϕ = permutedims(vs_ϕ, [2,3,1])
             
            shp = size(permuted_vs_ϕ)
            xv = reshape( R * cos.(permuted_vs_ϕ) .* cos.(permuted_vs_λ) , 1, shp...)
            yv = reshape( R * cos.(permuted_vs_ϕ) .* sin.(permuted_vs_λ) , 1, shp...)
            zv = reshape( R * sin.(permuted_vs_ϕ)                        , 1, shp...)

            ps = vcat(xv, yv, zv)
            u1 = ps[:, :, :, 2] - ps[:, :, :, 1]
            u2 = ps[:, :, :, 3] - ps[:, :, :, 2]
            u3 = ps[:, :, :, 4] - ps[:, :, :, 3]
            u4 = ps[:, :, :, 1] - ps[:, :, :, 4]

            ds1 = mapslices(norm, u1; dims=[1,])[1, :, :]
            ds2 = mapslices(norm, u2; dims=[1,])[1, :, :]
            ds3 = mapslices(norm, u3; dims=[1,])[1, :, :]
            ds4 = mapslices(norm, u4; dims=[1,])[1, :, :]

            #println("size of ds1: ", size(ds1))

            #=
            @time for j = 1:Ny
            for i = 1:Nx
                for k = 1:4
                    _λ = vs_λ[k, i, j]
                    _ϕ = vs_ϕ[k, i, j]
                    ps[1, k] = R * cos(_ϕ) * cos(_λ)
                    ps[2, k] = R * cos(_ϕ) * sin(_λ)
                    ps[3, k] = R * sin(_ϕ)
                end
                u1[:, i, j] = ps[:, 2] - ps[:, 1]
                u2[:, i, j] = ps[:, 3] - ps[:, 2]
                u3[:, i, j] = ps[:, 4] - ps[:, 3]
                u4[:, i, j] = ps[:, 1] - ps[:, 4]
                ds1[i, j] = norm(u1)
                ds2[i, j] = norm(u2)
                ds3[i, j] = norm(u3)
                ds4[i, j] = norm(u4)
            end
            end
            =#
            # Then, compute the angle α
            function computeα(
                _λ :: Float64,
                _ϕ :: Float64,
                grid_north :: Array{Float64, 1},
            )

                true_north  = [ - sin(_ϕ) * cos(_λ), - sin(_ϕ) * sin(_λ),   cos(_ϕ) ]
                true_east   = [ - sin(_λ)          ,   cos(_λ)          ,   0.0     ]
                # true_upward = [   cos(_ϕ) * cos(_λ),   cos(_ϕ) * sin(_λ),   sin(_ϕ) ]

                # α is the same "ANGLE" variable in POP2 output. 
                # It is the angle between true east and grid east (true north and grid north)
                # positive means the grid is rotated counterclockwise looking top-down
                cos_α =   (grid_north ⋅ true_north) / norm(grid_north) 
                sin_α = - (grid_north ⋅ true_east)  / norm(grid_north)

                # Prevent some numerical errors such that
                # cos_α is not within [-1, 1]
                if abs(cos_α) > 1.0
                    if abs(cos_α) - 1 < 1e-3   # some tolerance
                        cos_α = 1.0 * sign(cos_α)
                    else
                        throw(ErrorException("Angle calculation goes wrong. Please check.")) 
                    end
                end

                _α = acos(cos_α) # acos only gives number within [0, π]

                if sin_α < 0.0 # This means α is in [π, 2π]
                    _α = 2*π - _α
                end

                return _α, cos(_α), sin(_α)
            end


            #println("Compute grid curved angle on T grid")
            for i=1:Nx, j=1:Ny
                
                α_t[i, j], cosα_t[i, j], sinα_t[i, j] = computeα(
                    λ_t[i, j],
                    ϕ_t[i, j],
                    u2[:, i, j] + (- u4[:, i, j]),
                )

            end
            
            #println("Compute grid curved angle on UV grid")
            for i=1:Nx, j=1:Ny+1
                
                if j == 1
                    _λ = vs_λ[1, i, 1]
                    _ϕ = vs_ϕ[1, i, 1]
                    grid_north = - u4[:, i, 1]

                elseif j == Ny+1
                    _λ = vs_λ[4, i, Ny]
                    _ϕ = vs_ϕ[4, i, Ny]
                    grid_north = - u4[:, i, Ny]

                else
                    _λ = vs_λ[1, i, j]
                    _ϕ = vs_ϕ[1, i, j]
                    grid_north = - u4[:, i, j-1] - u4[:, i, j]
                end

                α_uv[i, j], cosα_uv[i, j], sinα_uv[i, j] = computeα(
                    _λ,
                    _ϕ,
                    grid_north,
                )
                
            end

            cvt13 = (arr, nx, ny) -> repeat(reshape(arr, :, 1, 1), outer=(1, nx, ny))
            cvt23 = (arr, nz) -> repeat( reshape(arr, 1, size(arr)...), outer=(nz, 1, 1) )

            λ_T  = cvt23(λ_t, Nz)
            λ_UV = zeros(Float64, Nx, Ny+1)
            λ_UV[:, 1:Ny] = vs_λ[1, :, :]
            λ_UV[:, end ] = vs_λ[4, :, end]
            λ_UV = cvt23(λ_UV, Nz) 

            ϕ_T = cvt23(ϕ_t, Nz)
            ϕ_UV = zeros(Float64, Nx, Ny+1)
            ϕ_UV[:, 1:Ny] = vs_ϕ[1, :, :]
            ϕ_UV[:, end] = vs_ϕ[4, :, end]

            ϕ_U = (ϕ_UV[:, 1:Ny] + ϕ_UV[:, 2:Ny+1]) / 2.0
            ϕ_V = (ϕ_UV + circshift(ϕ_UV, (-1, 0))) / 2.0

            ϕ_U  = cvt23(ϕ_U,  Nz) 
            ϕ_V  = cvt23(ϕ_V,  Nz) 
            ϕ_UV = cvt23(ϕ_UV, Nz) 

            z_t = ( z_w[1:end-1] + z_w[2:end] ) / 2.0
            z_T = cvt13(z_t, Nx, Ny)
            z_W = cvt13(z_w, Nx, Ny)
            z_UV = cvt13(z_t, Nx, Ny+1)
            z_U = z_T
            z_V = z_UV

            #println("Nz, Nx, Ny = ", (Nz, Nx, Ny))

            Δx_t = (ds1 + ds3) / 2.0
            Δx_T = cvt23(Δx_t, Nz)
            Δx_U = (Δx_T + circshift(Δx_T, (0, 1, 0)))/2.0
            Δx_v = zeros(Float64, Nx, Ny+1)
            Δx_v[:, 1:Ny] = ds1
            Δx_v[:, Ny+1] = ds4[:, end]
            Δx_V = cvt23(Δx_v, Nz)
            Δx_W = cvt23(Δx_t, Nz+1)
            Δx_uv = zeros(Float64, Nx, Ny+1)
            Δx_uv[:, 1:Ny] = (ds1 + circshift(ds1, (1, 0))) / 2.0
            Δx_uv[:, Ny+1] = (ds3[:, Ny] + circshift(ds3, (1, 0))[:, Ny]) / 2.0
            Δx_UV = cvt23(Δx_uv, Nz)

            Δy_t = (ds2 + ds4) / 2.0
            Δy_T = cvt23(Δy_t, Nz)
            Δy_U = cvt23(ds4, Nz)
            Δy_v = zeros(Float64, Nx, Ny+1)
            Δy_v[:, 2:Ny] = (Δy_t[:, 1:Ny-1] + Δy_t[:, 2:Ny]) / 2.0
            Δy_v[:, 1]    = Δy_v[:, 2]
            Δy_v[:, Ny+1] = Δy_v[:, Ny]
            Δy_V = cvt23(Δy_v, Nz)
            Δy_W = cvt23(Δy_t, Nz+1)
            Δy_uv = zeros(Float64, Nx, Ny+1)
            Δy_uv[:, 2:Ny] = (ds4[:, 1:Ny-1] + ds4[:, 2:Ny]) / 2.0
            Δy_uv[:, 1] = Δy_uv[:, 2]
            Δy_uv[:, Ny+1] = Δy_uv[:, Ny]
            Δy_UV = cvt23(Δy_uv, Nz)


            Δz_t = z_w[1:end-1] - z_w[2:end]
            Δz_T = cvt13(Δz_t, Nx, Ny)
            Δz_U = cvt13(Δz_t, Nx, Ny)
            Δz_V = cvt13(Δz_t, Nx, Ny+1)
            Δz_w = zeros(Float64, Nz+1)
            Δz_w[2:Nz] = (Δz_t[1:end-1] + Δz_t[2:end]) / 2.0
            Δz_w[1] = Δz_w[2]
            Δz_w[end] = Δz_w[end-1]
            Δz_W  = cvt13(Δz_w, Nx, Ny)
            Δz_UV = cvt13(Δz_t, Nx, Ny+1)

            f_T = 2Ω * sin.(ϕ_T)

            if sub_yrng == Colon()
                sub_yrng = 1:Ny
            end

            new_Ny = length(sub_yrng)
            sub_yrng_ext = (sub_yrng[1]):(sub_yrng[end]+1)
          
            return new(
                R,
                Ω,
                Nx,
                new_Ny,
                Nz,
                ds1[:, sub_yrng],
                ds2[:, sub_yrng],
                ds3[:, sub_yrng],
                ds4[:, sub_yrng],
                α_t[:, sub_yrng],
                cosα_t[:, sub_yrng],
                sinα_t[:, sub_yrng],
                α_uv[:, sub_yrng_ext],
                cosα_uv[:, sub_yrng_ext],
                sinα_uv[:, sub_yrng_ext],
                λ_T[:, :, sub_yrng],
                λ_UV[:, :, sub_yrng_ext],
                ϕ_T[:, :, sub_yrng],
                ϕ_U[:, :, sub_yrng],
                ϕ_V[:, :, sub_yrng_ext],
                ϕ_UV[:, :, sub_yrng_ext],
                z_T[:, :, sub_yrng],
                z_U[:, :, sub_yrng],
                z_V[:, :, sub_yrng_ext],
                z_W[:, :, sub_yrng],
                z_UV[:, :, sub_yrng_ext],


                Δx_T[:, :, sub_yrng],
                Δx_U[:, :, sub_yrng],
                Δx_V[:, :, sub_yrng_ext],
                Δx_W[:, :, sub_yrng],
                Δx_UV[:, :, sub_yrng_ext],

                Δy_T[:, :, sub_yrng],
                Δy_U[:, :, sub_yrng],
                Δy_V[:, :, sub_yrng_ext],
                Δy_W[:, :, sub_yrng],
                Δy_UV[:, :, sub_yrng_ext],

                Δz_T[:, :, sub_yrng],
                Δz_U[:, :, sub_yrng],
                Δz_V[:, :, sub_yrng_ext],
                Δz_W[:, :, sub_yrng],
                Δz_UV[:, :, sub_yrng_ext],
            )
     
        end
    end


    function project(
        gi    :: CurvilinearSphericalGrid,
        ivf_x :: AbstractArray{Float64, 3},     # input vector field east
        ivf_y :: AbstractArray{Float64, 3};     # input vector field north
        grid      :: Symbol,
        direction :: Symbol = :Forward,
    )

        ovf_x = similar(ivf_x)
        ovf_y = similar(ivf_y)

        project!(
            gi,
            ivf_x, ivf_y,
            ovf_x, ovf_y;
            direction=direction,
            grid=grid,
        )

        return ovf_x, ovf_y
    end


    function project!(
        gi    :: CurvilinearSphericalGrid,
        ivf_x :: AbstractArray{Float64, 3},     # input vector field east
        ivf_y :: AbstractArray{Float64, 3},     # input vector field north
        ovf_x :: AbstractArray{Float64, 3},    # output vector field east
        ovf_y :: AbstractArray{Float64, 3};    # output vector field north
        grid      :: Symbol,
        direction :: Symbol = :Forward,
    )

        if grid == :T
            COSα = gi.cosα_t
            SINα = gi.sinα_t
        elseif grid == :UV
            COSα = gi.cosα_uv
            SINα = gi.sinα_uv
        else
            throw(ErrorException("Unrecognized grid: " * string(grid)))
        end
        if direction == :Forward   # from outside world onto dispalced pole grid
            
            for k=1:size(ivf_x)[1]
                ovf_x_view = view(ovf_x, k, :, :)
                ovf_y_view = view(ovf_y, k, :, :)
                ivf_x_view = view(ivf_x, k, :, :)
                ivf_y_view = view(ivf_y, k, :, :)

                @. ovf_x_view =   ivf_x_view * COSα + ivf_y_view * SINα
                @. ovf_y_view = - ivf_x_view * SINα + ivf_y_view * COSα
            end
        elseif direction == :Backward   # from displaced pole grid onto outside world
            for k=1:size(ivf_x)[1]
                ovf_x_view = view(ovf_x, k, :, :)
                ovf_y_view = view(ovf_y, k, :, :)
                ivf_x_view = view(ivf_x, k, :, :)
                ivf_y_view = view(ivf_y, k, :, :)

                @. ovf_x_view = ivf_x_view * COSα - ivf_y_view * SINα
                @. ovf_y_view = ivf_x_view * SINα + ivf_y_view * COSα
            end
        else
            throw(ErrorException("Unsupported direction: " * string(direction)))
        end

    end

    function genGrid(
        gf       :: GridFile,
        z_w      :: AbstractArray{Float64, 1};
        sub_yrng :: Union{Colon, UnitRange} = Colon(),
    )

        local gi

        if typeof(gf) <: CurvilinearSphericalGridFile
            
            gi = CurvilinearSphericalGrid(;
                R=gf.R,
                Ω=gf.Ω,
                Nx=gf.Nx,
                Ny=gf.Ny,
                λ_t=gf.xc,
                ϕ_t=gf.yc,
                vs_λ=gf.xv,
                vs_ϕ=gf.yv,
                z_w=z_w,
                angle_unit=:deg,
                sub_yrng=sub_yrng,
            )
           
        end


        return gi
    end

end

=#

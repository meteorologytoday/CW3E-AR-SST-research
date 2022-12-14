if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end


module MITgcmTools

    using PyCall
    using Formatting
    using ..CoordinateModule
    
    function py2jl(x::AbstractArray)

        dim_len = length(size(x))
        permute = [dim_len - (i-1) for i=1:dim_len]

        return permutedims(x, permute)
    end



    function genDataDictionary(
        data      :: AbstractArray,
        metadata  :: Dict,
        copy_data :: Bool = false,
    )

        println(metadata)

        dimlist = metadata["dimlist"]
        Ndims = length(dimlist)
        
        spatial_selector = Tuple([ Colon() for i=1:Ndims ])
        
        data_dict = Dict()
        for (i, fldname) in enumerate(metadata["fldlist"])
            trimmed_fldname = fldname |> lstrip |> rstrip

            if copy_data
                _data = data[spatial_selector..., i]
            else
                _data = view(data, spatial_selector..., i)
            end

            data_dict[trimmed_fldname] = _data
        end
        
        return data_dict

    end
     
    #=

    Matthew only outputs grid box lengths for T, U, and V grids.
    The UV grid is omitted

    =#
    function readMITgcmGrid_MM(
        input_dir :: String = ".";
        DXG :: String = "DXG",
        DYG :: String = "DYG",
        DXC :: String = "DXC",
        DYC :: String = "DYC",
        DRF :: String = "DRF",
        DRC :: String = "DRC",
        RAC :: String = "RAC",
        XG  :: String = "XG",
        YG  :: String = "YG",
        XC  :: String = "XC",
        YC  :: String = "YC",
        RF  :: String = "RF",
        RC  :: String = "RC",
        maskInC :: String = "maskInC",
        maskInS :: String = "maskInS",
        maskInW :: String = "maskInW",
        verbose :: Bool = false,
        lev :: Union{Colon, UnitRange} = Colon(),
    )

        m = pyimport("MITgcmutils")


        c = Dict()

        for (varname, loaded_varname) in Dict(
            :DXG => DXG,
            :DYG => DYG,
            :DXC => DXC,
            :DYC => DYC,
            :DRC => DRC,
            :DRF => DRF,
            :RAC => RAC,
            :XG => XG,
            :YG => YG,
            :XC => XC,
            :YC => YC,
            :RF => RF,
            :RC => RC,
            :maskInC => maskInC,
            :maskInS => maskInS,
            :maskInW => maskInW,
        )
            c[varname] = m.mds.rdmds(joinpath(input_dir, loaded_varname))

            dim_len = length(size(c[varname]))

            permute = [dim_len - (i-1) for i=1:dim_len]

            verbose && println("$varname => ", size(c[varname]))
            c[varname] = permutedims(c[varname], permute)
        end

        if lev == Colon()
            lev = 1:length(c[:DRF])
        end

        lev_T = lev
        lev_W = lev[1]:(lev[end]+1)

        # Subset the z-coordinate     
        c[:DRF] = c[:DRF][lev_T]
        c[:DRC] = c[:DRC][lev_T]
        c[:RC] = c[:RC][lev_T]
        c[:RF] = c[:RF][lev_W]


        Nx, Ny = size(c[:DXG])
        Nz = length(c[:DRF])
        elt = eltype(c[:DXG])

        println(format("(Nx, Ny, Nz) = ({:d}, {:d}, {:d})", Nx, Ny, Nz))

        # ===== [ ?? ] ===== 

        ??_T  = zeros(elt, Nx, Ny, 1)
        ??_U  = zeros(elt, Nx+1, Ny, 1)
        ??_V  = zeros(elt, Nx, Ny+1, 1)
        ??_W  = zeros(elt, Nx, Ny, 1)
        ??_UV = zeros(elt, Nx+1, Ny+1, 1)

        ??_T[:, :, 1] = c[:XC]

        ??_U[1:Nx, :, 1] = c[:XG]
        ??_U[Nx+1, :, 1] = 2 * ??_U[Nx, :, 1] - ??_U[Nx-1, :, 1] # This is the best I can do

        ??_V[:, 2:Ny, 1] = (??_T[:, 1:Ny-1, :, 1] + ??_T[:, 2:Ny, 1] ) / 2
        ??_V[:,    1, 1] = ??_V[:,  2, 1]
        ??_V[:, Ny+1, 1] = ??_V[:, Ny, 1]

        ??_W .= ??_T
     

        # ===== [ ?? ] ===== 

        ??_T  = zeros(elt, Nx, Ny, 1)
        ??_U  = zeros(elt, Nx+1, Ny, 1)
        ??_V  = zeros(elt, Nx, Ny+1, 1)
        ??_W  = zeros(elt, Nx, Ny, 1)
        ??_UV = zeros(elt, Nx+1, Ny+1, 1)

        ??_T[:, :, 1] = c[:YC]

        ??_V[:, 1:Ny, 1] = c[:YG]
        ??_V[:, Ny+1, 1] = 2 * ??_V[:, Ny, 1] - ??_V[: , Ny-1, 1] # This is the best I can do

        ??_U[2:Nx, :, 1] = (??_T[1:Nx-1, :, 1] + ??_T[2:Nx, :, 1] ) / 2
        ??_U[1,    :, 1] = ??_U[2,  :, 1]
        ??_U[Nx+1, :, 1] = ??_U[Nx, :, 1]

        ??_W .= ??_T
     
        # ===== [ z ] ===== 

        z_T  = zeros(elt, 1, 1, Nz)
        z_U  = zeros(elt, 1, 1, Nz)
        z_V  = zeros(elt, 1, 1, Nz)
        z_UV  = zeros(elt, 1, 1, Nz)
        z_W  = zeros(elt, 1, 1, Nz+1)

        z_T[1, 1, :] = c[:RC]
        z_U  .= z_T
        z_V  .= z_T
        z_UV .= z_T
        z_W[1, 1, :] = c[:RF]



        # ===== [ mask ] ===== 

        mask_T = zeros(elt, Nx, Ny, 1)
        mask_T[:, :, 1] = c[:maskInC]

        # ===== [ ??x ] ===== 

        ??x_T  = zeros(elt, Nx, Ny, 1)
        ??x_U  = zeros(elt, Nx+1, Ny, 1)
        ??x_V  = zeros(elt, Nx, Ny+1, 1)
        ??x_W  = zeros(elt, Nx, Ny, 1)
        ??x_UV = zeros(elt, Nx+1, Ny+1, 1)

        ??x_V[:, 1:Ny, 1] = c[:DXG]
        ??x_V[:, Ny+1, 1] = c[:DXG][:, Ny] # Repeat is the best I can do

        # It should be DXF and DYF but it is not output
        ??x_T[:, :, 1] = (??x_V[:, 1:Ny, 1] + ??x_V[:, 2:Ny+1, 1]) / 2
        ??x_W .= ??x_T

        ??x_U[1:Nx, :, 1] = c[:DXC]
        ??x_U[Nx+1, :, 1] = c[:DXC][Nx, :] # Repeat is the best I can do
        
        ??x_UV[:, 2:Ny, 1] = (??x_U[:, 1:Ny-1, 1] + ??x_U[:, 2:Ny, 1]) / 2 
        ??x_UV[:,   1, 1] = ??x_UV[:,     2, 1]
        ??x_UV[:, end, 1] = ??x_UV[:, end-1, 1]



        # ===== [ ??y ] ===== 

        ??y_T  = zeros(elt, Nx, Ny, 1)
        ??y_U  = zeros(elt, Nx+1, Ny, 1)
        ??y_V  = zeros(elt, Nx, Ny+1, 1)
        ??y_W  = zeros(elt, Nx, Ny, 1)
        ??y_UV = zeros(elt, Nx+1, Ny+1, 1)

        ??y_U[1:Nx, :, 1] = c[:DYG]
        ??y_U[Nx+1, :, 1] = c[:DYG][Nx, :] # Repeat is the best I can do

        # It should be DXF and DYF but it is not output
        ??y_T[:, :, 1] = (??y_U[1:Nx, :, 1] + ??y_U[2:Nx+1, :, 1]) / 2
        ??y_W .= ??y_T

        ??y_V[:, 1:Ny, 1] = c[:DYC]
        ??y_V[:, Ny+1, 1] = c[:DYC][:, Ny] # Repeat is the best I can do
        
        ??y_UV[:, 2:Ny, 1] = (??y_U[:, 1:Ny-1, 1] + ??y_U[:, 2:Ny, 1]) / 2 
        ??y_UV[:,   1, 1] = ??y_UV[:,     2, 1]
        ??y_UV[:, end, 1] = ??y_UV[:, end-1, 1]


        # ===== [ ??z ] ===== 

        ??z_T  = zeros(elt, 1, 1, Nz)
        ??z_U  = zeros(elt, 1, 1, Nz)
        ??z_V  = zeros(elt, 1, 1, Nz)
        ??z_W  = zeros(elt, 1, 1, Nz+1)
        ??z_UV = zeros(elt, 1, 1, Nz)

        ??z_T[1, 1, :] = c[:DRF]
        ??z_U .= ??z_T
        ??z_V .= ??z_T
        ??z_UV .= ??z_T 
        
        ??z_W[1, 1, 2:Nz] = c[:DRC][1:Nz-1]
        ??z_W[1, 1, 1]    = ??z_W[1, 1, 2]
        ??z_W[1, 1, Nz+1] = ??z_W[1, 1, Nz]



        # ===== [ ??a ] =====
        ??a_T  = zeros(elt, Nx, Ny, 1)
        ??a_T[:, :, 1] = c[:RAC]

        gd = CoordinateModule.Grid(

            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
        
            R = 6.371e6,
            
            ??_T  = ??_T,
            ??_U  = ??_U,
            ??_V  = ??_V,
            ??_W  = ??_W,
            ??_UV = ??_UV,

            ??_T  = ??_T,
            ??_U  = ??_U,
            ??_V  = ??_V,
            ??_W  = ??_W,
            ??_UV = ??_UV,

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

            ??x_T  = ??x_T,
            ??x_U  = ??x_U,
            ??x_V  = ??x_V,
            ??x_W  = ??x_W,
            ??x_UV = ??x_UV,

            ??y_T  = ??y_T,
            ??y_U  = ??y_U,
            ??y_V  = ??y_V,
            ??y_W  = ??y_W,
            ??y_UV = ??y_UV,

            ??z_T  = ??z_T,
            ??z_U  = ??z_U,
            ??z_V  = ??z_V,
            ??z_W  = ??z_W,
            ??z_UV = ??z_UV,

            ??a_T  = ??a_T,
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

        ??_U = zeros(eltype(c[:XG]), length(c[:XG])+1, 1, 1)
        ??_U[1:end-1, 1, 1] = c[:XG]
        ??_U[end, 1, 1] = c[:XG][end] + (c[:XG][end] - c[:XG][end-1])
     
        ??_V = zeros(eltype(c[:YG]), 1, length(c[:YG])+1, 1)
        ??_V[1, 1:end-1, 1] = c[:YG]
        ??_V[1,     end, 1] = c[:YG][end] + (c[:YG][end] - c[:YG][end-1])


    end



            

    =#

    function postprocessRdmds(bundle::Tuple)

        data     = bundle[1]
        itrs     = bundle[2]
        metadata = bundle[3]

        data = data |> py2jl
        data_dict = genDataDictionary(data, metadata)
        

        return data_dict, itrs, metadata
    end


end

if !(:CoordinateModule in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "CoordinateModule.jl")))
end

using PyCall
using Formatting

function readMITgcmData(
)

    mitgcm = pyimport("MITgcmutils")
    data_2D = mitgcm.mds.rdmds("$data_dir/diag_2D", iters,    lev=collect(lev)) |> py2jl

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
    lev :: Union{Colon(), UnitRange} = Colon(),
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

    # ===== [ λ ] ===== 

    λ_T  = zeros(elt, Nx, Ny, 1)
    λ_U  = zeros(elt, Nx+1, Ny, 1)
    λ_V  = zeros(elt, Nx, Ny+1, 1)
    λ_W  = zeros(elt, Nx, Ny, 1)
    λ_UV = zeros(elt, Nx+1, Ny+1, 1)

    λ_T[:, :, 1] = c[:XC]

    λ_U[1:Nx, :, 1] = c[:XG]
    λ_U[Nx+1, :, 1] = 2 * λ_U[Nx, :, 1] - λ_U[Nx-1, :, 1] # This is the best I can do

    λ_V[:, 2:Ny, 1] = (λ_T[:, 1:Ny-1, :, 1] + λ_T[:, 2:Ny, 1] ) / 2
    λ_V[:,    1, 1] = λ_V[:,  2, 1]
    λ_V[:, Ny+1, 1] = λ_V[:, Ny, 1]

    λ_W .= λ_T
 

    # ===== [ ϕ ] ===== 

    ϕ_T  = zeros(elt, Nx, Ny, 1)
    ϕ_U  = zeros(elt, Nx+1, Ny, 1)
    ϕ_V  = zeros(elt, Nx, Ny+1, 1)
    ϕ_W  = zeros(elt, Nx, Ny, 1)
    ϕ_UV = zeros(elt, Nx+1, Ny+1, 1)

    ϕ_T[:, :, 1] = c[:YC]

    ϕ_V[:, 1:Ny, 1] = c[:YG]
    ϕ_V[:, Ny+1, 1] = 2 * ϕ_V[:, Ny, 1] - ϕ_V[: , Ny-1, 1] # This is the best I can do

    ϕ_U[2:Nx, :, 1] = (ϕ_T[1:Nx-1, :, 1] + λ_T[2:Nx, :, 1] ) / 2
    ϕ_U[1,    :, 1] = ϕ_U[2,  :, 1]
    ϕ_U[Nx+1, :, 1] = ϕ_U[Nx, :, 1]

    ϕ_W .= ϕ_T
 
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

    # ===== [ Δx ] ===== 

    Δx_T  = zeros(elt, Nx, Ny, 1)
    Δx_U  = zeros(elt, Nx+1, Ny, 1)
    Δx_V  = zeros(elt, Nx, Ny+1, 1)
    Δx_W  = zeros(elt, Nx, Ny, 1)
    Δx_UV = zeros(elt, Nx+1, Ny+1, 1)

    Δx_V[:, 1:Ny, 1] = c[:DXG]
    Δx_V[:, Ny+1, 1] = c[:DXG][:, Ny] # Repeat is the best I can do

    # It should be DXF and DYF but it is not output
    Δx_T[:, :, 1] = (Δx_V[:, 1:Ny, 1] + Δx_V[:, 2:Ny+1, 1]) / 2
    Δx_W .= Δx_T

    Δx_U[1:Nx, :, 1] = c[:DXC]
    Δx_U[Nx+1, :, 1] = c[:DXC][Nx, :] # Repeat is the best I can do
    
    Δx_UV[:, 2:Ny, 1] = (Δx_U[:, 1:Ny-1, 1] + Δx_U[:, 2:Ny, 1]) / 2 
    Δx_UV[:,   1, 1] = Δx_UV[:,     2, 1]
    Δx_UV[:, end, 1] = Δx_UV[:, end-1, 1]



    # ===== [ Δy ] ===== 

    Δy_T  = zeros(elt, Nx, Ny, 1)
    Δy_U  = zeros(elt, Nx+1, Ny, 1)
    Δy_V  = zeros(elt, Nx, Ny+1, 1)
    Δy_W  = zeros(elt, Nx, Ny, 1)
    Δy_UV = zeros(elt, Nx+1, Ny+1, 1)

    Δy_U[1:Nx, :, 1] = c[:DYG]
    Δy_U[Nx+1, :, 1] = c[:DYG][Nx, :] # Repeat is the best I can do

    # It should be DXF and DYF but it is not output
    Δy_T[:, :, 1] = (Δy_U[1:Nx, :, 1] + Δy_U[2:Nx+1, :, 1]) / 2
    Δy_W .= Δy_T

    Δy_V[:, 1:Ny, 1] = c[:DYC]
    Δy_V[:, Ny+1, 1] = c[:DYC][:, Ny] # Repeat is the best I can do
    
    Δy_UV[:, 2:Ny, 1] = (Δy_U[:, 1:Ny-1, 1] + Δy_U[:, 2:Ny, 1]) / 2 
    Δy_UV[:,   1, 1] = Δy_UV[:,     2, 1]
    Δy_UV[:, end, 1] = Δy_UV[:, end-1, 1]


    # ===== [ Δz ] ===== 

    Δz_T  = zeros(elt, 1, 1, Nz)
    Δz_U  = zeros(elt, 1, 1, Nz)
    Δz_V  = zeros(elt, 1, 1, Nz)
    Δz_W  = zeros(elt, 1, 1, Nz+1)
    Δz_UV = zeros(elt, 1, 1, Nz)

    Δz_T[1, 1, :] = c[:DRF]
    Δz_U .= Δz_T
    Δz_V .= Δz_T
    Δz_UV .= Δz_T 
    
    Δz_W[1, 1, 2:Nz] = c[:DRC][1:Nz-1]
    Δz_W[1, 1, 1]    = Δz_W[1, 1, 2]
    Δz_W[1, 1, Nz+1] = Δz_W[1, 1, Nz]



    # ===== [ Δa ] =====
    Δa_T  = zeros(elt, Nx, Ny, 1)
    Δa_T[:, :, 1] = c[:RAC]

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

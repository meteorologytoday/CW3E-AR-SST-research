using PyCall

function loadMITGCM(
    input_dir :: String = ".";
    DXG :: String = "DXG",
    DYG :: String = "DYG",
    DXC :: String = "DXC",
    DYC :: String = "DYC",
)

    m = pyimport("MITgcmutils")


    coor = Dict()

    for varname, loaded_varname in Dict(
        :XG => XG,
        :YG => YG,
        :RF => RF,
    )
        coor[varname] = m.mds.rdmds(joinpath(input_dir, varname))[:]
    end

    z_W = reshape(coor[:RF], 1, 1, :)

    λ_U = zeros(eltype(coor[:XG]), length(coor[:XG])+1, 1, 1)
    λ_U[1:end-1, 1, 1] = coor[:XG]
    λ_U[end, 1, 1] = coor[:XG][end] + (coor[:XG][end] - coor[:XG][end-1])
 
    ϕ_V = zeros(eltype(coor[:YG]), 1, length(coor[:YG])+1, 1)
    ϕ_V[1, 1:end-1, 1] = coor[:YG]
    ϕ_V[1,     end, 1] = coor[:YG][end] + (coor[:YG][end] - coor[:YG][end-1])


end



        




function loadMITGCM(
    input_dir :: String = ".";
    XG :: String = "XG",
    YG :: String = "YG",
    RF :: String = "RF",
)

    m = pyimport("MITgcmutils")


    coor = Dict()

    for varname, loaded_varname in Dict(
        :XG => XG,
        :YG => YG,
        :RF => RF,
    )
        coor[varname] = m.mds.rdmds(joinpath(input_dir, varname))[:]
    end

    z_W = reshape(coor[:RF], 1, 1, :)

    λ_U = zeros(eltype(coor[:XG]), length(coor[:XG])+1, 1, 1)
    λ_U[1:end-1, 1, 1] = coor[:XG]
    λ_U[end, 1, 1] = coor[:XG][end] + (coor[:XG][end] - coor[:XG][end-1])
 
    ϕ_V = zeros(eltype(coor[:YG]), 1, length(coor[:YG])+1, 1)
    ϕ_V[1, 1:end-1, 1] = coor[:YG]
    ϕ_V[1,     end, 1] = coor[:YG][end] + (coor[:YG][end] - coor[:YG][end-1])


end



        



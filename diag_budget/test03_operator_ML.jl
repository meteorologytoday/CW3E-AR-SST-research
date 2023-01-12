
include("Operators_ML.jl")
include("TestTools.jl")

using PyCall
using JLD2



println("*** Testing Operators_ML ***")


raw_λ_U = collect(range(0.0, 50, length=11))
raw_ϕ_V = collect(range(-30, 30, length=13))
raw_z_W = collect(range(0.0, -500.0, length=51))

coo = TestTools.genGrid(raw_λ_U, raw_ϕ_V, raw_z_W)

Nx, Ny, Nz = coo.gd.Nx, coo.gd.Ny, coo.gd.Nz

TEMP = zeros(Float64, Nx, Ny, Nz)
SALT = zeros(Float64, Nx, Ny, Nz)
ONES =  ones(Float64, Nx, Ny, Nz)

TEMP[:, :, 1:5]   .= 30.0
TEMP[:, :, 6:end] .= 20.0

SALT[:, :, 1:5]   .= 35.0
SALT[:, :, 6:end] .= 36.0

#ONES[:, :, 5] .= 2


println("Detect MLD")
h, bundle = Operators_ML.detectMLD(TEMP, SALT, coo)


println("Integrate ONES")
ML∫ONESdz = Operators_ML.sT_ML∫dz_T(ONES, coo; h=h, do_avg=true)

println("===== h =====")
println(h)

println("===== ML∫ONESdz =====")
println(ML∫ONESdz)

println("ρ")
println(bundle["ρ"][1, 1, :])   
 




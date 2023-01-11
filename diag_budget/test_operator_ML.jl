using PyCall
using JLD2

fill_value = 1e20

include("Operators_ML.jl")

println("*** Testing Operators_ML ***")

coo = MITgcmTools.readMITgcmGrid_MM(grid_dir, verbose=true, lev=lev)

#iters = 386400
iters = 388800
data_3D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_state", iters, lev=collect(lev), returnmeta=true))
data_2D, _, _ = MITgcmTools.postprocessRdmds(mitgcm.mds.rdmds("$data_dir/diag_2D", iters, returnmeta=true))

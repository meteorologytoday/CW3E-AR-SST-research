using PyCall
using JLD2

include("read_mitgcm_grid.jl")

println("*** Testing mitgcm ***")


operator_file = "grid_and_mtxops.jld2"


println("Loading operators...")

if isfile(operator_file)
    
    println("File $operator_file exist! Load it now...")
    
    d = JLD2.load(operator_file)

    gsp = d["gsp"]
    amo = d["amo"]

else
    
    println("File $operator_file does not exist. Generate it now...")
    gsp = readMITgcmGrid_MM("/data/SO2/SWOT/GRID/BIN")
    amo = MatrixOperators.AdvancedMatrixOperators(gsp=gsp)

    println("Saving operators...")
    JLD2.save(operator_file, "gsp", gsp, "amo", amo)
    println("Done.")
end

println("Loading operators complete.")



mitgcm = pyimport("MITgcmutils")

# Read data
# Read grid
# Compute h
# Compute T_ml T'
# Compute v_ml v'
# Compute gradients
# Compute vertical means
# Compute convergences
# Compute residue
# Output diagnostics

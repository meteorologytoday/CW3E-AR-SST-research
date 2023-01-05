

module CoordinateModule


    include("Grid.jl")
    include("GridSpacing.jl")


    mutable struct Coordinate
        gd :: Grid
        gsp :: GridSpacing
    end

end

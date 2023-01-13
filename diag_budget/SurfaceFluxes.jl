module SurfaceFluxes


    ζ_1 = 0.6  # meters
    ζ_2 = 20.0 # meters
    frac_1 = 0.62
    frac_2 = 1.0 - frac_lw
    cutoff_depth = 200.0 

    SFCFLX_shape(z) = z .>= 0
    RADFLX_shape(z) = ( frac_1 * exp.(z/ζ_1) + frac_2 * exp.(z/ζ_2) ) * (z .>=  -cutoff_depth)




end

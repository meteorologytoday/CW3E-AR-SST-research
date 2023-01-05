Base.@kwdef struct Grid

        Nx    :: Integer
        Ny    :: Integer
        Nz    :: Integer

        R     :: Float64
        
        λ_T  :: AbstractArray{Float64, 3}
        λ_U  :: AbstractArray{Float64, 3}
        λ_V  :: AbstractArray{Float64, 3}
        λ_W  :: AbstractArray{Float64, 3}
        λ_UV :: AbstractArray{Float64, 3}
         
        ϕ_T :: AbstractArray{Float64, 3}
        ϕ_U :: AbstractArray{Float64, 3}
        ϕ_V :: AbstractArray{Float64, 3}
        ϕ_W :: AbstractArray{Float64, 3}
        ϕ_UV :: AbstractArray{Float64, 3}
        
        z_T :: AbstractArray{Float64, 3}
        z_U :: AbstractArray{Float64, 3}
        z_V :: AbstractArray{Float64, 3}
        z_W :: AbstractArray{Float64, 3}
        z_UV :: AbstractArray{Float64, 3}
        

        mask_T :: AbstractArray{Float64, 3}
        

end



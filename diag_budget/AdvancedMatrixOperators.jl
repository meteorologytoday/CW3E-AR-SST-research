mutable struct AdvancedMatrixOperators

    bmo       :: BasicMatrixOperators
    gsp        :: GridSpacing

    T_DIVx_U    :: AbstractArray{Float64, 2}
    T_DIVy_V    :: AbstractArray{Float64, 2}
    T_DIVz_W    :: AbstractArray{Float64, 2}
    
    U_∂x_T      :: AbstractArray{Float64, 2}
    V_∂y_T      :: AbstractArray{Float64, 2}
    W_∂z_T      :: AbstractArray{Float64, 2}

    T_∂x_U      :: AbstractArray{Float64, 2}
    T_∂y_V      :: AbstractArray{Float64, 2}
    T_∂z_W      :: AbstractArray{Float64, 2}
    
    U_∂y_UV     :: AbstractArray{Float64, 2}
    V_∂x_UV     :: AbstractArray{Float64, 2}

    T_CURLx_U   :: AbstractArray{Float64, 2}
    T_CURLy_V   :: AbstractArray{Float64, 2}

    T_CURLx_T   :: AbstractArray{Float64, 2}
    T_CURLy_T   :: AbstractArray{Float64, 2}

    U_interp_T :: AbstractArray{Float64, 2}  # interpolation of T grid onto U grid
    V_interp_T :: AbstractArray{Float64, 2}  # interpolation of T grid onto V grid
    W_interp_T :: AbstractArray{Float64, 2}  # interpolation of T grid onto W grid

    T_interp_UV :: AbstractArray{Float64, 2}  
    U_interp_UV :: AbstractArray{Float64, 2}  
    V_interp_UV :: AbstractArray{Float64, 2}  

    UV_interp_U :: AbstractArray{Float64, 2}  
    UV_interp_V :: AbstractArray{Float64, 2}  
 
    T_mask_T       :: AbstractArray{Float64, 2}
    U_mask_U       :: AbstractArray{Float64, 2}
    V_mask_V       :: AbstractArray{Float64, 2}
    W_mask_W       :: AbstractArray{Float64, 2}
    T_bordermask_T :: AbstractArray{Float64, 2}

    T_Δx_T :: AbstractArray{Float64, 2}
    T_Δy_T :: AbstractArray{Float64, 2}
    T_Δz_T :: AbstractArray{Float64, 2}
 
    U_Δx_U :: AbstractArray{Float64, 2}
    U_Δy_U :: AbstractArray{Float64, 2}
    U_Δz_U :: AbstractArray{Float64, 2}
 
    V_Δx_V :: AbstractArray{Float64, 2}
    V_Δy_V :: AbstractArray{Float64, 2}
    V_Δz_V :: AbstractArray{Float64, 2}
 
    W_Δx_W :: AbstractArray{Float64, 2}
    W_Δy_W :: AbstractArray{Float64, 2}
    W_Δz_W :: AbstractArray{Float64, 2}
 
    UV_Δx_UV :: AbstractArray{Float64, 2}
    UV_Δy_UV :: AbstractArray{Float64, 2}
    UV_Δz_UV :: AbstractArray{Float64, 2}
 
    T_invΔx_T :: AbstractArray{Float64, 2}
    T_invΔy_T :: AbstractArray{Float64, 2}
    T_invΔz_T :: AbstractArray{Float64, 2}
 
    U_invΔx_U :: AbstractArray{Float64, 2}
    U_invΔy_U :: AbstractArray{Float64, 2}
    U_invΔz_U :: AbstractArray{Float64, 2}
 
    V_invΔx_V :: AbstractArray{Float64, 2}
    V_invΔy_V :: AbstractArray{Float64, 2}
    V_invΔz_V :: AbstractArray{Float64, 2}
 
    W_invΔx_W :: AbstractArray{Float64, 2}
    W_invΔy_W :: AbstractArray{Float64, 2}
    W_invΔz_W :: AbstractArray{Float64, 2}
 
    UV_invΔx_UV :: AbstractArray{Float64, 2}
    UV_invΔy_UV :: AbstractArray{Float64, 2}
    UV_invΔz_UV :: AbstractArray{Float64, 2}
    
    T_Δvol_T    :: AbstractArray{Float64, 2}
    T_invΔvol_T :: AbstractArray{Float64, 2}
    
    function AdvancedMatrixOperators(;
        gsp           :: GridSpacing,
        mask_T         :: Union{Nothing, AbstractArray{Float64, 3}} = nothing,
        bmo :: Union{Nothing, BasicMatrixOperators} = nothing,
    )

        cvt_diagm = (x, d) -> spdiagm( 0 => view(repeat(x, outer=d), :) )


        Nx = gsp.Nx
        Ny = gsp.Ny
        Nz = gsp.Nz

        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt3_diagm_T = (x,) -> spdiagm( 0 => view(repeat(x, outer=(1, 1, Nz)), :) )
        cvt3_diagm_U = cvt3_diagm_T
        cvt3_diagm_V = cvt3_diagm_T
        cvt3_diagm_UV = cvt3_diagm_T
        cvt3_diagm_W = (x,) -> spdiagm( 0 => view(repeat(x, outer=(1, 1, Nz+1)), :) )


        cvt3_diagm_T_z  = (x,) -> spdiagm( 0 => view(repeat(x, outer=(Nx, Ny, 1)), :) )
        cvt3_diagm_U_z  = (x,) -> spdiagm( 0 => view(repeat(x, outer=(Nx+1, Ny, 1)), :) )
        cvt3_diagm_V_z  = (x,) -> spdiagm( 0 => view(repeat(x, outer=(Nx, Ny+1, 1)), :) )
        cvt3_diagm_UV_z = (x,) -> spdiagm( 0 => view(repeat(x, outer=(Nx+1, Ny+1, 1)), :) )
        cvt3_diagm_W_z  = cvt3_diagm_T_z





        if bmo == nothing 
            println("Construct BMO")
            @time bmo = BasicMatrixOperators(Nx=Nx, Ny=Ny, Nz=Nz)
        end

        if mask_T == nothing
            mask_T = ones(Float64, bmo.T_dim...)
        end
        if length(mask_T) != bmo.T_pts
            throw(ErrorException("Length of topo.mask_T does not conform"))
        end

       
        mask3_flat = view(mask_T,  :)

        onV_if_unblocked_north_onT = bmo.V_S_T  * mask3_flat
        onV_if_unblocked_south_onT = bmo.V_N_T  * mask3_flat
        onU_if_unblocked_east_onT  = bmo.U_W_T  * mask3_flat
        onU_if_unblocked_west_onT  = bmo.U_E_T  * mask3_flat
        onW_if_unblocked_up_onT    = bmo.W_DN_T * mask3_flat
        onW_if_unblocked_dn_onT    = bmo.W_UP_T * mask3_flat

        V_mask = onV_if_unblocked_north_onT .* onV_if_unblocked_south_onT
        U_mask = onU_if_unblocked_east_onT  .* onU_if_unblocked_west_onT
        W_mask = onW_if_unblocked_up_onT    .* onW_if_unblocked_dn_onT

        T_mask_T = spdiagm(0 => mask3_flat)
        V_mask_V = spdiagm(0 => V_mask)
        U_mask_U = spdiagm(0 => U_mask)
        W_mask_W = spdiagm(0 => W_mask)

        # border in horizontal direction
        T_bordermask_T = spdiagm( 0 => ( 
               (bmo.T_N_T  * mask3_flat)
            .* (bmo.T_S_T  * mask3_flat)
            .* (bmo.T_E_T  * mask3_flat)
            .* (bmo.T_W_T  * mask3_flat)
        ))

        # ===== [ BEGIN face area and lengths on U V ] =====

        T_Δx_T = (gsp.Δx_T |>  cvt3_diagm_T)
        T_Δy_T = (gsp.Δy_T |>  cvt3_diagm_T)
        T_Δz_T = (gsp.Δz_T |>  cvt3_diagm_T_z)
 
        U_Δx_U = (gsp.Δx_U |>  cvt3_diagm_U)
        U_Δy_U = (gsp.Δy_U |>  cvt3_diagm_U)
        U_Δz_U = (gsp.Δz_U |>  cvt3_diagm_U_z)
 
        V_Δx_V = (gsp.Δx_V |>  cvt3_diagm_V)
        V_Δy_V = (gsp.Δy_V |>  cvt3_diagm_V)
        V_Δz_V = (gsp.Δz_V |>  cvt3_diagm_V_z)
 
        W_Δx_W = (gsp.Δx_W |>  cvt3_diagm_W)
        W_Δy_W = (gsp.Δy_W |>  cvt3_diagm_W)
        W_Δz_W = (gsp.Δz_W |>  cvt3_diagm_W_z)
 
        UV_Δx_UV = (gsp.Δx_UV |>  cvt3_diagm_UV)
        UV_Δy_UV = (gsp.Δy_UV |>  cvt3_diagm_UV)
        UV_Δz_UV = (gsp.Δz_UV |>  cvt3_diagm_UV_z)
 
        T_invΔx_T = (gsp.Δx_T.^(-1) |>  cvt3_diagm_T)
        T_invΔy_T = (gsp.Δy_T.^(-1) |>  cvt3_diagm_T)
        T_invΔz_T = (gsp.Δz_T.^(-1) |>  cvt3_diagm_T_z)
 
        U_invΔx_U = (gsp.Δx_U.^(-1) |>  cvt3_diagm_U)
        U_invΔy_U = (gsp.Δy_U.^(-1) |>  cvt3_diagm_U)
        U_invΔz_U = (gsp.Δz_U.^(-1) |>  cvt3_diagm_U_z)
 
        V_invΔx_V = (gsp.Δx_V.^(-1) |>  cvt3_diagm_V)
        V_invΔy_V = (gsp.Δy_V.^(-1) |>  cvt3_diagm_V)
        V_invΔz_V = (gsp.Δz_V.^(-1) |>  cvt3_diagm_V_z)
 
        W_invΔx_W = (gsp.Δx_W.^(-1) |>  cvt3_diagm_W)
        W_invΔy_W = (gsp.Δy_W.^(-1) |>  cvt3_diagm_W)
        W_invΔz_W = (gsp.Δz_W.^(-1) |>  cvt3_diagm_W_z)
 
        UV_invΔx_UV = (gsp.Δx_UV.^(-1) |>  cvt3_diagm_UV)
        UV_invΔy_UV = (gsp.Δy_UV.^(-1) |>  cvt3_diagm_UV)
        UV_invΔz_UV = (gsp.Δz_UV.^(-1) |>  cvt3_diagm_UV_z)

        T_Δa_T    = gsp.Δa_T         |> cvt3_diagm_T
        T_invΔa_T = (gsp.Δa_T.^(-1)) |> cvt3_diagm_T

        T_Δvol_T = T_Δa_T * T_Δz_T
        T_invΔvol_T = T_invΔa_T * T_invΔz_T

        # Δz is special. Need to clear NaN
        function clearNaN!(m)
            for i = 1:length(m.nzval)
                if isnan(m.nzval[i])
                    m.nzval[i] = 0
                end
            end
            dropzeros!(m)
        end

        # ===== [ END face area and lengths on T ] =====
        
        #println("Making derivatives") 

        # ===== [ BEG making matrix ] =====
        # MAGIC!!

        T_DIVx_U = T_mask_T * T_invΔvol_T * (bmo.T_W_U  - bmo.T_E_U  ) * U_Δz_U * U_Δy_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V = T_mask_T * T_invΔvol_T * (bmo.T_S_V  - bmo.T_N_V  ) * V_Δz_V * V_Δx_V  ; dropzeros!(T_DIVy_V);
        T_DIVz_W = T_mask_T * T_invΔvol_T * (bmo.T_DN_W - bmo.T_UP_W ) * W_Δx_W * W_Δy_W  ; dropzeros!(T_DIVz_W);

        U_∂x_T = U_mask_U * U_invΔx_U * (bmo.U_W_T  - bmo.U_E_T)                 ; dropzeros!(U_∂x_T);
        V_∂y_T = V_mask_V * V_invΔy_V * (bmo.V_S_T  - bmo.V_N_T)                 ; dropzeros!(V_∂y_T);
        W_∂z_T = W_mask_W * W_invΔz_W * (bmo.W_DN_T - bmo.W_UP_T)                ; dropzeros!(W_∂z_T);

        T_∂x_U  = T_mask_T * T_invΔx_T * ( bmo.T_W_U - bmo.T_E_U )               ; dropzeros!(T_∂x_U);
        T_∂y_V  = T_mask_T * T_invΔy_T * ( bmo.T_S_V - bmo.T_N_V )               ; dropzeros!(T_∂y_V);
        T_∂z_W  = T_mask_T * T_invΔz_T * ( bmo.T_DN_W - bmo.T_UP_W )             ; dropzeros!(T_∂z_W);

        U_∂y_UV = U_mask_U * U_invΔy_U * (bmo.U_S_UV  - bmo.U_N_UV)              ; dropzeros!(U_∂y_UV);
        V_∂x_UV = V_mask_V * V_invΔx_V * (bmo.V_W_UV  - bmo.V_E_UV)              ; dropzeros!(V_∂x_UV);



        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse
            for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] ./= wgt
                end
            end
          
            return dropzeros(transpose(m_t) |> sparse)
        end

        #println("Making interpolations part 1") 
        ones_T  = ones(Float64, bmo.T_pts)
        ones_U  = ones(Float64, bmo.U_pts)
        ones_V  = ones(Float64, bmo.V_pts)
        ones_UV = ones(Float64, bmo.UV_pts)
        
        #println("Making interpolations part 2") 

        U_interp_T = (bmo.U_W_T + bmo.U_E_T) * T_mask_T
        U_interp_T = selfDivision(U_interp_T, ones_T)

        V_interp_T = (bmo.V_S_T + bmo.V_N_T) * T_mask_T
        V_interp_T = selfDivision(V_interp_T, ones_T)

        W_interp_T = (bmo.W_DN_T + bmo.W_UP_T) * T_mask_T
        W_interp_T = selfDivision(W_interp_T, ones_T)

        # Notice the configuration might make T grid value
        # be produced on a masked grid if it is surrounded
        # by four valid UV pts. So, we need to put a T mask
        # at the end
        T_interp_UV = T_mask_T * (bmo.T_E_U + bmo.T_W_U) * (bmo.U_S_UV + bmo.U_N_UV)
        T_interp_UV = selfDivision(T_interp_UV, ones_UV)
        
        U_interp_UV = U_mask_U * (bmo.U_S_UV + bmo.U_N_UV)
        U_interp_UV = selfDivision(U_interp_UV, ones_UV)

        V_interp_UV = V_mask_V * (bmo.V_W_UV + bmo.V_E_UV)
        V_interp_UV = selfDivision(V_interp_UV, ones_UV)

        UV_interp_U = (bmo.UV_N_U + bmo.UV_S_U)
        UV_interp_U = selfDivision(UV_interp_U, ones_U)

        UV_interp_V = (bmo.UV_E_V + bmo.UV_W_V)
        UV_interp_V = selfDivision(UV_interp_V, ones_V)


        #sumΔvol_T = reshape(ones(Float64, bmo.T_pts), 1, :) * T_mask_T * T_Δvol_T
        
        # CURL operators needs interpolation in practical usage

        # used for ∂v/∂x
        T_CURLx_U =   T_bordermask_T * T_invΔa_T * (bmo.T_W_U - bmo.T_E_U) * U_Δy_U
        T_CURLy_V = - T_bordermask_T * T_invΔa_T * (bmo.T_S_V - bmo.T_N_V) * V_Δx_V

        T_CURLx_T = T_CURLx_U * U_interp_T
        T_CURLy_T = T_CURLy_V * V_interp_T


        return new(

            bmo,
            gsp,

            T_DIVx_U,
            T_DIVy_V,
            T_DIVz_W,
            
            U_∂x_T,
            V_∂y_T,
            W_∂z_T,

            T_∂x_U,
            T_∂y_V,
            T_∂z_W,
            
            U_∂y_UV,
            V_∂x_UV,

            T_CURLx_U,
            T_CURLy_V,

            T_CURLx_T,
            T_CURLy_T,

            U_interp_T,
            V_interp_T,
            W_interp_T,

            T_interp_UV,
            U_interp_UV,
            V_interp_UV,

            UV_interp_U,
            UV_interp_V,
         
            T_mask_T,
            U_mask_U,
            V_mask_V,
            W_mask_W,
            T_bordermask_T,

            T_Δx_T,
            T_Δy_T,
            T_Δz_T,
         
            U_Δx_U,
            U_Δy_U,
            U_Δz_U,
         
            V_Δx_V,
            V_Δy_V,
            V_Δz_V,
         
            W_Δx_W,
            W_Δy_W,
            W_Δz_W,
         
            UV_Δx_UV,
            UV_Δy_UV,
            UV_Δz_UV,
         
            T_invΔx_T,
            T_invΔy_T,
            T_invΔz_T,
         
            U_invΔx_U,
            U_invΔy_U,
            U_invΔz_U,
         
            V_invΔx_V,
            V_invΔy_V,
            V_invΔz_V,
         
            W_invΔx_W,
            W_invΔy_W,
            W_invΔz_W,
         
            UV_invΔx_UV,
            UV_invΔy_UV,
            UV_invΔz_UV,
            
            T_Δvol_T,
            T_invΔvol_T,
            
        )
    end
end

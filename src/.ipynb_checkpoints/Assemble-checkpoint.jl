###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################


    
function Massemble!(NGLL, NelX, NelY, dxe, dye, ThickX, 
                   ThickY, rho1, vs1, rho2, vs2, iglob, 
                   M, x, y, jac)


    xgll, wgll, H = GetGLL(NGLL)
    wgll2 = wgll*wgll';

    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    
    vso = zeros(NGLL, NGLL)
    vs = zeros(NGLL-1, NGLL)
    dx = zeros(NGLL-1, NGLL)
    muMax = 0
    dt = Inf
    
    # damage zone index
    damage_idx = zeros(Int, NelX*NelY)

    @inbounds @fastmath for ey = 1:NelY
        @inbounds @fastmath for ex = 1:NelX

            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]

            # Properties of heterogeneous medium
            if ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)
                damage_idx[eo] = eo
                rho[:,:] .= rho2
                mu[:,:] .= rho2*vs2^2
            else
                rho[:,:] .= rho1
                mu[:,:] .= rho1*vs1^2
            end

            if muMax < maximum(maximum(mu))
                muMax = maximum(maximum(mu))
            end

            # Diagonal Mass Matrix
            M[ig] .+= wgll2.*rho*jac

            # Local contributions to the stiffness matrix
            #  W[:,:,eo] .= wgll2.*mu;
            
            # Set timestep
            vso .= sqrt.(mu./rho)
            
            if dxe<dye
                vs .= max.(vso[1:NGLL-1,:], vso[2:NGLL,:])
                dx .= repeat( diff(xgll)*0.5*dxe, 1, NGLL)
            else
                vs .= max.(vso[:,1:NGLL-1], vso[:,2:NGLL])'
                dx .= repeat( diff(xgll)*0.5*dye, 1, NGLL)
            end
            
            dtloc = dx./vs
            dt = minimum( push!(dtloc[1:end], dt) )

        end
    end

    return M,dt, muMax, damage_idx[damage_idx .> 0]
end

# Material properties for a narrow rectangular damaged zone of half-thickness ThickY and depth ThickX 
function material_properties(NelX, NelY,NGLL, dxe, dye, ThickX, ThickY, wgll2, rho1, vs1, rho2, vs2)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    
    # damage zone index
    #  damage_idx = zeros(Int, NelX*NelY)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            
            # Properties of heterogeneous medium
            if ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)
                mu .= rho2*vs2^2
                #  damage_idx[eo] = eo
            else
                mu .= rho1*vs1^2
            end
            W[:,:,eo] = wgll2.*mu
        end
    end
    W
end



# Material properties for a trapezium damage zone
# Linear function: shape of trapezoid
function line(x,y)
    P1 = [0 3e3]
    P2 = [-8e3 1.5e3]

    f = (y - P2[2]) - ((P1[2]-P2[2])/(P1[1]-P2[1]))*(x - P2[1])

    return f
end

# Set up trapezoidal rigidity
function rigid(x,y)
    # Rigidity: host rock and fault zone
    rho1::Float64 = 2670
    vs1::Float64 = 3464
    
    rho2 = 0.6*rho1
    vs2 = 0.6*vs1
    rho3 = 0.8*rho1
    vs3 = 0.8*vs1
    
    rhoglob::Array{Float64} = zeros(length(x))
    vsglob::Array{Float64} = zeros(length(x))

    for i = 1:length(x)
        if x[i] > -8e3
            if line(x[i],y[i]) < 0
                rhoglob[i] = rho3
                vsglob[i] = vs3
            else
                rhoglob[i] = rho1
                vsglob[i] = vs1
            end
        else
            rhoglob[i] = rho1
            vsglob[i] = vs1
        end

    end

    for i = 1:length(x)
        if y[i]<0.25e3
            rhoglob[i] = rho2
            vsglob[i] = vs2
        end
    end

    
    return rhoglob, vsglob
end

function mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x,y, wgll2)
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    rhoglob, vsglob = rigid(x,y)
    muglob = rhoglob.*(vsglob.^2)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]
            
            mu[:,:] = muglob[ig]
            rho[:,:] = rhoglob[ig]
            
            W[:,:,eo] = wgll2.*mu
            M[ig] .+= wgll2.*rho*jac
        end
    end
    return M,W
end

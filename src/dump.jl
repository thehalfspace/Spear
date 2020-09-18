##########################################
## trying out stuff here. this file is not
## important for simulations
##########################################


## Testing elemental k matrix
# Single element stiffness
function stiff_element(NGLL, NelX, NelY, nglob, iglob, dxe, dye)
    xgll, wgll, H = GetGLL(NGLL)
    Ht = H'
    wgll2 = wgll*wgll'

    # Jacobians
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    c1::Float64 = jac/dx_dxi^2
    c2::Float64 = jac/dy_deta^2
    
    rho1::Float64 = 2500
    vs1::Float64 = 0.6*3464
    mu = 20
    Nel = 600;
    Ke2 = zeros(NGLL,NGLL,NGLL,NGLL)

    u = rand(5,5)

    W = wgll2.*mu
    del = Matrix{Float64}(I,NGLL,NGLL)  # identity matrix
    nn = 1
    n = 0; q=0; w=0 # iterator
    term1 = 0; term2 = 0
    for i in 1:5
        for j in 1:5
            term1 = 0; term2 = 0
            for k in 1:5
                for l in 1:5
                    term1 = 0; term2 = 0
                    for p in 1:5
                        term1 += del[i,k]*W[k,p]*(jac/dy_deta^2)*H[j,p]*H[l,p]
                        term2 += del[j,l]*W[p,j]*(jac/dx_dxi^2)*H[i,p]*H[k,p]
                    end
                    Ke2[i,j,k,l] = term1 + term2
                end
            end
        end
    end
   

    Ke = reshape(Ke2,NGLL*NGLL,NGLL*NGLL)


    # Calculate Ku for one element
    wloc = wgll2.*mu
    d_xi = Ht*u
    d_eta = u*H

    d_xi = H*(wloc.*d_xi)
    d_eta = (wloc.*d_eta)*Ht

    Ku = c1*d_xi + c2*d_eta
   
    return Ku, reshape(Ke*u[:],NGLL,NGLL)
        
end


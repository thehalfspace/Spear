###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################

# Gaussian function
function gauss(x, mu, sigma)
    return ((x .- mu)./(2*sigma)).^2
end

# Debug the setup
function rigid(x,y)

    # Rigidity: host rock and fault zone
    muhost = P.rho1*P.vs1^2
    mufz = P.rho2*P.vs2^2

    # Gaussian fault zone mean and std
    meanx = 0
    meany = 0
    sigx = (P.LX - P.ThickX)/3
    sigy = P.ThickY/3

    muglob = (mufz-muhost)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ muhost
    
    vsglob = (P.vs2-P.vs1)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ P.vs1
    rhoglob = (P.rho2-P.rho1)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ P.rho1
    
    return muglob, vsglob, rhoglob
end
    
function assemble(P::parameters, iglob, M, W, x, y)


    xgll, wgll, H = GetGLL(P.NGLL)
    wgll2 = wgll*wgll';

    rho::Matrix{Float64} = zeros(P.NGLL, P.NGLL)
    mu::Matrix{Float64} = zeros(P.NGLL, P.NGLL)
    
    vso = zeros(P.NGLL, P.NGLL)
    vs = zeros(P.NGLL-1, P.NGLL)
    dx = zeros(P.NGLL-1, P.NGLL)
    muMax = 0
    dt = Inf
    
    muglob, vsglob, rhoglob = rigid(x,y)

    #  # Rigidity: host rock and fault zone
    #  muhost = P.rho1*P.vs1^2
    #  mufz = P.rho2*P.vs2^2

    #  # Gaussian fault zone mean and std
    #  meanx = -P.LX
    #  meany = 0
    #  sigx = (P.LX - P.ThickX)/3
    #  sigy = P.ThickY/3

    for ey = 1:P.NelY
        for ex = 1:P.NelX

            eo = (ey-1)*P.NelX + ex
            ig = iglob[:,:,eo]

            mu[:,:] = muglob[ig]
            rho[:,:] = rhoglob[ig]

            if muMax < maximum(maximum(mu))
                muMax = maximum(maximum(mu))
            end

            # Diagonal Mass Matrix
            M[ig] .+= wgll2.*rho*P.jac

            # Local contributions to the stiffness matrix
            W[:,:,eo] .= wgll2.*mu;
            
            # Set timestep
            vso .= sqrt.(mu./rho)
            
            if P.dxe<P.dye
                vs .= max.(vso[1:P.NGLL-1,:], vso[2:P.NGLL,:])
                dx .= repeat( diff(xgll)*0.5*P.dxe, 1, P.NGLL)
            else
                vs .= max.(vso[:,1:P.NGLL-1], vso[:,2:P.NGLL])'
                dx .= repeat( diff(xgll)*0.5*P.dye, 1, P.NGLL)
            end
            
            dtloc = dx./vs
            dt = minimum( push!(dtloc[1:end], dt) )

        end
    end

    return M, W, dt, muMax
end


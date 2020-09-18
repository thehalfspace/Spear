###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################

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
    rho2 = 0.6*P.rho1
    vs2 = 0.6*P.vs1
    rho3 = 0.8*P.rho1
    vs3 = 0.8*P.rho1
    
    rhoglob::Array{Float64} = zeros(length(x))
    vsglob::Array{Float64} = zeros(length(x))

    for i = 1:length(x)
        if x[i] > -8e3
            if line(x[i],y[i]) < 0
                rhoglob[i] = rho3
                vsglob[i] = vs3
            else
                rhoglob[i] = P.rho1
                vsglob[i] = P.vs1
            end
        else
            rhoglob[i] = P.rho1
            vsglob[i] = P.vs1
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
    
function assemble(P::parameters, iglob, M, W, x, y)

    xgll, wgll, H = GetGLL(P.NGLL)
    wgll2 = wgll*wgll';
    
    rhoglob, vsglob = rigid(x,y)
    muglob = rhoglob.*(vsglob.^2)

    rho::Matrix{Float64} = zeros(P.NGLL, P.NGLL)
    mu::Matrix{Float64} = zeros(P.NGLL, P.NGLL)
    
    vso = zeros(P.NGLL, P.NGLL)
    vs = zeros(P.NGLL-1, P.NGLL)
    dx = zeros(P.NGLL-1, P.NGLL)
    muMax = 0
    dt = Inf

    # Rigidity: host rock and fault zone

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

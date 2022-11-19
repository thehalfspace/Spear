# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end


# Compute rate-state friciton with depth
function fricDepth(FltX)
    
    FltNglob = length(FltX)
    
    # Friction with depth
    cca::Array{Float64} = repeat([0.010], FltNglob)
    ccb::Array{Float64} = repeat([0.015], FltNglob)

    a_b = cca - ccb
    fP1 = [-0.005, 0e3]
    fP2 = [-0.005, -15e3]
    fP3 = [0.010, -18e3]
    fP4 = [0.010, -40e3]

    fric_depth1 = findall(abs.(FltX) .<= abs(fP2[2]))
    fric_depth2 = findall(abs(fP2[2]) .< abs.(FltX) .<= abs(fP3[2]))
    fric_depth3 = findall(abs(fP3[2]) .< abs.(FltX) .<= abs(fP4[2]))
    fric_depth4 = findall(abs.(FltX) .> abs(fP4[2]))

    a_b[fric_depth1] .= Int1D(fP1, fP2, FltX[fric_depth1])
    a_b[fric_depth2] .= Int1D(fP2, fP3, FltX[fric_depth2])
    a_b[fric_depth3] .= Int1D(fP3, fP4, FltX[fric_depth3])
    a_b[fric_depth4] .= 0.010

    #  cca[fric_depth4] .= Int1D(fP4, fP5, FltX[fric_depth4]) .+ 0.0001
    cca .= ccb .+ a_b
    #  ccb .= cca .- a_b

    return cca, ccb
end



# Effective normal stress
function SeffDepth(FltX)

    FltNglob = length(FltX)

    Seff::Array{Float64} = repeat([50e6], FltNglob)
    # sP1 = [10e6 0]
    # sP2 = [50e6 -2e3]
    # Seff_depth = findall(abs.(FltX) .<= abs(sP2[2]))
    # Seff[Seff_depth] = Int1D(sP1, sP2, FltX[Seff_depth])

    return Seff
end


# Shear stress
function tauDepth(FltX, rho, vs)

    FltNglob = length(FltX)

    Seff = 50e6
    a_max = 0.025
    cca = 0.01
    ccb = 0.015
    V = 1.0e-9
    Vo = 1.0e-6
    fo = 0.6

    mu = rho*vs^2

    tauo = Seff*a_max*asinh(V/(2*Vo)*exp((fo + ccb*log(Vo/V))/a_max))

    return tauo.*ones(FltNglob)
end

# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end


# Compute rate-state friciton with depth
function fricDepth(FltX)
    
    FltNglob = length(FltX)
    
    # Friction with depth
    cca::Array{Float64} = repeat([0.015], FltNglob)
    ccb::Array{Float64} = repeat([0.019], FltNglob)

    a_b = cca - ccb
    fP1 = [-0.003, 0e3]
    fP2 = [-0.0041, -2e3]
    fP3 = [-0.0041, -14e3]
    fP4 = [0.015, -17e3]
    fP5 = [0.024, -24e3]

    fric_depth1 = findall(abs.(FltX) .<= abs(fP2[2]))
    fric_depth2 = findall(abs(fP2[2]) .< abs.(FltX) .<= abs(fP3[2]))
    fric_depth3 = findall(abs(fP3[2]) .< abs.(FltX) .<= abs(fP4[2]))
    fric_depth4 = findall(abs(fP4[2]) .< abs.(FltX) .<= abs(fP5[2]))
    fric_depth5 = findall(abs.(FltX) .> abs(fP5[2]))

    a_b[fric_depth1] .= Int1D(fP1, fP2, FltX[fric_depth1])
    a_b[fric_depth2] .= Int1D(fP2, fP3, FltX[fric_depth2])
    a_b[fric_depth3] .= Int1D(fP3, fP4, FltX[fric_depth3])
    a_b[fric_depth4] .= Int1D(fP4, fP5, FltX[fric_depth4])
    a_b[fric_depth5] .= 0.0047

    #  cca[fric_depth4] .= Int1D(fP4, fP5, FltX[fric_depth4]) .+ 0.0001
    cca .= ccb .+ a_b
    #  ccb .= cca .- a_b

    return cca, ccb
end



# Effective normal stress
function SeffDepth(FltX)

    FltNglob = length(FltX)

    Seff::Array{Float64} = repeat([50e6], FltNglob)
    sP1 = [10e6 0]
    sP2 = [50e6 -2e3]
    Seff_depth = findall(abs.(FltX) .<= abs(sP2[2]))
    Seff[Seff_depth] = Int1D(sP1, sP2, FltX[Seff_depth])

    return Seff
end


# Shear stress
function tauDepth(FltX)

    FltNglob = length(FltX)

    tauo::Array{Float64} = repeat([22.5e6], FltNglob)
    tP1 = [0.01e6 0]
    tP2 = [30e6 -2e3]
    #  tP2 = [30e6 -0.5e3]
    tP3 = [30e6 -14e3]
    tP4 = [22.5e6 -17e3]
    tP5 = [22.5e6 -24e3]

    tau_depth1 = findall(abs.(FltX) .<= abs(tP2[2]))
    tau_depth2 = findall(abs(tP2[2]) .< abs.(FltX) .<= abs(tP3[2]))
    tau_depth3 = findall(abs(tP3[2]) .< abs.(FltX) .<= abs(tP4[2]))
    tau_depth4 = findall(abs(tP4[2]) .< abs.(FltX) .<= abs(tP5[2]))

    tauo[tau_depth1] = Int1D(tP1, tP2, FltX[tau_depth1])
    tauo[tau_depth2] = Int1D(tP2, tP3, FltX[tau_depth2])
    tauo[tau_depth3] = Int1D(tP3, tP4, FltX[tau_depth3])
    tauo[tau_depth4] = Int1D(tP4, tP5, FltX[tau_depth4])

    return tauo
end

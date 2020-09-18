exp1(x::Float64) = ccall(:exp, Float64, (Float64,), x)
log1(x::Float64) = ccall(:log, Float64, (Float64,), x)

# IDstate functions
function IDS!(xLf, Vo, psi, dt, Vf, cnd, IDstate = 2)
    #= compute slip-rates on fault based on different
       formulations =#

    if IDstate == 1
        psi1 = psi + dt*((Vo./xLf).*exp1(-psi) - abs(Vf)./xLf)

    elseif IDstate == 2
        VdtL = abs(Vf)*dt/xLf
        if VdtL < cnd
            psi1 = log1( exp1(psi-VdtL) + Vo*dt/xLf -
                        0.5*Vo*abs(Vf)*dt*dt/(xLf^2))
        else
            psi1 = log1(exp1(psi-VdtL) + (Vo/abs(Vf))*(1-exp1(-VdtL)))
        end

    elseif IDstate == 3
        psi1 = exp1(-abs(Vf)*dt/xLf) * log1(abs(Vf)/Vo) + 
        exp1(-abs(Vf)*dt/xLf)*psi + log1(Vo/abs(Vf))

        if ~any(imag(psi1)) == 0
            return
        end
    end

    return psi1

end

# On fault slip rates
function IDS2!(xLf, Vo, psi, psi1, dt, Vf, Vf1, IDstate = 2)
            
    if IDstate == 1
        psi2 = psi + 0.5*dt*( (Vo/xLf)*exp1(-psi) - abs(Vf)/xLf 
                                + (Vo/xLf)*exp1(-psi1) - abs(Vf1)/xLf )

    elseif IDstate == 2
        VdtL = 0.5*abs(Vf1 + Vf)*dt/xLf

        if VdtL < 1e-6
            psi2 = log1( exp1(psi-VdtL) + Vo*dt/xLf -
                            0.5*Vo*0.5*abs(Vf1 + Vf)*dt*dt/(xLf^2))
        else
            psi2 = log1(exp1(psi-VdtL) + 
                            (Vo/(0.5*abs(Vf + Vf1)))*(1-exp1(-VdtL)))
        end

    elseif IDstate == 3
        psi2 = exp1(-0.5*abs(Vf + Vf1)*dt/xLf) * log1(0.5*abs(Vf + Vf1)/Vo) + 
                        exp1(-0.5*abs(Vf + Vf1)*dt/xLf)*psi 
                        + log1(Vo/(-0.5*abs(Vf + Vf1)) )
    end

    return psi2
end

# Slip rates on fault for quasi-static regime
function slrFunc!(P::params_farray, NFBC, FltNglob, psi, psi1, Vf, Vf1, IDstate, tau1, dt)
    
    tauAB::Vector{Float64} = zeros(FltNglob)

    #  temp::Float64 = 0.

    for j = NFBC:FltNglob

        #  temp = 0.
        psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)

        tauAB[j] = tau1[j] + P.tauo[j]
        fa = tauAB[j]/(P.Seff[j]*P.cca[j])
        help = -(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j]
        help1 = exp(help + fa)
        help2 = exp(help - fa)
        Vf1[j] = P.Vo[j]*(help1 - help2)
    end

return psi1, Vf1

end

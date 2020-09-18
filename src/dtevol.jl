############################################
#   Compute the timestep for next iteration
############################################

function dtevol!(dt, dtmin, XiLf, FaultNglob, NFBC, Vf, isolver)
    
    dtmax::Int = 50 * 24 * 60*60		# 5 days
    dtincf::Float64 = 1.2

    if isolver == 1

        # initial value of dt
        dtnx = dtmax

        # Adjust the timestep according to cell velocities and slip
        for i = NFBC:FaultNglob 

            if abs(Vf[i])*dtmax > XiLf[i]
                dtcell = XiLf[i]/abs(Vf[i])

                if dtcell < dtnx
                    dtnx = dtcell
                end
            end
        end

        if dtmin > dtnx
            dtnx = dtmin
        end

        if dtnx > dtincf*dt
            dtnx = dtincf*dt
        end

        dt = dtnx

    elseif isolver == 2
        
        dt = dtmin
    end

    return dt

end

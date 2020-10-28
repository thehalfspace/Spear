###############################################################################
#
#	SPECTRAL ELEMENT METHOD FOR EARTHQUAKE CYCLE SIMULATION
#
#   Written in: Julia 1.0
#
#	Created: 09/18/2022
#   Author: Prithvi Thakur (Original code by Kaneko et al. 2011)
#
#	Adapted from Kaneko et al. (2011)
#	and J.P. Ampuero's SEMLAB
#
###############################################################################

# Healing exponential function
function healing2(t,tStart,dam)
    """ hmax: coseismic damage amplitude
        r: healing rate (0.05 => 80 years to heal completely)
                        (0.5 => 15 years to heal completely)
                        (0.8 => 8 years to heal completely)
    """
    hmax = 0.05
    r =  0.7   # 1/1.5

    hmax*(1 .- exp.(-r*(t .- tStart)/P[1].yr2sec)) .+ dam
end

function main(P)

    # P[1] = integer
    # P[2] = float
    # P[3] = float array
    # P[4] = integer array
    # P[5] = ksparse
    # P[6] = damage_idx


    #  W_orig = W[:,:,damage_idx]
    #  damage_amount::Float64 = 1.0

    # Shear modulus ratio of damage/host rock
    alphaa = 0.80

    # Time solver variables
    dt::Float64 = P[2].dt0
    dtmin::Float64 = dt
    half_dt::Float64 = 0.5*dtmin
    half_dt_sq::Float64 = 0.5*dtmin^2

    # dt modified slightly for damping
    if P[2].ETA != 0
        dt = dt/sqrt(1 + 2*P[2].ETA)
    end

    # Initialize kinematic field: global arrays
    d::Vector{Float64} = zeros(P[1].nglob)
    v::Vector{Float64} = zeros(P[1].nglob)
    v .= 0.5e-3
    a::Vector{Float64} = zeros(P[1].nglob)

    #.....................................
    # Stresses and time related variables
    #.....................................
    tau::Vector{Float64} = zeros(P[1].FltNglob)
    FaultC::Vector{Float64} = zeros(P[1].FltNglob)
    Vf::Vector{Float64} =  zeros(P[1].FltNglob)
    Vf1::Vector{Float64} = zeros(P[1].FltNglob)
    Vf2::Vector{Float64} = zeros(P[1].FltNglob)
    Vf0::Vector{Float64} = zeros(length(P[4].iFlt))
    FltVfree::Vector{Float64} = zeros(length(P[4].iFlt))
    psi::Vector{Float64} = zeros(P[1].FltNglob)
    psi0::Vector{Float64} = zeros(P[1].FltNglob)
    psi1::Vector{Float64} = zeros(P[1].FltNglob)
    psi2::Vector{Float64} = zeros(P[1].FltNglob)
    tau1::Vector{Float64} = zeros(P[1].FltNglob)
    tau2::Vector{Float64} = zeros(P[1].FltNglob)
    tau3::Vector{Float64} = zeros(P[1].FltNglob)


    # Initial state variable
    psi = P[3].tauo./(P[3].Seff.*P[3].ccb) - P[3].fo./P[3].ccb - (P[3].cca./P[3].ccb).*log.(2*v[P[4].iFlt]./P[3].Vo)
    psi0 .= psi[:]

    isolver::Int = 1

    # Some more initializations
    r::Vector{Float64} = zeros(P[1].nglob)
    beta_::Vector{Float64} = zeros(P[1].nglob)
    alpha_::Vector{Float64} = zeros(P[1].nglob)

    F::Vector{Float64} = zeros(P[1].nglob)
    dPre::Vector{Float64} = zeros(P[1].nglob)
    vPre::Vector{Float64} = zeros(P[1].nglob)
    dd::Vector{Float64} = zeros(P[1].nglob)
    dnew::Vector{Float64} = zeros(length(P[4].FltNI))


    # Save output variables at certain timesteps: define those timesteps
    tvsx::Float64 = 2e-0*P[1].yr2sec  # 2 years for interseismic period
    tvsxinc::Float64 = tvsx

    tevneinc::Float64 = 0.1    # 0.5 second for seismic period
    delfref = zeros(P[1].FltNglob)

    # Iterators
    idelevne::Int= 0
    tevneb::Float64= 0.
    tevne::Float64= 0.
    ntvsx::Int= 0
    nevne::Int= 0
    slipstart::Int= 0
    idd::Int = 0
    it_s = 0; it_e = 0
    rit = 0

    v = v[:] .- 0.5*P[2].Vpl
    Vf = 2*v[P[4].iFlt]
    iFBC::Vector{Int64} = findall(abs.(P[3].FltX) .> 24e3)
    NFBC::Int64 = length(iFBC) + 1
    Vf[iFBC] .= 0.


    v[P[4].FltIglobBC] .= 0.

    # on fault and off fault stiffness
    Ksparse = P[5]

    # Intact rock stiffness
    Korig = copy(Ksparse)   # K original

    # Linear solver stuff
    kni = -Ksparse[P[4].FltNI, P[4].FltNI]
    nKsparse = -Ksparse
    
    # algebraic multigrid preconditioner
    ml = ruge_stuben(kni)
    p = aspreconditioner(ml)
    tmp = copy(a)


    # faster matrix multiplication
     #  Ksparse = Ksparse'
     #  nKsparse = nKsparse'
     #  kni = kni'

    #  Ksparse = ThreadedMul(Ksparse)
    #  nKsparse = ThreadedMul(nKsparse)
    #  kni = ThreadedMul(kni)


    # Damage evolution stuff
    did = P[10]
    dam = alphaa

    # Save parameters to file
    open(string(out_dir,"params.out"), "w") do io
        write(io, join(P[3].Seff/1e6, " "), "\n")
        write(io, join(P[3].tauo/1e6, " "), "\n")
        write(io, join(-P[3].FltX/1e3, " "), "\n")
        write(io, join(P[3].cca, " "), "\n")
        write(io, join(P[3].ccb, " "), "\n")
        write(io, join(P[3].xLf, " "), "\n")
    end


    # Open files to begin writing
    open(string(out_dir,"stress.out"), "w") do stress
    open(string(out_dir,"sliprate.out"), "w") do sliprate
    open(string(out_dir,"slip.out"), "w") do slip
    open(string(out_dir,"delfsec.out"), "w") do dfsec
    open(string(out_dir,"delfyr.out"), "w") do dfyr
    open(string(out_dir,"event_time.out"), "w") do event_time
    open(string(out_dir,"event_stress.out"), "w") do event_stress
    open(string(out_dir,"coseismic_slip.out"), "w") do dfafter
    open(string(out_dir,"time_velocity.out"), "w") do Vf_time

    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0.
    Vfmax = 0.
    
    tStart2 = dt
    tStart = dt
    tEnd = dt
    taubefore = P[3].tauo
    tauafter = P[3].tauo
    delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t 
    hypo = 0.

    while t < P[1].Total_time
        it = it + 1
        t = t + dt

        if isolver == 1

            vPre .= v
            dPre .= d

            Vf0 .= 2*v[P[4].iFlt] .+ P[2].Vpl
            Vf  .= Vf0

            for p1 = 1:2

                # Compute the on-Fault displacement
                F .= 0.
                F[P[4].iFlt] .= dPre[P[4].iFlt] .+ v[P[4].iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[P[4].FltNI]


                # Solve d = K^-1F by MGCG
                rhs = (mul!(tmp,Ksparse,F))[P[4].FltNI]
                #  rhs = (Ksparse*F)[P[4].FltNI]

                # direct inversion
                #  dnew = -(kni\rhs)

                # mgcg
                dnew = cg!(dnew, kni, rhs, Pl=p, tol=1e-6)

                # update displacement on the medium
                d[P[4].FltNI] .= dnew

                # make d = F on the fault
                d[P[4].iFlt] .= F[P[4].iFlt]

                # Compute on-fault stress
                a .= 0.
                mul!(a,Ksparse,d)
                #   a = Ksparse*d

                # Enforce K*d to be zero for velocity boundary
                a[P[4].FltIglobBC] .= 0.

                tau1 .= -a[P[4].iFlt]./P[3].FltB

                # Function to calculate on-fault sliprate
                psi1, Vf1 = slrFunc!(P[3], NFBC, P[1].FltNglob, psi, psi1, Vf, Vf1, P[1].IDstate, tau1, dt)

                Vf1[iFBC] .= P[2].Vpl
                Vf .= (Vf0 + Vf1)/2
                v[P[4].iFlt] .= 0.5*(Vf .- P[2].Vpl)

            end

            psi .= psi1[:]
            tau .= tau1[:]
            tau[iFBC] .= 0.
            Vf1[iFBC] .= P[2].Vpl

            v[P[4].iFlt] .= 0.5*(Vf1 .- P[2].Vpl)
            v[P[4].FltNI] .= (d[P[4].FltNI] .- dPre[P[4].FltNI])/dt

            # Line 731: P_MA: Omitted
            a .= 0.
            d[P[4].FltIglobBC] .= 0.
            v[P[4].FltIglobBC] .= 0.

            #---------------
            # Healing stuff: Ignore for now
            # --------------
            if it > 3
            #  if t > 12*P[1].yr2sec
                alphaa = healing2(t, tStart2, dam)
                #  alphaa[it] = αD(t, tStart2, dam)

                for id in did
                    Ksparse[id] = alphaa*Korig[id]
                end

                #  println("alpha healing = ", alphaa[it])

                # Linear solver stuff
                kni = -Ksparse[P[4].FltNI, P[4].FltNI]
                nKsparse = -Ksparse
                # multigrid
                ml = ruge_stuben(kni)
                p = aspreconditioner(ml)

                # faster matrix multiplication
                #  Ksparse = Ksparse'
                #  nKsparse = nKsparse'
                #  kni = kni'
            end

        
        # If isolver != 1, or max slip rate is < 10^-2 m/s
        else

            dPre .= d
            vPre .= v

            # Update
            d .= d .+ dt.*v .+ (half_dt_sq).*a

            # Prediction
            v .= v .+ half_dt.*a
            a .= 0.

            # Internal forces -K*d[t+1] stored in global array 'a'
            mul!(a,nKsparse,d)
            #   a = nKsparse*d

            # Enforce K*d to be zero for velocity boundary
            a[P[4].FltIglobBC] .= 0.

            # Absorbing boundaries
            a[P[4].iBcL] .= a[P[4].iBcL] .- P[3].BcLC.*v[P[4].iBcL]
            a[P[4].iBcT] .= a[P[4].iBcT] .- P[3].BcTC.*v[P[4].iBcT]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[P[4].iFlt] .+ 2*half_dt*a[P[4].iFlt]./P[3].M[P[4].iFlt]
            Vf .= 2*vPre[P[4].iFlt] .+ P[2].Vpl


            # Sliprate and NR search
            psi1, Vf1, tau1, psi2, Vf2, tau2 = FBC!(P[1].IDstate, P[3], NFBC, P[1].FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)

            tau .= tau2 .- P[3].tauo
            tau[iFBC] .= 0.
            psi .= psi2
            a[P[4].iFlt] .= a[P[4].iFlt] .- P[3].FltB.*tau
            ########## End of fault boundary condition ##############


            # Solve for a_new
            a .= a./P[3].M

            # Correction
            v .= v .+ half_dt*a

            v[P[4].FltIglobBC] .= 0.
            a[P[4].FltIglobBC] .= 0.

            #### Line 861: Omitting P_Ma


        end # of isolver if loop

        Vfmax = 2*maximum(v[P[4].iFlt]) .+ P[2].Vpl

        #-----
        # Output the variables before and after events
        #-----
        if Vfmax > 1.01*P[2].Vthres && slipstart == 0
            it_s = it_s + 1
            delfref = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            
            slipstart = 1

            tStart = t
            taubefore = (tau +P[3].tauo)./1e6

            vhypo, indx = findmax(2*v[P[4].iFlt] .+ P[2].Vpl)
            hypo = P[3].FltX[indx]

        end
        if Vfmax < 0.99*P[2].Vthres && slipstart == 1
            it_e = it_e + 1
            delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t .- delfref
            
            tEnd = t 
            tauafter = (tau +P[3].tauo)./1e6
            
            # Save start and end time and stress
            write(event_time, join(hcat(tStart,tEnd, -hypo), " "), "\n")
            write(event_stress, join(hcat(taubefore, tauafter), " "), "\n")
            write(dfafter, join(delfafter, " "), "\n")
            
            slipstart = 0
            
            # at the end of each earthquake, the shear wave velocity in the damaged zone reduces by 10%

                # Time condition of 10 years
                #  if t > 10*P[1].yr2sec

                    #  use this for no permanent damage
                    alphaa = 0.8
                    dam = alphaa


                    #  Use this for permanent damage
                    #  alphaa = alphaa - 0.05
                    #  dam = alphaa
                    #  if dam < 0.60
                        #  alphaa = 0.60
                        #  dam = 0.60
                    #  end

                    tStart2 = t

                    for id in did
                        Ksparse[id] = alphaa*Korig[id]
                    end

                    # Linear solver stuff
                    kni = -Ksparse[P[4].FltNI, P[4].FltNI]
                    nKsparse = -Ksparse
                    # multigrid
                    ml = ruge_stuben(kni)
                    p = aspreconditioner(ml)

                #  end

                println("alphaa = ", alphaa)

            #  end

        end
        


        #-----
        # Output the variables certain timesteps: 2yr interseismic, 1 sec dynamic
        #-----
        if t > tvsx
            ntvsx = ntvsx + 1
            idd += 1
            #  write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
            write(dfyr, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")

            tvsx = tvsx + tvsxinc
        end

        if Vfmax > P[2].Vevne
            if idelevne == 0
                nevne = nevne + 1
                idd += 1
                idelevne = 1
                tevneb = t
                tevne = tevneinc

                #  write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
            end

            if idelevne == 1 && (t - tevneb) > tevne
                nevne = nevne + 1
                idd += 1

                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end

        current_sliprate = 2*v[P[4].iFlt] .+ P[2].Vpl

        # Output timestep info on screen
        if mod(it,500) == 0
            @printf("Time (yr) = %1.5g\n", t/P[1].yr2sec) 
            #  println("Vfmax = ", maximum(current_sliprate))
        end


        # Write stress, sliprate, slip to file every 10 timesteps
        if mod(it,10) == 0
            write(sliprate, join(2*v[P[4].iFlt] .+ P[2].Vpl, " "), "\n")
            write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
        end

        # Determine quasi-static or dynamic regime based on max-slip velocity
        #  if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
        if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
            isolver = 1
        else
            isolver = 2
        end

        # Write max sliprate and time
        write(Vf_time, join(hcat(t,Vfmax,Vf[end], alphaa), " "), "\n")

        # Compute next timestep dt
        dt = dtevol!(dt , dtmin, P[3].XiLf, P[1].FltNglob, NFBC, current_sliprate, isolver)


    end # end of time loop

    # close files
    end
    end
    end
    end
    end
    end
    end
    end
    end


end


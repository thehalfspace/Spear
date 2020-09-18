###############################################################################
#
#	SPECTRAL ELEMENT METHOD FOR EARTHQUAKE CYCLE SIMULATION
#
#   Written in: Julia 1.0
#
#	Created: 06/20/2018
#   Author: Prithvi Thakur (Original code by Kaneko et al. 2011)
#
#	Adapted from Kaneko et al. (2011)
#	and J.P. Ampuero's SEMLAB
#
#   CHANGELOG:
#       * 06-15-2019: Assemble stiffness as sparse matrix
#       * 04-20-2019: Implement multithreading for element calculations
#       * 02-19-2019: Remove all the shared array stuff
#       * 10-24-2018: Correct the free surface boundary condition
#       * 08-27-2018: Remove the dependency on Parameters.jl
#       * 08-26-2018: Using distributed for loop in PCG and NRsearch
#       * 08-24-2018: Create a separate function for NRsearch loop: FBC()
#       * 08-20-2018: Use JLD2 to store data instead of JLD
#       * 08-14-2018: Modify script to automatically make plots directory
#                     and save.
#
#       (Old stuff: I don't remember the dates (08/2017-08/2018))
#       * Change functions to adapt Julia 1.0 changes
#       * Move the cumulative slip calculation outside the time loop
#       * Add scripts to compute the earthquake magnitude
#       * Add plots script for various plotting functions
#       * Implemented elastic halfspace
#       * Setup for a shallow fault zone
#       * Organize everything into structs and functions
#       * Interpolation for initial stress and friction in halfspace
#       * Add separate files for parameters, initial conditions, functions
###############################################################################
#  include("$(@__DIR__)/damageEvol.jl")	    #	Set Parameters

 # Threaded matrix multiplication
 import Base: eltype, size
 #  import LinearAlgebra: A_mul_B!
 using Base.Threads

 struct ThreadedMul{Tv,Ti}
         A::SparseMatrixCSC{Tv,Ti}
 end

 function LinearAlgebra.mul!(y::AbstractVector, M::ThreadedMul, x::AbstractVector)
     @threads for i = 1 : M.A.n
          _threaded_mul!(y, M.A, x, i)
     end
      y
 end

 @inline function _threaded_mul!(y, A::SparseMatrixCSC{Tv}, x, i) where {Tv}
     s = zero(Tv)
     @inbounds for j = A.colptr[i] : A.colptr[i + 1] - 1
         s += A.nzval[j] * x[A.rowval[j]]
     end

     @inbounds y[i] = s
     y
 end
 eltype(M::ThreadedMul) = eltype(M.A)
 size(M::ThreadedMul, I...) = size(M.A, I...)

# Damage multiplier
#  d1 = P[3].vs1
#  d2 = P[3].vs2


#  α_is(t, co) = (log10(t/P[1].yr2sec)/log10(1.0e4 - t/P[1].yr2sec) + co)^2
#  α_s(t, co) = (co*(1 - 1.0e-3*t))^2


# Damage multiplier function
function αD(α0, t, tStart, co, isolver)
    if isolver == 1
        aa = (log10(t/P[1].yr2sec)/log10(1.0e4 - t/P[1].yr2sec)) + co

        if aa < α0
            return α0
        elseif aa > 1
            return 1
        else
            return aa
        end

    #  elseif isolver == 2
        #  aa = (co - 1.0e-3*(t - tStart))

        #  if aa < 0.6
            #  return 0.6
        #  elseif aa > 1
            #  return 1
        #  else
            #  return aa
        #  end
    else
        return 1
    end

end


# Save output to file dynamically
file  = jldopen("$(@__DIR__)/data/test02.jld2", "w")

function main(P)

    # P[1] = integer
    # P[2] = float
    # P[3] = float array
    # P[4] = integer array
    # P[5] = ksparse
    # P[6] = damage_idx


    #  W_orig = W[:,:,damage_idx]
    #  damage_amount::Float64 = 1.0


    #  wgll2::Array{Float64,2} = S.wgll*S.wgll';

    # Some Stuff
    res = 2
    LX::Int = 48e3  # depth dimension of rectangular domain
    LY::Int = 30e3 # off fault dimenstion of rectangular domain

    NelX::Int = 30*res # no. of elements in x
    NelY::Int = 20*res # no. of elements in y

    dxe::Float64 = LX/NelX #	Size of one element along X
    dye::Float64 = LY/NelY #	Size of one element along Y


    rho1 = 2670
    vs1 = 3464

    ThickX = LX - 24e3
    ThickY::Float64 = ceil(5e3/dye)*dye   # ~ 0.25*2 km wide

    iglob = P[6]
    ngll = P[7]
    wgll2 = P[8]
    nglob = P[9]


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

    # Skip lines 486-490
    # Skip lines 492-507: Outloc1, 2, variables.

    # Some more initializations
    r::Vector{Float64} = zeros(P[1].nglob)
    beta_::Vector{Float64} = zeros(P[1].nglob)
    alpha_::Vector{Float64} = zeros(P[1].nglob)

    F::Vector{Float64} = zeros(P[1].nglob)
    dPre::Vector{Float64} = zeros(P[1].nglob)
    vPre::Vector{Float64} = zeros(P[1].nglob)
    dd::Vector{Float64} = zeros(P[1].nglob)
    dnew::Vector{Float64} = zeros(length(P[4].FltNI))

    # Preallocate variables with unknown size
    #  seismic_stress, seismic_slipvel, seismic_slip
    #  index_eq
    #  is_stress, is_slipvel, is_slip
    #  dSeis, vSeis, aSeis
    #  tStart, tEnd
    #  taubefore, tauafter, delfafter
    #  hypo, time_, Vfmax
    nseis = length(P[4].out_seis)

    output = results(zeros(P[1].FltNglob, 8000), zeros(P[1].FltNglob, 8000),
                     zeros(P[1].FltNglob, 8000),
                     zeros(10000),
                     zeros(P[1].FltNglob, 2000), zeros(P[1].FltNglob, 2000),
                     zeros(P[1].FltNglob, 2000),
                     zeros(80000,nseis), zeros(80000,nseis), zeros(80000,nseis),
                     zeros(400), zeros(400),
                     zeros(P[1].FltNglob, 400), zeros(P[1].FltNglob, 400),
                     zeros(P[1].FltNglob, 400), zeros(400), zeros(700000),
                     zeros(700000))

    # Save output variables at certain timesteps: define those timesteps
    tvsx::Float64 = 2*P[1].yr2sec  # 2 years for interseismic period
    tvsxinc::Float64 = tvsx

    tevneinc::Float64 = 0.5    # 0.5 second for seismic period
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
    # Indices of Stiffness matrix corresponding to the damaged zone
    #  Ksparse = rcmpermute(P[5])

    # Initial intact rigidity ratio b/w host rock and damage
    co = 1

    # Linear solver stuff
    kni = -Ksparse[P[4].FltNI, P[4].FltNI]
    nKsparse = -Ksparse
    # multigrid
    ml = ruge_stuben(kni)
    p = aspreconditioner(ml)
    tmp = copy(a)

    # faster matrix multiplication
    #   Ksparse = Ksparse'
    #   nKsparse = nKsparse'
    #   kni = kni'

    Ksparse = ThreadedMul(Ksparse)
    nKsparse = ThreadedMul(nKsparse)
    kni = ThreadedMul(kni)

    # Temporary Debugging variables: CLEAN UP LATER
    alphaa = ones(1000000)

    tStart = dt


    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0.

    while t < P[1].Total_time
        it = it + 1
        t = t + dt

        if it == 1
            #  alphaa[it] = 1
        end

        output.time_[it] = t

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


        #----
        # Output variables at different depths for every timestep
        # Omitted the part of code from line 871 - 890, because I
        # want to output only certain variables each timestep
        #----


        #-----
        # Output the variables before and after events
        #-----
        if Vfmax > 1.01*P[2].Vthres && slipstart == 0
            it_s = it_s + 1
            delfref = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            slipstart = 1
            output.tStart[it_s] = output.time_[it]
            tStart = output.time_[it]
            output.taubefore[:,it_s] = (tau +P[3].tauo)./1e6
            vhypo, indx = findmax(2*v[P[4].iFlt] .+ P[2].Vpl)
            output.hypo[it_s] = P[3].FltX[indx]

        end
        if Vfmax < 0.99*P[2].Vthres && slipstart == 1
            it_e = it_e + 1
            output.delfafter[:,it_e] = 2*d[P[4].iFlt] .+ P[2].Vpl*t .- delfref
            output.tauafter[:,it_e] = (tau + P[3].tauo)./1e6
            output.tEnd[it_e] = output.time_[it]
            slipstart = 0

        end
        #-----
        # Output the variables certain timesteps: 2yr interseismic, 1 sec dynamic
        #-----
        if output.time_[it] > tvsx
            ntvsx = ntvsx + 1
            idd += 1
            output.is_slip[:,ntvsx] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            output.is_slipvel[:,ntvsx] = 2*v[P[4].iFlt] .+ P[2].Vpl
            output.is_stress[:,ntvsx] = (tau + P[3].tauo)./1e6
            output.index_eq[idd] = 1

            tvsx = tvsx + tvsxinc
        end

        if Vfmax > P[2].Vevne
            if idelevne == 0
                nevne = nevne + 1
                idd += 1
                idelevne = 1
                tevneb = t
                tevne = tevneinc

                output.seismic_slip[:,nevne] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
                output.seismic_slipvel[:,nevne] = 2*v[P[4].iFlt] .+ P[2].Vpl
                output.seismic_stress[:,nevne] = (tau + P[3].tauo)./1e6
                output.index_eq[idd] = 2
            end

            if idelevne == 1 && (t - tevneb) > tevne
                nevne = nevne + 1
                idd += 1

                output.seismic_slip[:,nevne] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
                output.seismic_slipvel[:,nevne] = 2*v[P[4].iFlt] .+ P[2].Vpl
                output.seismic_stress[:,nevne] = (tau + P[3].tauo)./1e6
                output.index_eq[idd] = 2
                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end

        # Output timestep info on screen
        if mod(it,500) == 0
            @printf("Time (yr) = %1.5g\n", t/P[1].yr2sec)
        end


        # output variables at prescribed locations every 10 timesteps
        if mod(it,10) == 0
            rit += 1
            output.dSeis[rit,:] = d[P[4].out_seis]
            output.vSeis[rit,:] = v[P[4].out_seis]
            output.aSeis[rit,:] = a[P[4].out_seis]
        end

        # Determine quasi-static or dynamic regime based on max-slip velocity
        if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
            isolver = 1
        else
            isolver = 2
        end

        output.Vfmax[it] = Vfmax

        current_sliprate = 2*v[P[4].iFlt] .+ P[2].Vpl

        # Compute next timestep dt
        dt = dtevol!(dt , dtmin, P[3].XiLf, P[1].FltNglob, NFBC, current_sliprate, isolver)

    end # end of time loop

    # Remove zeros from preallocated vectors
    output.seismic_stress   = output.seismic_stress[:,1:nevne]
    output.seismic_slipvel  = output.seismic_slipvel[:,1:nevne]
    output.seismic_slip     = output.seismic_slip[:,1:nevne]
    output.index_eq         = output.index_eq[1:idd]
    output.is_stress        = output.is_stress[:,1:ntvsx]
    output.is_slipvel       = output.is_slipvel[:,1:ntvsx]
    output.is_slip          = output.is_slip[:,1:ntvsx]
    output.dSeis            = output.dSeis[1:rit,:]
    output.vSeis            = output.vSeis[1:rit,:]
    output.aSeis            = output.aSeis[1:rit,:]
    output.tStart           = output.tStart[1:it_s]
    output.tEnd             = output.tEnd[1:it_e]
    output.taubefore        = output.taubefore[:,1:it_s]
    output.tauafter         = output.tauafter[:,1:it_e]
    output.delfafter        = output.delfafter[:,1:it_e]
    output.hypo             = output.hypo[1:it_s]
    output.time_            = output.time_[1:it]
    output.Vfmax            = output.Vfmax[1:it]
    alphaa                   = alphaa[1:it]

    #  return d, v, a, 2*v[P[4].iFlt] .+ P[2].Vpl
    return output, alphaa
    file["O"] = output

    close(file)

end

mutable struct results
    seismic_stress::Matrix{Float64}
    seismic_slipvel::Matrix{Float64}
    seismic_slip::Matrix{Float64}
    index_eq::Vector{Float64}
    is_stress::Matrix{Float64}
    is_slipvel::Matrix{Float64}
    is_slip::Matrix{Float64}
    dSeis::Matrix{Float64}
    vSeis::Matrix{Float64}
    aSeis::Matrix{Float64}
    tStart::Vector{Float64}
    tEnd::Vector{Float64}
    taubefore::Matrix{Float64}
    tauafter::Matrix{Float64}
    delfafter::Matrix{Float64}
    hypo::Vector{Float64}
    time_::Vector{Float64}
    Vfmax::Vector{Float64}
end

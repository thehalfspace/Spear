#################################
# READ OUTPUT FROM SIMULATION
#################################

mutable struct results
    seismic_stress::Array{Float64,2}
    seismic_slipvel::Array{Float64,2}
    seismic_slip::Array{Float64,2}
    index_eq::Array{Float64}
    is_stress::Array{Float64,2}
    is_slipvel::Array{Float64,2}
    is_slip::Array{Float64,2}
    dSeis::Matrix{Float64}
    vSeis::Matrix{Float64}
    aSeis::Matrix{Float64}
    tStart::Array{Float64}
    tEnd::Array{Float64}
    taubefore::Array{Float64,2}
    tauafter::Array{Float64,2}
    delfafter::Array{Float64,2}
    hypo::Array{Float64}
    time_::Array{Float64}
    Vfmax::Array{Float64}
end

struct params_int{T<:Int}
    # Domain size
    Nel::T
    FltNglob::T

    # Time parameters
    yr2sec::T
    Total_time::T
    IDstate::T
    
    # Fault setup parameters
    nglob::T

end

struct params_float{T<:AbstractFloat}
    # Jacobian for global -> local coordinate conversion
    #  jac::T
    #  coefint1::T
    #  coefint2::T
    ETA::T

    # Earthquake parameters
    Vpl::T
    Vthres::T
    Vevne::T

    # Setup parameters
    dt0::T
end

struct params_farray{T<:Array{Float64}}
    fo::T
    Vo::T
    xLf::T
    
    M::T

    BcLC::T
    BcTC::T

    FltB::T
    FltZ::T
    FltX::T

    cca::T
    ccb::T
    Seff::T
    tauo::T
    XiLf::T
    #  diagKnew::T

    xout::T
    yout::T
end

struct params_iarray{T<:Array{Int}}
    iFlt::T
    iBcL::T
    iBcT::T
    FltIglobBC::T
    FltNI::T
    out_seis::T
end

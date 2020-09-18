#######################################################
# CALCULATE CUMULATIVE SLIP AT CERTAIN TIME INTERVALS
#######################################################

function cumSlip(Slip, SlipVel, time_)
    # Output stress, slip, sliprate on fault every certain interval

    yr2sec = 365*24*60*60 
    tvsx = 2*yr2sec
    tvsxinc = tvsx

    tevneinc = 1
    Vevne = 0.001

    ntvsx = 0
    nevne = 0
    idelevne = 0
    tevneb = 0
    tevne = 0
    
    Vfmax = maximum(SlipVel, dims = 1)[:]
    it = length(time_)

    FltNglob = length(SlipVel[:,1])
    
    delfsec::Array{Float64} = zeros(FltNglob, it)
    delf5yr::Array{Float64} = zeros(FltNglob, it)
  
    for i = 1:it
        if time_[i] > tvsx
            ntvsx = ntvsx + 1
            
            delf5yr[:,ntvsx] = Slip[:,i]
            #Vf5yr[:,ntvsx] = 2*v[S.iFlt] .+ P.Vpl
            #Tau5yr[:,ntvsx] = (tau + tauo)./1e6
            
            tvsx = tvsx +tvsxinc
        end
        
        if Vfmax[i] > Vevne 
            if idelevne == 0
                nevne = nevne + 1
                idelevne = 1
                tevneb = time_[i]
                tevne = tevneinc

                delfsec[:,nevne] = Slip[:,i]
                #Vfsec[:,nevne] = 2*v[S.iFlt] .+ P.Vpl
                #Tausec[:,nevne] = (tau + tauo)./1e6
            end

            if idelevne == 1 && (time_[i] - tevneb) > tevne
                nevne = nevne + 1
                
                delfsec[:,nevne] = Slip[:,i]
                #Vfsec[:,nevne] = 2*v[S.iFlt] .+ P.Vpl
                #Tausec[:,nevne] = (tau + tauo)./1e6

                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end
    end 

    delfsec = delfsec[:,1:nevne]
    delf5yr = delf5yr[:,1:ntvsx]
    return delfsec, delf5yr
end

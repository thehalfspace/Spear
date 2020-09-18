#####################################
# DAMAGE EVOLUTION IN TIME
#####################################

function damage_indx!(ThickX, ThickY, Wid, dxe, dye, NGLL, NelX, NelY, iglob)

    ww::Matrix{Float64} = zeros(NGLL, NGLL)
    Ke2::Array{Float64,4} = zeros(NGLL,NGLL,NGLL,NGLL)
    Ke::Array{Float64,3} = zeros(NGLL*NGLL,NGLL*NGLL, Nel)

    @inbounds for eo in 1:Nel
        Ke2 .= 0.

        for i in 1:NGLL, j in 1:NGLL
            for k in 1:NGLL, l in 1:NGLL
                Ke2[i,j,k,l] = 1
            end
        end
        Ke[:,:,eo] = reshape(Ke2,NGLL*NGLL,NGLL*NGLL)
    end

    return FEsparse(NelX*NelY, Ke, iglob)

end

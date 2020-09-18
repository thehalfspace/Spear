# Function to find the mesh node that is closest to the requested location

function FindNearestNode(xin, yin, X, Y)
    nseis = length(xin)
    dist = zeros(nseis)
    iglob::Vector{Int} = zeros(nseis)

    for k = 1:nseis
        iglob[k] = argmin( (X .- xin[k]).^2 + (Y .- yin[k]).^2 )
        dist[k] = minimum( (X .- xin[k]).^2 + (Y .- yin[k]).^2 )
    end

    dist = sqrt.(dist)
    xout = X[iglob]
    yout = Y[iglob]

    return xout, yout, iglob, dist

end

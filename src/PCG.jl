################################################
#                                              
#   SOLVE FOR DISPLACEMENT USING PRECONDITIONED 
#           CONJUGATE GRADIENT METHOD          
#                                              
################################################

function PCG!(P::params_float, Nel::Int, diagKnew::Array{Float64}, dnew::Array{Float64}, F::Array{Float64}, iFlt::Array{Int},FltNI::Array{Int}, H::Array{Float64,2}, Ht::Array{Float64,2}, iglob::Array{Int,3}, nglob::Int, W::Array{Float64,3}, a_elem::Array{Float64}, Conn)
    
    a_local::Array{Float64} = zeros(nglob)
    dd_local::Array{Float64} = zeros(nglob)
    p_local::Array{Float64} = zeros(nglob)
   
    a_elem = element_computation!(P, iglob, F, H, Ht, W, Nel)
    Fnew = -mul!(a_local, Conn, a_elem[:])[FltNI]
    
    dd_local[FltNI] .= dnew
    dd_local[iFlt] .= 0.

    a_local[:] .= 0.
    
    a_elem = element_computation!(P, iglob, dd_local, H, Ht, W, Nel)
        
    anew = mul!(a_local, Conn, a_elem[:])[FltNI]

    # Initial residue
    rnew = Fnew - anew
    znew = rnew./diagKnew
    pnew = znew
    p_local[:] .= 0.
    p_local[FltNI] = pnew

    @inbounds for n = 1:8000
        anew[:] .= 0.
        a_local[:] .= 0.
        
        a_elem = element_computation!(P, iglob, p_local, H, Ht, W, Nel)
        anew = mul!(a_local, Conn, a_elem[:])[FltNI]

        alpha_ = znew'*rnew/(pnew'*anew)
        dnew  .+= alpha_*pnew
        rold = rnew
        zold = znew
        rnew = rold - alpha_*anew
        znew = rnew./diagKnew
        beta_ = znew'*rnew/(zold'*rold)
        pnew = znew + beta_*pnew
        p_local[:] .= 0.
        p_local[FltNI] = pnew

        if norm(rnew)/norm(Fnew) < 1e-5
            break;
        end

        if n == 8000 || norm(rnew)/norm(Fnew) > 1e10
            print(norm(rnew)/norm(Fnew))
            println("\nn = ", n)

            #filename = string(dir, "/data", name, "pcgfail.jld2")
            #@save filename dnew rnew Fnew
            @error("PCG did not converge")
            return
        end
    end

    return dnew
end


# Multi-threading
function element_computation!(P::params_float, iglob::Array{Int,3}, F_local::Array{Float64}, H::Array{Float64,2}, Ht::Array{Float64,2}, W::Array{Float64,3}, Nel)
    a_local = zeros(size(F_local))
    a_elem = zeros(size(iglob))
    Threads.@threads for tid in 1:Threads.nthreads()
        len = div(Nel, Threads.nthreads())
        domain = ((tid-1)*len + 1):tid*len

        @inbounds @simd for eo in domain
            ig = iglob[:,:,eo]
            Wlocal = W[:,:,eo]
            locall = F_local[ig]
            a_elem[:,:,eo] =  P.coefint1*H*(Wlocal.*(Ht*locall)) + P.coefint2*(Wlocal.*(locall*H))*Ht
        end
    end
    return a_elem
end


# Assembly of global stiffness matrix as a sparse matrix

function stiffness_assembly(NGLL, NelX, NelY, dxe,dye, nglob, iglob, W) 
    xgll, wgll, H::SMatrix{NGLL,NGLL,Float64} = GetGLL(NGLL)
    wgll2::SMatrix{NGLL,NGLL,Float64} = wgll*wgll'

    
    #  ig::Matrix{Int64} = zeros(NGLL,NGLL)  # iterator

    #  W = material_properties(NelX, NelY,NGLL,dxe, dye, ThickX, ThickY, wgll2, rho1, rho2, vs1, vs2)
    Ke = K_element(W, dxe, dye, NGLL, H, NelX*NelY)
    #  Ks22 = assembley(Ke, iglob, NelX*NelY, nglob)
    K = FEsparse(NelX*NelY, Ke, iglob)


    return dropzeros!(K)
    #  return rcmpermute(dropzeros!(K))
end

function FEsparse(Nel, Ke, iglob)
    K = SparseMatrixCOO()
    # using FEMSparse
    for eo in 1:Nel
        FEMSparse.assemble_local_matrix!(K, vec(iglob[:,:,eo]),
                                         vec(iglob[:,:,eo]), Ke[:,:,eo])
    end

    SparseMatrixCSC(K)
end

function K_element(W, dxe, dye, NGLL, H, Nel)
    # Jacobians
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta

    ww::Matrix{Float64} = zeros(NGLL, NGLL)
    Ke2::Array{Float64,4} = zeros(NGLL,NGLL,NGLL,NGLL)
    Ke::Array{Float64,3} = zeros(NGLL*NGLL,NGLL*NGLL, Nel)
    #  ig::Matrix{Int64} = zeros(NGLL,NGLL)  # iterator
    
    #  term1::Float64 = 0.; term2::Float64 = 0.
    del = Matrix{Float64}(I,NGLL,NGLL)  # identity matrix

        @inbounds for eo in 1:Nel
            Ke2 .= 0.

            ww = W[:,:,eo]
            term1 = 0.; term2 = 0.
            for i in 1:NGLL, j in 1:NGLL
                for k in 1:NGLL, l in 1:NGLL
                    term1 = 0.; term2 = 0.
                    for p in 1:NGLL
                        term1 += del[i,k]*ww[k,p]*(jac/dy_deta^2)*H[j,p]*H[l,p]
                        term2 += del[j,l]*ww[p,j]*(jac/dx_dxi^2)*H[i,p]*H[k,p]
                    end
                    Ke2[i,j,k,l] = term1 + term2
                end
            end
            Ke[:,:,eo] = reshape(Ke2,NGLL*NGLL,NGLL*NGLL)
        end
    Ke
end

# My naive approach
#  function assembley_2(Ke, iglob, Nel, nglob)
    #  Ksparse::SparseMatrixCSC{Float64} = spzeros(nglob,nglob) 
    #  for eo in 1:Nel
        #  ig = iglob[:,:,eo]
        #  Ksparse[vec(ig),vec(ig)] += Ke[:,:,eo]
    #  end

    #  Ksparse
#  end


# faster assembly apprach: just as fast as FESparse
#  function assembley(Ke, iglob, Nel, nglob)
    #  #  Ksparse::SparseMatrixCSC{Float64} = spzeros(nglob,nglob) 
    #  I = Vector{Int}(undef, length(Ke))
    #  J = Vector{Int}(undef, length(Ke))
    #  V = Vector{Float64}(undef, length(Ke))
    #  ct = 1

    #  for eo in 1:Nel
        #  v = view(iglob,:,:,eo)
        #  #  v = iglob[:,:,eo][:]
        #  for j in 1:length(v)
            #  for i in 1:length(v)
                #  I[ct] = v[i]
                #  J[ct] = v[j]
                #  V[ct] = Ke[i, j, eo]
                #  ct += 1
                #  #  Ksparse[vec(ig),vec(ig)] += Ke[:,:,eo]
            #  end
        #  end
    #  end

    #  return sparse(I,J,V,nglob,nglob,+)
#  end


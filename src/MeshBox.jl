###############################################################################
#		Spectral Element Mesh for Rectangular Box, with internal 
#		Gauss-Legendre-Lobatto (GLL) dub-grids.
#
#	INPUT:	LX = x-dimension
#			LY = y-dimension
#			NELX = no. of elements in x
#			NELY = no. of elements in y
#			NGLL = no. of GLL nodes


#	OUTPUT:	iglob[ngll, ngll, NELX, NELY] = maps local to global
#											numbering
#			I = iglob[i,j,e] is the global node index of the 
#							 (i,j)th GLL node internal to the
#							 e-th element.

#			Elements are numbered row by row from from bottom-left
#			to top-right. The table iglob is needed to assemble 
#			global data from local data.

#			x[:] = global x coordinates of GLL nodes, starting at 0
#			y[:] = global y coordinates of GLL nodes, starting at 0
###############################################################################


function MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)

	XGLL = GetGLL(NGLL)[1]

	iglob = zeros(Int, NGLL, NGLL, Nel)
	nglob = FltNglob*(NelY*(NGLL-1) + 1)

    x::Vector{Float64} = zeros(nglob)
    y::Vector{Float64} = zeros(nglob)

	et = 0
	last_iglob = 0

	ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)
	igL = reshape(collect(1:NGLL*(NGLL-1)), NGLL-1, NGLL) # Left edge
	igB = reshape(collect(1:NGLL*(NGLL-1)), NGLL, NGLL-1) # Bottom edge
	igLB = reshape(collect(1:(NGLL-1)*(NGLL-1)), NGLL-1, NGLL-1) # rest of the elements

	xgll = repeat(0.5*(1 .+ XGLL), 1, NGLL)
	ygll = dye*xgll'
	xgll = dxe*xgll


	@inbounds for ey = 1:NelY         # number of x elements
		@inbounds for ex = 1:NelX     # number of y elements

			et = et + 1

			# Redundant nodes at element edges
            
            # NGLL = number of GLL nodes per element
            
			if et == 1
				ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)
			else
				if ey ==1 # Bottom Row
					ig[1,:] = iglob[NGLL, :, et-1]   # Left edge
					ig[2:end, :] = last_iglob .+ igL # The rest

				elseif ex == 1	# Left Column
					ig[:,1] = iglob[:,NGLL,et-NelX]	# Bottom edge
					ig[:,2:end] = last_iglob .+ igB 	# The rest

				else 			# Other Elements
					ig[1,:] = iglob[NGLL, :, et-1]	# Left edge
					ig[:,1] = iglob[:, NGLL, et-NelX]# Bottom edge
					ig[2:end, 2:end] = last_iglob .+ igLB
				end
			end

			iglob[:,:,et] = ig
			last_iglob = ig[NGLL, NGLL]

			# Global coordinates of computational nodes
			@inbounds x[ig] .= dxe*(ex-1) .+ xgll
			@inbounds y[ig] .= dye*(ey-1) .+ ygll

		end
	end

	return iglob, x, y
end


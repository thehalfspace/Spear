########################################################
#
#	FUNCTION TO BUILD THE 2D BOUNDARIES FOR SEM
#
########################################################

function BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, 
                         dy_deta, dx_dxi, wgll, iglob, side)

	# INPUT: 
	#		wgll = GLL weights (see GetGLL)
	#		NelX, NelY = no. of elements
	#		iglob = global to local index
	#		jac1D = line Jacobian
	#		side = 'L', 'R', 'T', 'B'


	if side == 'L'
		eB = collect(0:NelY-1)*NelX .+ 1
		igll = 1
		jgll = collect(1:NGLL)
        jac1D = dy_deta
        impedance = rho1*vs1

	elseif side == 'R'
		eB = collect(0:NelY-1)*NelX .+ NelX
		igll = NGLL
		jgll = collect(1:NGLL)
        jac1D = dy_deta
        impedance = rho1*vs1

	elseif side == 'T'
		eB = (NelY-1)*NelX .+ collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = NGLL
        jac1D = dx_dxi
        impedance = rho1*vs1

	else 
		eB = collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = 1
        jac1D = dx_dxi
        impedance = 1
	end

	NelB = length(eB)
	ng = NelB*(NGLL-1) .+ 1
	iB = zeros(Int, ng)
	B = zeros(ng)
	jB = zeros(NGLL, NelB)

	for e=1:NelB
		ip = (NGLL-1)*(e-1) .+ collect(1:NGLL)
		iB[ip] = iglob[igll, jgll, eB[e]]
		jB[:,e] = ip
		B[ip] .+= jac1D*wgll*impedance
	end

	return B, iB
end

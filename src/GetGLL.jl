##################################################################
# This function obtains the weights and derivatives of 
# Gauss-Lobatto-Legendre grid points from external file.
# 
#
# INPUT:
# 		NGLL = P+1: P is the polynomial degree for interpolation
#
#
# OUTPUT:
#		x[n]: coordinates of GLL points in the reference segment [-1,1]
#		w[n]: quadrature weights
#		h[n,n]: derivatives of Lagrange polynomials at the GLL nodes
#		h[i,j] = h'_i (x[j])
using DelimitedFiles

function GetGLL(ngll)

	a = ngll
	data = readdlm("$(@__DIR__)/gll_xwh/gll_$(lpad(a, 2, "0")).tab", header=false)

	x = data[1,:]
	w = data[2,:]
	h = data[3:end,:]

	return x, w, h

end

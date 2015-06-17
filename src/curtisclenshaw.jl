module CurtisClenshaw
const libwget 	= :_Z7w_get_lPdiiPiS0_S_iS_S_S_
const libwgetinv=:_Z11w_get_inv_lPdiiPiS0_S_iS_S_S_
const libinterp = :_Z15sparse_interp_lPdiS_iiPiS0_S_iS_S_S_

Mi(i::Int) = (i==1) ? 1 : 2^(i-1)+1

function dMi(i::Int64)
	if (i==1)
		return 1
	elseif (i==2)
		return 2
	else
		return 2^(i-2)
	end
end

xi(i::Int,j::Int) = (i==1) ? 0.5 : (j-1)/(Mi(i)-1.0)


function dxi(i::Int64,j::Int64)
	if (i==1)
		return 0.5
	elseif (i==2)
		if (j==1)
			return 0.0;
		else
			return 1.0;
		end
	else
		return ((j)*2.0-1.0)/(Mi(i)-1.0)
	end
end

Xi(i::Int) = (i==1) ? [0.5] : collect(linspace(0,1,Mi(i)))


function dXi(i::Int64)
	dM 	= dMi(i)
	M  	= Mi(i)
	X  	= Array(Float64,dM)
	if (i==1)
		X = 0.5
	elseif (i==2)
		X = [0.0,1.0]
	else
		for ii = 1:dM
			X[ii] = ((ii)*2.0-1.0)/(M-1.0)
		end
	end
	return X
end

function getind(grid::Array{Float64},q::Int)
	dim = size(grid,2)
	nG = size(grid,1)
	ind = fill!(Array(Int64,nG,dim),0)#zeros(nG,dim)
	for i = 1:nG*dim
		ind[i]=q
		for ddi=q+dim:-1:2
			if mod(grid[i],1/(Mi(ddi)-1))==0.0
				ind[i] = ddi
			end
		end
		if grid[i]==0.5
			ind[i] = 1
		end
	end
	return ind
end


function basis_func(x::Float64,xij::Float64,mi::Int32)
	if (mi==1)
		return 1.0
	elseif (abs(x-xij)<(1.0/(mi-1.0)))
		return (1.0-(mi-1.0)*abs(x-xij))
	else
		return 0.0
	end
end

include("all.jl")


end



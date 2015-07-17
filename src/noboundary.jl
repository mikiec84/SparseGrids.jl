module NoBoundary
@windows ? (begin
const libwget 	= :w_get_nb
const libwgetinv= :w_get_inv_nb
const libinterp = :sparse_interp_nb
end
:begin
const libwget 	= :_Z8w_get_nbPdiiPiS0_S_iS_S_S_
const libwgetinv= :_Z12w_get_inv_nbPdiiPiS0_S_iS_S_S_
const libinterp = :_Z16sparse_interp_nbPdiS_iiPiS0_S_iS_S_S_
end)

Mi(i::Int) = 2^i-1
dMi(i::Int) = Mi(i)-Mi(i-1)
xi(i::Int,j::Int) = j/(Mi(i)+1)
dxi(i::Int,j::Int) = dXi(i)[j]
Xi(i::Int64)=[xi(i,j) for j = 1:Mi(i)]
dXi(i::Int) = setdiff(Xi(i),Xi(i-1))

function getind(grid::Array{Float64},q::Int64)
	dim = size(grid,2)
	nG = size(grid,1)
	ind = zeros(Int,nG,dim)
	for i = 1:nG*dim
		ind[i]=q
		for ddi=q+dim:-1:2
			if mod(grid[i],1/(Mi(ddi)+1))==0.0
				ind[i] = ddi
			end
		end
		if grid[i]==0.5
			ind[i] = 1
		end
	end
	return ind
end

nbfirstXi(mi) = 1/(mi+1);
nblastXi(mi) = mi/(mi+1);

function basis_func(x::Float64, xij::Float64, mi)
    if mi==1
        return 1.0
    else
        if xij==nbfirstXi(mi)

            if (x<2./(mi+1.0))
                return 2.-(mi+1.)*x;
            else
                return 0.0
            end
        end
        if xij==nblastXi(mi)

            if (x>(mi-1.)/(mi+1.))
                return (mi+1.)*x -mi + 1.0
            else
                return 0.0
            end
        else
            if (abs(x-xij)<1/(mi+1.))
                return 1.-(mi+1.)*abs(x-xij)
            else
                return 0.0
            end
        end
    end
end


include("all.jl")


end

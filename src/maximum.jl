module Maximum
@windows ? (begin
const libwget 	 = :w_get_max
const libwgetinv = :w_get_inv_max
const libinterp  = :sparse_interp_max
end
:begin
const libwget 	= :_Z9w_get_maxPdiiPiS0_S_iS_S_S_
const libwgetinv = :_Z13w_get_inv_maxPdiiPiS0_S_iS_S_S_
const libinterp = :_Z17sparse_interp_maxPdiS_iiPiS0_S_iS_S_S_
end)


Mi(i::Int) = if i>0 2^i+1 else 0 end
dMi(i::Int) = Mi(i)-Mi(i-1)
xi(i::Int,j::Int) = (j-1)/(Mi(i)-1)
dxi(i::Int,j::Int) = dXi(i)[j]
Xi(i::Int64)=[xi(i,j) for j = 1:Mi(i)]
dXi(i::Int) = setdiff(Xi(i),Xi(i-1))

function getind(grid::Array{Float64},q::Int64)
    ugrid = sort(unique(grid))
    uind = zeros(Int,size(ugrid))
    ind = zeros(Int,size(grid))
    for i = 1:length(ugrid)
        for q1 = q+size(grid,2):-1:1
            if in(ugrid[i],dXi(q1))
                uind[i]=q1
            end
        end
    end
    for i = 1:length(ugrid)
        ind[grid.==ugrid[i]]=uind[i]
    end
    return ind
end

function basis_func(x::Float64, xij::Float64, mi)
    if abs(x-xij)<1/(mi-1)
        return 1.0-(mi-1)*abs(x-xij)
    else
        return 0.0
    end
end


include("all.jl")


end


function XtoX(X::Array{Float64,1},bnds1::Array{Float64,1},bnds2::Array{Float64,1})
	lb1 = minimum(bnds1)
	ub1 = maximum(bnds1)
	lb2 = minimum(bnds2)
	ub2 = maximum(bnds2)
	out = Array(Float64,length(X))
	for i = 1:length(X)
		out[i] = (X[i]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
	end
	return out
end

function nXtoX(X::Array{Float64,2},bnds1::Array{Float64,2},bnds2::Array{Float64,2})
	N,R = size(X)
	out = Array(Float64,size(X))
	for r = 1:R
		lb1 = minimum(bnds1[:,r])
		ub1 = maximum(bnds1[:,r])
		lb2 = minimum(bnds2[:,r])
		ub2 = maximum(bnds2[:,r])
		for i = 1:N
			out[i,r] = (X[i,r]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
		end
	end
	return out
end
nXtoX(X::Array{Float64,1},bnds1::Array{Float64,2},bnds2::Array{Float64,2})=XtoX(X,bnds1[:],bnds2[:])

UtoX(U::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(U,[0.0,1.0],bnds)
XtoU(X::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(X,bnds,[0.0,1.0])

nUtoX(U::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(U,[zeros(1,size(U,2));ones(1,size(U,2))],bnds)
nXtoU(X::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(X,bnds,[zeros(1,size(X,2));ones(1,size(X,2))])
nXtoU(X::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(X,bnds[:],[0.0,1.0])
nUtoX(U::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(U,[0.0,1.0],bnds[:])


function hsh(x::Vector{Int})
	h = x[1]*17
	for i = 2:length(x)
		h+=x[i]
		h*=17
	end
	h
end

function hsh(x::Array{Int,2})
	D = size(x,2)
	H = zeros(Int,size(x,1))
	for i = 1:size(x,1)
		H[i] = x[i,1]*17
		for j = 2:D
			H[i] += x[i,j]
			H[i] *= 17
		end
	end
	return H
end


function swap!(x::Vector{Int},i::Int,j::Int)
	ta 	= x[i]
	x[i]=x[j]
	x[j]=ta
	return nothing
end

function next_permute!(seq::Vector{Int})
	n = length(seq)
	j = n-1
	while seq[j]>=seq[j+1]
		j-=1
		if j==0
			return true
		end
	end
	l = n
	while seq[j]>=seq[l]
		l-=1
	end
	swap!(seq,j,l)

	seq[j+1:n]=reverse(seq[j+1:n])
	return false
end

function get_perms(seq::Vector{Int})
	a = sort(seq)
	out = Array(Vector{Int},1)
	out[1] = [a;]
	if length(seq)==1
		return out
	end
	done = false
	while !done
		done=next_permute!(a)
		push!(out,[a;])
	end
	if out[end]==out[end-1]
		out = out[1:end-1]
	end
	return out::Vector{Vector{Int}}
end

function comb(D::Int,Q::Int)
	if D ==1
		return [Q]
	end
	Pq = Combinatorics.integer_partitions(Q)
	out = Array(Vector{Int},0)

	for p in Pq
		if length(p)==D
			pPq=get_perms(p)
			for j = 1:length(pPq)
				push!(out,pPq[j])
			end
		end
	end

	return out
end


ndgrid(v::AbstractVector) = (copy(v),)

function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid{T}(vs::AbstractVector{T}...)
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array(T, sz),n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

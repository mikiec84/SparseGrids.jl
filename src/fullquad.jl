module FullQuad
module ModGrid
using Grid
end
ndgrid(v::AbstractVector) = copy(v)

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

include("conversion.jl")
Mi(i::Int) = (i==1) ? 1 : 2^(i-1)+1

type Grid
	d::Int64
	q::Int64
	n::Int64
	grid::Array{Float64,2}
	bounds::Array{Float64,2}
	active::Vector{Bool}
	level::Vector{Int}
end

Grid(D::Int,Q::Int,bounds = [zeros(1,D);ones(1,D)]) =	Grid(ones(Int,D)*Q,bounds)
Grid(D::Int,Q::Vector{Int},bounds = [zeros(1,D);ones(1,D)]) = Grid(Q,bounds)

function Grid(Q::Vector{Int},bounds::Array{Float64,2})
	if all(bounds.==0.0)
		bounds = [zeros(1,d);ones(1,d)]
	end
	X = ndgrid(ntuple(i->linspace(bounds[1,i],bounds[2,i],Mi(Q[i]+1)),length(Q))...)
	x = hcat([x[:] for x in X]...)
	D = length(Q)
	N = prod(map(Mi,Q+1))
	return Grid(D,maximum(Q),N,x,bounds,ones(Bool,size(x,1)),ones(Int,size(x,1)))
end



Base.values(G::Grid)= G.grid
Base.values(G::Grid,i::Int)= G.grid[:,i]

function interp(x1::Array{Float64},G::Grid,A::Array{Float64,1})
	@assert length(A) == G.n
	N = ntuple(i->length(unique(G.grid[:,i])),G.d)
	IF = ModGrid.Grid.InterpGrid(reshape(A,N),ModGrid.Grid.BCnearest,ModGrid.Grid.InterpQuadratic)
	x = nXtoX(x1,G.bounds,vcat(ones(1,G.d),hcat(N...)))
	return  Float64[getindex(IF,ntuple(j->x[i,j],G.d)...) for i = 1:size(x,1)]
end


end

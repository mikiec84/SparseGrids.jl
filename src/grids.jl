abstract BasisFunction
abstract Linear  <: BasisFunction
abstract Quadratic <: BasisFunction
abstract GaussianRadialBF <: BasisFunction

immutable NodeLocation
    l::UInt8
    j::UInt16
    dj::UInt16
end
NodeLocation(x::Float64) = NodeLocation(position(x)...)

type AdaptiveGrid
    nodeinfo
    active::BitArray{1}
    overlap::Vector{Vector{Bool}}
end

type NGrid{D,B}
    L::Vector{Int}
    bounds::Array{Float64,2}
    grid::Array{Float64,2}
	covers::Array{UInt16,2}
    covers_dM::Array{UInt16,2}
	covers_loc::Vector{UInt16}

    adapt::AdaptiveGrid

    IDs::Vector{Vector{Int}}
    Bs::Vector{Vector{Float64}}
end

function Base.kron(X::Vector{Vector{Float64}})
    D = length(X)
    lxs = Int[length(x) for x in X]
    G = zeros(prod(lxs),D)
    s=1
    for d=1:D
        snext = s*lxs[d]
        for j = 1:prod(lxs)
            G[j,d] = X[d][div(rem(j-1, snext), s)+1]
        end
        s = snext
    end
    return G
end

Base.kron(x::Vector{Float64}) = x

function TensorGrid(L::Vector{Int})
    G = ndgrid(Vector{Float64}[cc_g(i) for i in L]...)
    G = hcat([vec(g) for g in G]...)
end

function SmolyakSize(L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    m = Int[cc_dM(l) for l = 1:maximum(L)+D]
    S = 0
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                s=m[covering[1]]
                for i = 2:length(covering)
                    s*=m[covering[i]]
                end
                S+=s
            end
        end
    end
    return S
end

function SmolyakGrid(L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    dg = Vector{Float64}[cc_dg(l) for l = 1:maximum(L)+D]
    G = Array(Array{Float64},0)
    index = Array(Array{Int},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(dg[covering]))
                push!(index,repmat(covering',size(G[end],1)))
            end
        end
    end
    G = vcat(G...)::Array{Float64,2}
    index = vcat(index...)::Array{Int,2}
    return G,index
end

function level(x::Float64)
    l=0
    if x==0.5
        l =  1
    elseif x==0.0 || x==1.0
        l = 2
    else
        for l = 3:12
            mod(x,1/2^(l-1))==0.0 && break
        end
    end
    return l
end
level(G::NGrid) = level(G.grid)
level(X::Array{Float64}) = vec(sum(map(level,X),2)-size(X,2))


level_loc(lev::Vector) = [Int32[findfirst(lev.==l) for l = 0:maximum(lev)];Int32(length(lev)+1)]::Vector{Int32}

function nextid(ind)
    C = unique(ind,1)
    Cl = Vector{Int}[(1:size(ind,1))[vec(all(ind.==C[i,:]',2))] for i = 1:size(C,1)]
    itoCi = [findfirst([in(i,c) for c in Cl]) for i = 1:size(ind,1)]
    nextid=Int32[findfirst(itoCi.>itoCi[i]) for i = 1:size(ind,1)]
    nextid[nextid.==0]=(size(ind,1)+1)
    return nextid
end


"""
    NGrid(L,bounds,B)

# Construction
Construct a Smolyak grid.
The vector of integers L determines the maximum level of the grid in each
dimension. The bounds can be optionally changed from the default of [0,1]^D and
basis function B can be either Linear or Quadratic.

# Interpolation
Grid objects are callable taking two arguments. The first is a vector containing
function values at the grid nodes. The second array contains rows of points at
which the interpolant is to be evaluated.
"""
function NGrid{BT<:BasisFunction}(L::Vector{Int},bounds::Array{Float64} = [0,1]*ones(1,length(L));B::Type{BT}=Linear)
    grid,ind=SmolyakGrid(L)
	covers = map(UInt16,unique(ind,1))
    weights= zeros(size(grid,1))
    covers_dM = map(x->cc_dM(Int(x)),covers)
    covers_loc = zeros(UInt32,size(covers,1))
    hind = [hash(ind[i,:]) for i = 1:size(ind,1)]
    hcovers = [hash(covers[i,:]) for i = 1:size(covers,1)]
    for i=1:size(covers,1)
        covers_loc[i] = findfirst(hind,hcovers[i])
    end

    G = NGrid{length(L),B}(L,
                        bounds,
                        grid,
                        covers,
                        covers_dM,
                        covers_loc,
                        AdaptiveGrid([],zeros(Bool,size(grid,1)),[]),
                        Vector{Int}[Int[] for i = 1:size(grid,1)],
                        Vector{Float64}[Float64[] for i = 1:size(grid,1)]
                        )
    G.adapt.active = (level(G).==maximum(level(G)))
    buildW(G,hind,hcovers)
    return G
end

Base.show(io::IO,G::NGrid) = println(typeof(G),": $(size(G.grid,1))pts")
Base.length(G::NGrid) = size(G.grid,1)
Base.size(G::NGrid) = size(G.grid)
Base.values(G::NGrid) = nUtoX(G.grid,G.bounds)

# Base.getindex(G::NGrid,args...) = getindex(G.grid,args...)
Base.getindex(G::NGrid,args...) = getindex(G.grid,args...)

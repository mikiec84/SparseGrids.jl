abstract BasisFunction
type LinearBF  <: BasisFunction end
type QuadraticBF <: BasisFunction end
type GaussianRadialBF <: BasisFunction end


type NGrid{T,B}
    L::Vector{Int}
    bounds::Array{Float64,2}
    level::Vector{Int16}
    level_M::Array{Int16,2}
	level_loc::Vector{Int32}
    index::Array{Int16,2}
    grid::Array{Float64,2}
    nextid::Vector{Int32}
    nextc::Vector{Int32}
	coverings::Array{Int16,2}
    coverings_dM::Array{Int16,2}
	coveringsloc::Tuple{Vector{Int32},Vector{Int32}}
    coveringsid::Vector{Int32}
    IDs::Vector{Vector{Int}}
    Bs::Vector{Vector{Float64}}
    weights::Vector{Float64}
    active::BitArray{1}
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

function TensorGrid(gd::UnivariateGrid,L::Vector{Int})
    G = ndgrid(Vector{Float64}[gd.g(i) for i in L]...)
    G = hcat([vec(g) for g in G]...)
end

function SmolyakSize(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    m = Int[ug.dM(l) for l = 1:maximum(L)+D]
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

function SmolyakGrid(ug::UnivariateGrid,L::Vector{Int},R,mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    g = Vector{Float64}[ug.g(l) for l = 1:maximum(L)+D]
    dw = Vector{Float64}[cc_dsimpsonsw(l) for l = 1:maximum(L)+D]
    G = Array(Array{Float64},0)
    W = Array(Array{Float64},0)
    index = Array(Array{Int},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(g[covering]))
                tw = dw[covering]
                push!(W, vec(prod(kron(Vector{Float64}[tw[i]*R[i] for i = 1:D]),2)))
                push!(index,repmat(covering',size(G[end],1)))
            end
        end
    end
    G = vcat(G...)::Array{Float64,2}
    W = vcat(W...)
    index = vcat(index...)::Array{Int,2}

    uG = unique(G,1)
    uW = zeros(size(uG,1))
    uindex = zeros(Int16,size(uG))
    for i = 1:size(uG,1)
        ids = vec(all(G.==uG[i,:]',2))
        uW[i]=sum(W[ids])
        uindex[i,:]=index[findfirst(ids),:]
    end
    return uG,uindex,uW
end

level(index::Array{Int16,2}) = map(Int16,vec(sum(index,2)-size(index,2)))::Vector{Int16}
level(ug::UnivariateGrid,G::Array{Float64,2}) = level(index(ug,G))
level_loc(lev::Vector{Int16}) = [Int32[findfirst(lev.==l) for l = 0:maximum(lev)];Int32(length(lev)+1)]::Vector{Int32}

function nextid(::Type{CCGrid},ind::Array{Int16})
    C = unique(ind,1)
    Cl = Vector{Int}[(1:size(ind,1))[vec(all(ind.==C[i,:]',2))] for i = 1:size(C,1)]
    itoCi = [findfirst([in(i,c) for c in Cl]) for i = 1:size(ind,1)]
    nextid=Int32[findfirst(itoCi.>itoCi[i]) for i = 1:size(ind,1)]
    nextid[nextid.==0]=(size(ind,1)+1)
    return nextid
end

SparseGrids.nextid(::Type{CC2Grid},ind::Array{Int16}) = SparseGrids.nextid(CCGrid,ind)


function NGrid{T<:GridType,BT<:BasisFunction}(ug::UnivariateGrid{T},L::Vector{Int},bounds::Array{Float64} = [0,1]*ones(1,length(L));B::Type{BT}=LinearBF)
    R = vec(diff(bounds,1))
    grid,ind,weights = SmolyakGrid(ug,L,R)
	coverings = map(UInt16,unique(ind,1))
    nextc = Int[findfirst((sum(coverings,2)-length(L))[i]+1.==(sum(coverings,2)-length(L))) for i = 1:size(coverings,1)]
    nextc[nextc.==0] = size(coverings,1)+10
    cid     = Int32[findfirst(all(ind[i:i,:].== coverings,2)) for i = 1:size(grid,1)]
    c_dM = map(x->ug.dM(Int(x)),coverings)
	coveringsloc = (Int32[findfirst(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)],Int32[findlast(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)])
    lev = level(ind)
	lev_loc = level_loc(lev)
	lev_M = map(i->Int16(ug.M(Int(i))),ind)::Array{Int16,2}
	active = !(lev.<maximum(L))
	nid = nextid(T,ind)
    G = NGrid{T,B}(L,bounds,
                        lev,lev_M,lev_loc,ind,
                        grid,
                        nid,nextc,
                        coverings,c_dM,coveringsloc,cid,
                        Vector{Int}[Int[] for i = 1:size(grid,1)],
                        Vector{Float64}[Float64[] for i = 1:size(grid,1)],
                        weights,active)
    buildW(G)
    return G
end

Base.show(io::IO,G::NGrid) = println(typeof(G),"{$(length(G.L))-D,$(size(G.grid,1))n}")
Base.length(G::NGrid) = size(G.grid,1)
Base.size(G::NGrid) = size(G.grid)
Base.values(G::NGrid) = nUtoX(G.grid,G.bounds)

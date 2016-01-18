abstract BasisFunction
type LinearBF  <: BasisFunction end
type QuadraticBF <: BasisFunction end
type GaussianRadialBF <: BasisFunction end


type NGrid{T,B}
    L::Vector{Int}
    bounds::Array{Float64,2}
    level::Vector{Int16}
    index::Array{Int16,2}
    grid::Array{Float64,2}
    level_M::Array{Int16,2}
    active::BitArray{1}
    nextid::Vector{Int32}
    hashG::Vector{Int}
	coverings::Array{Int16,2}
	coveringsloc::Tuple{Vector{Int32},Vector{Int32}}
	level_loc::Vector{Int32}
    weights::Vector{Float64}
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

function SmolyakGrid(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
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
                push!(W, vec(prod(kron(dw[covering]),2)))
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

function nextid(::Type{MaxGrid},ind::Array{Int16})
	C = unique(ind,1)
	Cl = Vector{Int}[(1:size(ind,1))[vec(all(ind.==C[i,:]',2))] for i = 1:size(C,1)]
	itoCi = [findfirst([in(i,c) for c in Cl]) for i = 1:size(ind,1)]
    nextid = zeros(Int32,size(ind,1))
    for i = 1:size(ind,1)
        nextid[i] = in(1,ind[i,:]) ? i+1 : findfirst(itoCi.>itoCi[i])
    end
	nextid[nextid.==0]=(size(ind,1)+1)
	return nextid
end

function NGrid{T<:GridType,BT<:BasisFunction}(ug::UnivariateGrid{T},L::Vector{Int},bounds::Array{Float64} = [0,1]*ones(1,length(L));B::Type{BT}=LinearBF)
    grid,ind,weights = SmolyakGrid(ug,L)
    hashG=zeros(Int,size(grid,1))
    for i ∈ 1:size(grid,1)
        for d = 1:size(grid,2)
            hashG[i] += ug.itoi(Int(ind[i,d]),clamp(round(Int16,grid[i,d]*(ug.dM(Int(ind[i,d])))+1/2),1,ug.dM(Int(ind[i,d]))),Int(maximum(ind)))
            hashG[i] *= 17
        end
    end
	coverings = map(UInt16,unique(ind,1))
	coveringsloc = (Int32[findfirst(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)],Int32[findlast(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)])
    for i = 1:length(coveringsloc[1])
        id = coveringsloc[1][i]:coveringsloc[2][i]
        id1 = id[sortperm(hashG[id])]
        grid[id,:] = grid[id1,:]
        ind[id,:] = ind[id1,:]
        hashG[id] = hashG[id1]
        weights[id] = weights[id1]
    end
    lev = level(ind)
	lev_loc = level_loc(lev)
	lev_M = map(i->Int16(ug.M(Int(i))),ind)::Array{Int16,2}
	active = !(lev.<maximum(L))
	nid = nextid(T,ind)
    return  NGrid{T,B}(L,bounds,lev,ind,grid,lev_M,active,nid,hashG,coverings,coveringsloc,lev_loc,weights)
end

Base.show(io::IO,G::NGrid) = println(typeof(G),"{$(length(G.L))-D,$(size(G.grid,1))n}")
Base.length(G::NGrid) = size(G.grid,1)
Base.size(G::NGrid) = size(G.grid)
Base.values(G::NGrid) = nUtoX(G.grid,G.bounds)

abstract GridType
type CCGrid <:GridType end
type MaxGrid <:GridType end

abstract BasisFunction
type LinearBF  <: BasisFunction end
type QuadraticBF <: BasisFunction end

type UnivariateGrid{T<:GridType}
    M::Function
    iM::Function
    dM::Function
    g::Function
    dg::Function
	itoj::Function
end

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
end

const CC = UnivariateGrid{CCGrid}(cc_M,cc_iM,cc_dM,cc_g,cc_dg,cc_itoj)
const Max = UnivariateGrid{MaxGrid}(m_M,m_iM,m_dM,m_g,m_dg,m_itoj)


function Base.kron(ug::UnivariateGrid,L::Vector{Int})
    D = length(L)
    gs = map(ug.dg,L)
    lgs = Int[length(g) for g in gs]
    G = zeros(prod(lgs),D)
    s=1
    for d=1:D
        snext = s*lgs[d]
        for j = 1:prod(lgs)
            G[j,d] = gs[d][div(rem(j-1, snext), s)+1]
        end
        s = snext
    end
    return G
end

Base.kron(ug::UnivariateGrid,L::Int) = kron(ug,[L])

function TensorGrid(gd::UnivariateGrid,L::Vector{Int})
    G = ndgrid(Vector{Float64}[gd.g(i) for i in L]...)
    G = hcat([vec(g) for g in G]...)
end

function SmolyakGrid(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    G = Array(Array{Float64},0)
    index = Array(Array{Int},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(ug,covering))
                push!(index,repmat(covering',size(G[end],1)))
            end
        end
    end
    G = vcat(G...)::Array{Float64,2}
    index = vcat(index...)::Array{Int,2}
    index = map(Int16,index)
    return G::Array{Float64,2},index::Array{Int16,2}
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
    grid,ind = SmolyakGrid(ug,L)
	coverings = map(UInt16,unique(ind,1))
	coveringsloc = (Int32[findfirst(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)],Int32[findlast(all(coverings[i:i,:].==ind,2)) for i = 1:size(coverings,1)])
    lev = level(ind)
	lev_loc = level_loc(lev)
	lev_M = map(i->Int16(ug.M(Int(i))),ind)::Array{Int16,2}
	active = !(lev.<maximum(L))
	nid = nextid(T,ind)
	Gji = zeros(Int,size(grid))
	for d = 1:size(grid,2)
		for i ∈ 1:size(grid,1)
		   Gji[i,d] = ug.itoj(Int(ind[i,d]),clamp(round(Int16,grid[i,d]*(ug.dM(Int(ind[i,d])))+1/2),1,ug.dM(Int(ind[i,d]))),Int(maximum(ind)))
		end
	end
	hashG = vcat(Int[hsh(Gji[i,:]) for i = 1:size(grid,1)]...)::Vector{Int}
    return  NGrid{T,B}(L,bounds,lev,ind,grid,lev_M,active,nid,hashG,coverings,coveringsloc,lev_loc)
end

Base.show(io::IO,G::NGrid) = println(typeof(G),"{$(length(G.L))-D,$(size(G.grid,1))n}")
Base.length(G::NGrid) = size(G.grid,1)
Base.size(G::NGrid) = size(G.grid)
Base.values(G::NGrid) = nUtoX(G.grid,G.bounds)

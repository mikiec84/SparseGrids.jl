function shrink!{T<:GridType,BT}(G::NGrid{T,BT},id::BitArray{1})
    G.grid= G.grid[id,:]
    G.index= G.index[id,:]
    G.level= G.level[id]
    G.level_M= G.level_M[id,:]
    G.active= G.active[id]
    G.hashG= G.hashG[id]
    G.coverings = map(UInt16,unique(G.index,1))
	G.coveringsloc = (Int32[findfirst(all(G.coverings[i:i,:].==G.index,2)) for i = 1:size(G.coverings,1)],Int32[findlast(all(G.coverings[i:i,:].==G.index,2)) for i = 1:size(G.coverings,1)])
    G.level_loc = level_loc(G.level)
    G.nextid = nextid(T,G.index)
    return Void
end




function grow!{T<:GridType,BT<:BasisFunction}(G::NGrid{T,BT},id::Int,bounds::Vector{Int}=12*ones(Int,length(G.L)))
    !G.active[id] ? (return nothing) : nothing
    G.active[id] = false

    ug = T==CCGrid ? CC : Max
    x = G.grid[id,:]
    n = length(G.L)*2
    for d = 1:length(G.L)
        n += (x[d]==0.0 || x[d] == 1.0) ? -1 : 0
    end

    X,ind = SmolyakGrid(ug,bounds,G.level[id]+1:G.level[id]+1)

    id1 = BitArray(size(X,1))
    for i = 1:size(X,1)
        id1[i] = !any(all(G.grid.==X[i:i,:],2))
    end
    all(!id1) ? (return nothing) : nothing
    X = X[id1,:]
    ind = ind[id1,:]

    id1 = sortperm(Float64[norm(X[i,:]-vec(x)) for i = 1:size(X,1)])[1:min(n,size(X,1))]
    X = X[id1,:]
    ind = map(Int16,ind[id1,:])
    add!(G,X,ind)
    return nothing
end

function add!{T<:GridType,BT<:BasisFunction}(G::NGrid{T,BT},X::Array{Float64,2},ind::Array{Int16,2})
    ug = T==CCGrid ? CC : Max
    @assert size(X,2) == size(G.grid,2)

    G.grid=[G.grid;X]
    G.index = [G.index;ind]
    G.level = level(G.index)
    G.level_M = [G.level_M;map(i->Int16(ug.M(Int(i))),ind)]
    G.L = vec(maximum(G.index,1)-1)
    Gji = zeros(Int,size(X))
    for d = 1:size(X,2)
		for i ∈ 1:size(X,1)
		   Gji[i,d] = ug.itoj(Int(ind[i,d]),clamp(round(Int16,X[i,d]*(ug.dM(Int(ind[i,d])))+1/2),1,ug.dM(Int(ind[i,d]))),Int(maximum(ind)))
		end
	end
	hashG = vcat(Int[hsh(Gji[i,:]) for i = 1:size(X,1)]...)::Vector{Int}

    G.hashG = [G.hashG;hashG]
    G.active = [G.active;!BitArray(size(X,1))]

    id = sortperm(G.level)
    # t1 = sortrows([G.level G.index G.grid])
    id = sortperm(nextid(T,G.index))

    G.level = G.level[id]
    G.index = G.index[id,:]
    G.grid = G.grid[id,:]
    G.level_M = G.level_M[id,:]
    G.hashG = G.hashG[id]
    G.nextid= nextid(T,G.index)


    G.coverings =  map(UInt16,unique(G.index,1))
    G.coveringsloc = (Int32[findfirst(all(G.coverings[i:i,:].==G.index,2)) for i = 1:size(G.coverings,1)],Int32[findlast(all(G.coverings[i:i,:].==G.index,2)) for i = 1:size(G.coverings,1)])
    G.level_loc = level_loc(G.level)

    return nothing
end

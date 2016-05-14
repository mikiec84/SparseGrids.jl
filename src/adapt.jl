function addcover!{D,B}(G::NGrid{D,B},L::Vector{Int})
    any(all(G.covers.==L',2)) && return
    bf = B==Linear ? cc_bf_l : cc_bf_q
    x   = kron([SparseGrids.cc_dg(i) for i in L])
    ind = repmat(map(Int16,L'),size(x,1))
    level_M = map(i->Int16(SparseGrids.cc_M(level(i))),G.grid)

    for id = 1:size(x,1)
        push!(G.Bs,[])
        push!(G.IDs,Int[])
        for i = 1:length(G)
            b = 1.0
            for d = 1:length(G.L)
                b*=bf(x[id,d],G.grid[i,d],level_M[i,d])
            end
            b>0 && push!(G.Bs[end],-b)
            b>0 && push!(G.IDs[end],i)
        end
        push!(G.Bs[end],1.0)
        push!(G.IDs[end],length(G)+id)
    end
    G.grid = [G.grid; x]
    ind = map(level,G.grid)

    G.covers = [G.covers;L']
    G.covers_loc = Int32[findfirst(all(G.covers[i:i,:].==ind,2)) for i = 1:size(G.covers,1)]
    G.covers_dM = map(x->SparseGrids.cc_dM(Int(x)),G.covers)
    G.adapt.active = [G.adapt.active;zeros(Bool,size(x,1))]
    G.L = vec(maximum(G.covers,1))-1
    return
end

function Base.sort(G::NGrid)
    lev = level(G)
    try
        for i = 2:length(G)
            @assert lev[i-1]≤lev[i]
        end
    catch
        id = zeros(Int,0)
        for l = 0:maximum(lev)
            id = [id;find(lev.==l)]
        end
        G.grid = G.grid[id,:]
        ind = map(level,G.grid)
        G.covers = map(UInt16,unique(ind,1))
        G.covers_dM = map(x->cc_dM(Int(x)),G.covers)
        G.covers_loc = Int32[findfirst(all(G.covers[i:i,:].==ind,2)) for i = 1:size(G.covers,1)]
        SparseGrids.buildW(G)
    end
end

function iscoverfull(G::NGrid,L::Vector{Int})
    @assert length(L)==length(G.L)
    n = 0
    D= length(L)
    for i = 1:length(G)
        b = true
        for j = 1:D
            b = b && level(G.grid[i,j])==L[j]
        end
        b && (n+=1)
    end
    return n == prod(map(cc_dM,L))
end

looseid(G::NGrid) = G.covers_loc[end]+prod(G.covers_dM[end,:]):size(G.grid,1)

function grow!{D,BF}(G::NGrid{D,BF},id::Int,bounds::Vector{Int}=12*ones(Int,length(G.L)))
    !G.adapt.active[id] && (return nothing)
    G.adapt.active[id] = false
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)

    x = G.grid[id,:]
    n = D*2
    for d = 1:D
        n += (x[d]==0.0 || x[d] == 1.0) ? -1 : 0
    end
    Lid = sum(map(level,G.grid[id,:]))-D
    X,ind = SmolyakGrid(bounds,Lid+1:Lid+1)

    id1 = BitArray(size(X,1))
    for i = 1:size(X,1)
        id1[i] = !any(all(G[:,:].==X[i:i,:],2))
    end
    all(!id1) && (return nothing)
    X = X[id1,:]
    id1 = sortperm(Float64[norm(X[i,:]-vec(x)) for i = 1:size(X,1)])[1:min(n,size(X,1))]
    X = X[id1,:]

    Bs = Vector{Float64}[]
    IDs = Vector{Int}[]

    level_M = map(i->Int16(SparseGrids.cc_M(level(i))),G.grid)
    for j = 1:size(X,1)
        push!(Bs,[])
        push!(IDs,Int[])
        for i = 1:length(G)
            b = 1.0
            for d = 1:length(G.L)
                b*=bf(X[j,d],G.grid[i,d],level_M[i,d])
            end
            b>0 && push!(Bs[end],-b)
            b>0 && push!(IDs[end],i)
        end
        push!(Bs[end],1.0)
        push!(IDs[end],length(G)+j)
    end

    insertloc = size(G.grid,1)
    resort = length(looseid(G))>0 && Lid+1<maximum(level(G.grid[looseid(G),:]))
    if resort
        insertloc = findfirst(level(G.grid[looseid(G),:]).>Lid+1) + first(looseid(G)) - 2
    end

    G.grid  = [G.grid[1:insertloc,:];X;G.grid[insertloc+1:end,:]]
    G.Bs    = [G.Bs[1:insertloc];Bs;G.Bs[insertloc+1:end]]
    G.IDs    = [G.IDs[1:insertloc];IDs;G.IDs[insertloc+1:end]]

    G.adapt.active     = [G.adapt.active[1:insertloc];ones(Bool,size(X,1));G.adapt.active[insertloc+1:end]]

    if resort
        rebuildW(G,insertloc+size(X,1):length(G))
    end
    return
end

function rebuildW{D,BF}(G::NGrid{D,BF},idr)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    level_M = map(i->Int16(SparseGrids.cc_M(level(i))),G.grid)
    for j = idr
        G.Bs[j] = Float64[]
        G.IDs[j] = Int[]
        for i = 1:j-1
            b = 1.0
            for d = 1:D
                b*=bf(G.grid[j,d],G.grid[i,d],level_M[i,d])
            end
            b>0 && push!(G.Bs[j],-b)
            b>0 && push!(G.IDs[j],i)
        end
        push!(G.Bs[j],1.0)
        push!(G.IDs[j],length(G)+j)
    end
end


function getparents{D,BF}(G::NGrid{D,BF})
    nc = size(G.covers,1)
    topc = [findlast((G.covers[:,d].==maximum(G.covers[:,d]))) for d=1:D]
    topcp = UnitRange{Int}[G.covers_loc[c]:(c == nc ? length(G) : G.covers_loc[c+1]-1) for c in topc]
    # bT = Array(Bool,length(G),length(G))
    bT = BitArray(length(G),length(G))
    for i = 1:length(G)
        for j = 1:i
            x = intersect(box(G,i),box(G,j))
            bT[i,j] = x
            i!=j && (bT[j,i] = x)
        end
    end
    T = Vector{UInt32}[map(UInt32,find(bT[i,:])) for i = 1:length(G)]
    Ts = vcat([[length(T[i[1]]) length(i)] for i in topcp]...)
    Tp = vcat(vcat([T[i] for i in topcp]...)...)
    bT,T,Ts,Tp
end

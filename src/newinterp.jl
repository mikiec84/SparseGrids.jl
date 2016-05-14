type box{D}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function box{D,BF}(G::NGrid{D,BF},i::Int)
    lower = Array(Float64,D)
    upper = Array(Float64,D)
    for d = 1:D
        p = SparseGrids.position(G.grid[i,d])
        lower[d] = clamp(G.grid[i,d]-1/2^(p[1]-1),0,1)+eps(Float64)
        upper[d] = clamp(G.grid[i,d]+1/2^(p[1]-1),0,1)-eps(Float64)
    end
    return box{D}(lower,upper)
end

function Base.intersect{D}(a::box{D},b::box{D})
    nool = a.lower[1] > b.upper[1] ||
           b.lower[1] > a.upper[1]
    for d = 2:D
        nool = nool ||
           a.lower[d] > b.upper[d] ||
           b.lower[d] > a.upper[d]
       end
    return !nool
end

function getparents(G::NGrid)
    # T = zeros(Bool,length(G),length(G))
    T = Array(Bool,length(G),length(G))
    for i = 1:length(G)
        for j = 1:i
            x = intersect(box(G,i),box(G,j))
            T[i,j] = x
            i!=j && (T[j,i] = x)
        end
    end
    return Vector{Int}[find(T[i,:]) for i = 1:length(G)]
end


function Base.intersect(A::Vector{Int}, B::Vector{Int},n::Int)
    ret = zeros(Int,n)
    na,nb = length(A),length(B)
    cnt = 1
    for ia = 1:na
        a = A[ia]
        ib = ia
        b = B[ib]
        if a<=b
            while a!=b && ib<nb
                ib-=1
                b = B[ib]
            end
            a==B[ib] && (ret[cnt] = b;cnt+=1)

        else
            while a!=b && ib<nb
                ib+=1
                b = B[ib]
            end
            a==B[ib] && (ret[cnt] = b;cnt+=1)
        end
        cnt>n && break
    end
    ret
end


function interp1{D,BF}(G::NGrid{D,BF},T::Vector{Vector{Int}},A::Vector{Float64},x::Array{Float64})
    bf = (BF == Linear ? cc_bf_l : cc_bf_q)
    w = getW(G,A)
    y = zeros(size(x,1))
    D = length(G.L)
    mL      = maximum(G.L)
    lev = map(level,G[:,:])
    nc = size(G.covers,1)
    J         = ones(Int,mL+1,D)
    B         = zeros(mL+1,D)

    icmax = (sum(G.covers,2)-D).==maximum(sum(G.covers,2)-D)
    inmax = level(G).==maximum(level(G))
    lmax = [maximum(G.covers[:,d])::Int16 for d = 1:length(G.L)]

    topn = [find((map(level,G[:,d]).==lmax[d]).*inmax) for d=1:D]
    topc = [findfirst((G.covers[:,d].==lmax[d])) for d=1:D]
    for i = 1:size(x,1)
        for d = 1:D
            for l = 1:mL+1
                j     = clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
                B[l,d]     = bf(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
                J[l,d]  = j
            end
        end
        id = []
        for ii in topc
            b  = B[G.covers[ii,D],D]*B[G.covers[ii,1],1]
            id1 = J[G.covers[ii,D],D]-1
            for d = D-1:-1:2
                b*=B[G.covers[ii,d],d]
                id1 = id1*G.covers_dM[ii,d]+(J[G.covers[ii,d],d]-1)
            end
            push!(id,(J[G.covers[ii,1],1]-1)+G.covers_dM[ii,1]*id1+1+G.covers_loc[ii]-1)
        end
        T0 = intersect(T[id[1]],T[id[2]])
        for ii in 1:nc
            b = 1.0
            for d = 1:D
                b*=B[G.covers[ii,d],d]
            end
            y[i]+=b*w[T0[ii]]
        end
    end
    y
end

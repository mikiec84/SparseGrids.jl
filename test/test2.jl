import SparseGrids: getparents


@inline cc_M{T}(i::T)           = (i==1) ? T(1) : T(2)^(i-T(1))+T(1)
@inline function cc_bf_l(x::Float64,xij::Float64,mi)
    if (mi==1)
        return 1.0
    end
    dx = (1.0-(mi-1.0)*abs(x-xij))
    return (dx>0.0) ? dx : 0.0
end
@inline function cc_bf_q(x::Float64,xij::Float64,mi)
    if (mi==1)
        return 1.0
    end
    dx = 1.0-((mi-1.0)*(x-xij))^2
    return (dx>0.0) ?  dx : 0.0
end


function ainterp{D,BF}(G::NGrid{D,BF},A,x)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    nG = length(G)
    y = G(A,x)
    w = getW(G,A)
    for i = 1:length(y)
        for j = (G.covers_loc[end]+prod(G.covers_dM[end,:])):nG
            b = 1.0
            for d = 1:D
                b *= bf(x[i,d],G.grid[j,d],cc_M(level(G.grid[j,d])))
                b==0 && break
            end
            b>0 && (y[i] += w[j]*b)
        end
    end
    return y
end


G.covers_dM

function interp1{D,BF}(G::NGrid{D,BF},A,x)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    nG = length(G)
    w = getW(G,A)
    y = zeros(size(x,1))
    for i = 1:length(y)
        for j = 1:nG
            b = 1.0
            for d = 1:D
                b *= bf(x[i,d],G.grid[j,d],cc_M(level(G.grid[j,d])))
                b==0 && break
            end
            b>0 && (y[i] += w[j]*b)
        end
    end
    return y
end

function interp11{D,BF}(G::NGrid{D,BF},A,x)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    nG = length(G)
    w = getW(G,A)
    y = zeros(size(x,1))
    pos = map(SparseGrids.position,G.grid)
    mL      = 8
    B         = zeros(mL+1,D)

    for i = 1:length(y)
        for d = 1:D
            for l = 1:mL+1
                j     = clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
                B[l,d]     = bf(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
            end
        end
        for j = 1:nG
            b = 1.0
            for d = 1:D
                b *= B[pos[j,d][2]]
                b==0 && break
            end
            b>0 && (y[i] += w[j]*b)
        end
    end
    return y
end



function interp2{D,BF}(G::NGrid{D,BF},A,x,baT)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    nG = length(G)
    w = getW(G,A)
    y = zeros(size(x,1))
    id = BitArray(nG)
    for i = 1:length(y)
        for j = 1:nG
            id[j] = true
        end
        j = 1
        while j<=nG
            if id[j]
                b = 1.0
                for d = 1:D
                    b *= bf(x[i,d],G.grid[j,d],cc_M(level(G.grid[j,d])))
                    b==0 && break
                end
                if b>0
                    (y[i] += w[j]*b)
                    id = id & sub(baT,:,j)
                end
            end
            j+=1
        end
    end
    return y
end


function interp3{D,BF}(G::NGrid{D,BF},A,x,baT,jlim=D)
    bf = (BF==Linear ? cc_bf_l : cc_bf_q)
    nG = length(G)
    w = getW(G,A)
    y = zeros(size(x,1))
    id = BitArray(nG)
    for i = 1:length(y)
        for j = 1:nG
            id[j] = true
        end
        j = 1
        jit=0
        while j<=nG
            if id[j]
                b = 1.0
                for d = 1:D
                    b *= bf(x[i,d],G.grid[j,d],cc_M(level(G.grid[j,d])))
                    b==0 && break
                end
                if b>0
                    (y[i] += w[j]*b)

                    (jit<jlim) && (id = id & sub(baT,:,j))
                    jit+=1
                end
            end
            j+=1
        end
    end
    return y
end



G = NGrid([2,2,2,2,2])
f(x) = vec(prod(x[:,:].^5,2))
x = SparseGrids.TensorGrid([4,4,4,4,4])

for i = 1:30
    id = indmax(abs(getW(G,f(G))).*G.adapt.active)
    grow!(G,id)
end


(baT,T,Ts,Tp)=getparents(G);

ay = SparseGrids.ainterp(G,f(G),X)
y1 = interp1(G,f(G),X)
y11 = interp11(G,f(G),X)
y2 = interp2(G,f(G),X,baT)
y3 = interp3(G,f(G),X,baT,5)
@assert maximum(abs(ay-y1))<1e-15
@assert maximum(abs(ay-y2))<1e-15
@assert maximum(abs(ay-y3))<1e-15



# ab = @benchmark ainterp(G,f(G),X)
# a1 = @benchmark interp1(G,f(G),X)
# a2 = @benchmark interp2(G,f(G),X,baT)
# a3 = @benchmark interp3(G,f(G),X,baT)


nrep = 1
@time for i = 1:nrep SparseGrids.ainterp(G,f(G),X) end
@time for i = 1:nrep interp1(G,f(G),X) end
@time for i = 1:nrep interp11(G,f(G),X) end
@time for i = 1:nrep interp2(G,f(G),X,baT) end
@time for i = 1:nrep interp3(G,f(G),X,baT) end
@time for i = 1:nrep interp3(G,f(G),X,baT,1) end

@pview interp11(G,f(G),X)

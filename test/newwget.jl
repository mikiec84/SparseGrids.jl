import SparseGrids:cc_bf_l,cc_dM,cc_M,cc_dg


function buildW(G::NGrid)
    nc 		= size(G.coverings,1)
    mL      = maximum(G.L)
    D       = length(G.L)

    J 		= ones(Int,mL+1,D)
    B 		= zeros(mL+1,D)

    IDs = Vector{Int}[Int[] for i = 1:length(G)]
    Bs  = Vector{Float64}[Float64[] for i = 1:length(G)]

    for i = 1:length(G)
        for d = 1:D
            for l = 1:G.level[i]+1
                j 	= clamp(round(Int,G.grid[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
                B[l,d] 	= cc_bf_l(G.grid[i,d],cc_dg(l,j),Int16(cc_M(l)))
                J[l,d]  = j
            end
        end
        for ii = 1:G.coveringsid[i]
            b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
            id1 = J[G.coverings[ii,D],D]-1
            for d = D-1:-1:2
                b*=B[G.coverings[ii,d],d]
                id1 = id1*G.coverings_dM[ii,d]+(J[G.coverings[ii,d],d]-1)
            end
            id1=(J[G.coverings[ii,1],1]-1)+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
            if b>0
                push!(IDs[i],id1)
                push!(Bs[i],id1==i ? b : -b)
            end
        end
    end
    G.IDs=IDs
    G.Bs = Bs
    return
end


function getW1(G,A,IDs,Bs)
    w = zeros(length(G))
    for i = 1:length(G)
        id = IDs[i]
        B = Bs[i]
        for l ∈ 1:length(B)
            w[i]+= id[l]==i ? B[l]*A[id[l]] : B[l]*w[id[l]]
        end
    end
    return w
end

function getW2(G,A,IDs,Bs)
    w = zeros(length(G))
    for i = 1:length(G)
        for l ∈ 1:length(Bs[i])
            w[i]+= IDs[i][l]==i ? Bs[i][l]*A[IDs[i][l]] : Bs[i][l]*w[IDs[i][l]]
        end
    end
    return w
end

function getW3(G,A)
    w = copy(A)
    for i = 1:length(G)
        for l ∈ 1:length(G.Bs[i])-1
            w[i]+=  G.Bs[i][l]*w[G.IDs[i][l]]
        end
    end
    return w
end



G = NGrid(CC,[5,5,5,5,5])
A = rand(length(G))
IDs,Bs=buildW(G)
w1 = getW1(G,A,IDs,Bs)
w2 = getW2(G,A,IDs,Bs)
w3 = getW3(G,A,IDs,Bs)
W = SparseGrids.c_getW(G,A)
@assert all(abs(w1-W).<1e-12)
@assert all(abs(w2-W).<1e-12)
@assert all(abs(w3-W).<1e-12)





begin
@time for i = 1:10;SparseGrids.jl_getW(G,A);end;
@time for i = 1:100;getW1(G,A,IDs,Bs);end;
@time for i = 1:100;getW2(G,A,IDs,Bs);end;
@time for i = 1:100;getW3(G,A,IDs,Bs);end;
end

function jl_interpbig1(xi::Array{Float64},G::NGrid,A::Vector{Float64})
    w 		= getW(G,A)
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(x,1))
    nc 		= size(G.coverings,1)
    mL      = maximum(G.L)
    D       = length(G.L)

    for i = 1:size(x,1)
        J 		= ones(Int,mL+1,D)
        B 		= zeros(mL+1,D)
        id  	= zeros(Int,D)
        for d = 1:D
            for l = 1:mL+1
                j 	= clamp(round(Int,x[i,d]*(SparseGrids.cc_dM(l))+1/2),1,SparseGrids.cc_dM(l))
                B[l,d] 	= SparseGrids.cc_bf_l(x[i,d],SparseGrids.cc_dg(l,j),Int16(SparseGrids.cc_M(l)))
                J[l,d]  = j
            end
        end

        for ii = 1:nc
            b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
            id1 = J[G.coverings[ii,D],D]-1
            for d = D-1:-1:2
                b*=B[G.coverings[ii,d],d]
                id1 *= G.coverings_dM[ii,d]
                id1 += (J[G.coverings[ii,d],d]-1)
            end
            id1=(J[G.coverings[ii,1],1]-1)+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
            y[i]+=b*w[id1]
        end
    end
    y
end

function jl_interpbig2(xi::Array{Float64},G::NGrid,A::Vector{Float64})
    w 		= getW(G,A)
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(x,1))
    nc 		= size(G.coverings,1)
    mL      = maximum(G.L)
    D       = length(G.L)

    for i = 1:size(x,1)
        J 		= ones(Int,mL+1,D)
        B 		= zeros(mL+1,D)
        id  	= zeros(Int,D)
        for d = 1:D
            for l = 1:mL+1
                j 	= clamp(round(Int,x[i,d]*(SparseGrids.cc_dM(l))+1/2),1,SparseGrids.cc_dM(l))
                B[l,d] 	= SparseGrids.cc_bf_l(x[i,d],SparseGrids.cc_dg(l,j),Int16(SparseGrids.cc_M(l)))
                J[l,d]  = j-1
            end
        end

        for ii = 1:nc
            b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
            id1 = J[G.coverings[ii,D],D]
            for d = D-1:-1:2
                b*=B[G.coverings[ii,d],d]
                id1 *= G.coverings_dM[ii,d]
                id1 += J[G.coverings[ii,d],d]
            end
            id1=J[G.coverings[ii,1],1]+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
            y[i]+=b*w[id1]
        end
    end
    y
end

G = NGrid(CC,[5,5,5,2,2,2])
A = rand(length(G))


#
# b=@benchmark SparseGrids.jl_interpbig(G.grid,G,A);
# b1=@benchmark jl_interpbig1(G.grid,G,A);
# b2=@benchmark jl_interpbig2(G.grid,G,A);
b = [(tic();SparseGrids.jl_interpbig(G.grid,G,A);toq()) for i = 1:1000]
b1 = [(tic();jl_interpbig1(G.grid,G,A);toq()) for i = 1:1000]
b2 = [(tic();jl_interpbig2(G.grid,G,A);toq()) for i = 1:1000]

using HypothesisTests
# ApproximateTwoSampleKSTest(b1.samples.elapsed_times,b2.samples.elapsed_times)
ApproximateTwoSampleKSTest(b,b1)
ApproximateTwoSampleKSTest(b1,b2)
ApproximateTwoSampleKSTest(b,b2)

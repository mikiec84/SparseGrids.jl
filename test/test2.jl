function c_interpbigbase(xi::Array{Float64},G::NGrid{CCGrid,LinearBF},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z19interp_cc_l_bigbaselllllPdS_S_PsS0_PiS1_,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

function c_interpbigbase(xi::Array{Float64},G::NGrid{CCGrid,QuadraticBF},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z19interp_cc_q_bigbaselllllPdS_S_PsS0_PiS1_,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

function c_interpbig2(xi::Array{Float64},G::NGrid{CCGrid,QuadraticBF},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z16interp_cc_q_big2lllllPdS_S_PsS0_PiS1_,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

G = NGrid(CC,[4,4,4,1,1])
Gq = NGrid(CC,G.L,B=QuadraticBF)
A = rand(length(G))


function runtest(G,A)
    x = rand(5000,length(G.L))
    b,b1,b2=0.0,0.0,0.0
    for i = 1:1000
        tic(); SparseGrids.c_interpbig(x,G,A); b1+=toq()
        tic(); SparseGrids.c_interp(x,G,A); b+=toq()
    end
    println(b/b1)
end
runtest(Gq,A)



@assert maximum(abs(c_interpbigbase(G.grid,G,A)-A))<1e-15
@assert maximum(abs(c_interpbigbase(G.grid,Gq,A)-A))<1e-15
@assert maximum(abs(c_interpbig2(G.grid,Gq,A)-A))<1e-15

function runtest(G,A)
    x = rand(5000,length(G.L))
    b,b1,b2=0.0,0.0,0.0
    for i = 1:1000
        tic(); c_interpbigbase(x,G,A);   b+=toq()
        tic(); SparseGrids.c_interpbig(x,G,A); b1+=toq()
        tic(); c_interpbig2(x,G,A);   b2+=toq()
    end
    println(b/b1,"  ",b/b2)
end
runtest(Gq,A)





#
# println(mean(b)/mean(b1))
# println(mean(b)/mean(b2))
# using HypothesisTests
# # ApproximateTwoSampleKSTest(b1.samples.elapsed_times,b2.samples.elapsed_times)
# ApproximateTwoSampleKSTest(b,b1)
# ApproximateTwoSampleKSTest(b1,b2)
# ApproximateTwoSampleKSTest(b,b2)


nothing

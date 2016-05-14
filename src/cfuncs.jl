const lsparse = Pkg.dir("SparseGrids")*"/deps/libsparse.so"

###############################################################################
#  CC   Linear
###############################################################################


function c_interp{D}(G::NGrid{D,Linear},A::Vector{Float64},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(xi,1))
    w 		= jl_getWbig(G,A)
    ccall((:_Z8interp_lllllPdS_S_PsS0_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32}),
        length(G.L),maximum(G.L),size(x,1),size(G.covers,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc))
    return y
end

function c_interp{D}(G::NGrid{D,Linear},A::Array{Float64,2},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)

    ccall((:_Z8interp_lllllllPdS_S_PsS0_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),size(A,2),
        pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc))
    return y
end

function c_interp_a{D}(G::NGrid{D,Linear},A::Array{Float64,1},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)
    Mgrid = map(Int32,map(SparseGrids.cc_M,map(level,G.grid)))

    ccall((:_Z8interp_llllllPdS_S_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        D,maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),
        pointer(G.grid),pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc),pointer(Mgrid))
    return y
end

function c_interp_a{D}(G::NGrid{D,Linear},A::Array{Float64,2},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)
    Mgrid = map(Int32,map(SparseGrids.cc_M,map(level,G.grid)))

    ccall((:_Z8interp_lllllllPdS_S_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        D,maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),size(A,2),
        pointer(G.grid),pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc),pointer(Mgrid))
    return y
end

###############################################################################
#  CC   Quadratic
###############################################################################


function c_interp{D}(G::NGrid{D,Quadratic},A::Vector{Float64},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(xi,1))
    w 		= jl_getWbig(G,A)
    ccall((:_Z8interp_qllllPdS_S_PsS0_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32}),
        length(G.L),maximum(G.L),size(x,1),size(G.covers,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc))
    return y
end

function c_interp{D}(G::NGrid{D,Quadratic},A::Array{Float64,2},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)

    ccall((:_Z8interp_qllllllPdS_S_PsS0_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),size(A,2),
        pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc))
    return y
end

function c_interp_a{D}(G::NGrid{D,Quadratic},A::Array{Float64,1},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)
    Mgrid = map(Int32,map(SparseGrids.cc_M,map(level,G.grid)))

    ccall((:_Z8interp_qlllllPdS_S_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        D,maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),
        pointer(G.grid),pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc),pointer(Mgrid))
    return y
end


function c_interp_a{D}(G::NGrid{D,Quadratic},A::Array{Float64,2},xi::Array{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)
    Mgrid = map(Int32,map(SparseGrids.cc_M,map(level,G.grid)))

    ccall((:_Z8interp_qllllllPdS_S_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        D,maximum(G.L),size(G.grid,1),size(x,1),size(G.covers,1),size(A,2),
        pointer(G.grid),pointer(x),pointer(w),pointer(y),
        pointer(G.covers),pointer(G.covers_dM),pointer(G.covers_loc),pointer(Mgrid))
    return y
end








function cintersect(A::Vector{UInt32},B::Vector{UInt32})
    C  = Array(UInt32,length(A))
    n = ccall((:_Z9intersectPKjmS0_mPj,SparseGrids.lsparse),Csize_t,(Ptr{UInt32},Csize_t,Ptr{UInt32},Csize_t,Ptr{UInt32}),pointer(A),length(A),pointer(B),length(B),pointer(C))
    C[1:n]
end

function cintersect(A::Vector{UInt32},B::Vector{UInt32},C::Vector{UInt32})
    n = ccall((:_Z9intersectPKjmS0_mPj,SparseGrids.lsparse),Csize_t,(Ptr{UInt32},Csize_t,Ptr{UInt32},Csize_t,Ptr{UInt32}),pointer(A),length(A),pointer(B),length(B),pointer(C))
    return n
end

function cintersect(A::Vector{UInt32},B::Vector{UInt32},C::Vector{UInt32},na::Int,nb::Int)
    n = ccall((:_Z9intersectPKjmS0_mPj,SparseGrids.lsparse),Csize_t,(Ptr{UInt32},Csize_t,Ptr{UInt32},Csize_t,Ptr{UInt32}),pointer(A),na,pointer(B),nb,pointer(C))
    return n
end

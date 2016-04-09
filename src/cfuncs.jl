const lsparse = Pkg.dir("SparseGrids")*"/deps/libsparse.so"

###############################################################################
#  CC   Linear
###############################################################################

function c_getW(G::NGrid{Linear},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((:_Z6w_cc_liiiPdS_PsPiS_S_S_,lsparse),
        Void,
        (Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Int16},Ptr{Float64},
        Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G),size(G.grid,2),maximum(G.L),
        pointer(G.grid),pointer(A),pointer(G.level_M),pointer(G.level_loc),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end

function c_getW(G::NGrid{Linear},A::Array{Float64,2})
    Aold 	= zeros(size(A))
    w 		= zeros(size(A))
    ccall((:_Z10w_cc_l_arrPdiiiiPsPiS_S_S_,lsparse),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),size(A,2),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(w))
    return w
end

function c_interp(xi::Array{Float64},G::NGrid{Linear},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= c_getW(G,A)
    ccall((:_Z11interp_cc_liiiiPdPsS_S_S_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interp(xi::Array{Float64},G::NGrid{Linear},A::Array{Float64,2})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(x,1),size(A,2))
    w 		= c_getW(G,A)
    ccall((:_Z15interp_cc_l_arriiiiiPdPsS_S_S_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(A,2),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{Linear},A::Vector{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(xi,1))
    w 		= jl_getWbig(G,A)
    ccall((:_Z15interp_cc_l_biglllllPdS_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{Linear},A::Array{Float64,2})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)

    ccall((:_Z19interp_cc_l_big_arrllllllPdS_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),size(A,2),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

###############################################################################
#  CC   Quadratic
###############################################################################


function c_getW(G::NGrid{Quadratic},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((:_Z6w_cc_qPdiiiPsPiS_S_S_S_,lsparse),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end

function c_getW(G::NGrid{Quadratic},A::Array{Float64,2})
    Aold 	= zeros(size(A))
    w 		= zeros(size(A))
    ccall((:_Z10w_cc_q_arrPdiiiiPsPiS_S_S_,lsparse),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),size(A,2),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(w))
    return w
end

function c_interp(xi::Array{Float64},G::NGrid{Quadratic},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= c_getW(G,A)
    ccall((:_Z11interp_cc_qiiiiPdPsS_S_S_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interp(xi::Array{Float64},G::NGrid{Quadratic},A::Array{Float64,2})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(x,1),size(A,2))
    w 		= c_getW(G,A)
    ccall((:_Z15interp_cc_q_arriiiiiPdPsS_S_S_Pi,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(A,2),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{Quadratic},A::Vector{Float64})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(xi,1))
    w 		= jl_getWbig(G,A)
    ccall((:_Z15interp_cc_q_biglllllPdS_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{Quadratic},A::Array{Float64,2})
    x 		= clamp!(nXtoU(xi,G.bounds),0.0,1.0)
    y 		= zeros(size(x,1),size(A,2))
    w 		= jl_getWbig(G,A)
    ccall((:_Z19interp_cc_q_big_arrllllllPdS_S_PsS0_PiS1_,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),size(A,2),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]))
    return y
end

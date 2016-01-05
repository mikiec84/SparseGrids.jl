const cc_wget_l = :_Z6w_cc_lPdiiiPsPiS_S_S_S_
const cc_interp_l = :_Z11interp_cc_liiiiPdPsS_S_S_Pi
const cc_interp_l_big = :_Z15interp_cc_l_biglllllPdS_S_PsPiS1_Pl

const cc_wget_q = :_Z6w_cc_qPdiiiPsPiS_S_S_S_
const cc_interp_q = :_Z11interp_cc_qiiiiPdPsS_S_S_Pi
const cc_interp_q_big = :_Z15interp_cc_q_biglllllPdS_S_PsPiS1_Pl

const m_wget_l = :_Z5w_m_lPdiiiPsPiS_S_S_S_
const m_interp_l = :_Z10interp_m_liiiiPdPsS_S_S_Pi
const m_wget_q = :_Z5w_m_qPdiiiPsPiS_S_S_S_
const m_interp_q = :_Z10interp_m_qiiiiPdPsS_S_S_Pi


function c_getW(G::NGrid{CCGrid,LinearBF},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((cc_wget_l,"libsparse.so"),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end

function c_getW(G::NGrid{CCGrid,QuadraticBF},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((cc_wget_q,"libsparse.so"),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end

function c_getW(G::NGrid{MaxGrid,LinearBF},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((m_wget_l,"libsparse.so"),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end

function c_getW(G::NGrid{MaxGrid,QuadraticBF},A::Vector{Float64})
    Aold 	= zeros(length(G))
    dA 		= zeros(length(G))
    w 		= zeros(length(G))
    ccall((m_wget_q,"libsparse.so"),
        Void,
        (Ptr{Float64},Int32,Int32,Int32,
         Ptr{Int16},Ptr{Float64},Ptr{Float64},
         Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        pointer(G.grid),length(G),size(G.grid,2),maximum(G.L),
        pointer(G.level_M),pointer(G.level_loc),pointer(A),
        pointer(Aold),pointer(dA),pointer(w))
    return w
end


function c_interp(xi::Array{Float64},G::NGrid{CCGrid,LinearBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((cc_interp_l,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interp(xi::Array{Float64},G::NGrid{CCGrid,QuadraticBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((cc_interp_q,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interp(xi::Array{Float64},G::NGrid{MaxGrid,LinearBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((m_interp_l,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interp(xi::Array{Float64},G::NGrid{MaxGrid,QuadraticBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((m_interp_q,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Ptr{Float64},Ptr{Int16},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),pointer(G.grid),pointer(G.level_M),pointer(w),pointer(x),pointer(y),pointer(G.nextid))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{CCGrid,LinearBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((cc_interp_l_big,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int32},Ptr{Int32},Ptr{Int}),        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]),pointer(G.hashG))
    return y
end

function c_interpbig(xi::Array{Float64},G::NGrid{CCGrid,QuadraticBF},A::Vector{Float64})
    x 		= nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((cc_interp_q_big,"libsparse.so"),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int32},Ptr{Int32},Ptr{Int}),        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]),pointer(G.hashG))
    return y
end

__precompile__()
module SparseGrids
const GridSpecs = [(:CCGrid,:LinearBF,    :cc_bf_l, :cc_M, :cc_dM, :cc_itoj, :cc_dg)
			       (:CCGrid,:QuadraticBF, :cc_bf_q, :cc_M, :cc_dM, :cc_itoj, :cc_dg)
			       (:MaxGrid,:LinearBF,   :m_bf_l,  :m_M,  :m_dM,  :m_itoj,  :m_dg)
			       (:MaxGrid,:QuadraticBF,:m_bf_q,  :m_M,  :m_dM,  :m_itoj,  :m_dg)]

using Combinatorics
include("utils.jl")
include("univariategrids.jl")
include("grids.jl")
include("jlfuncs.jl")
include("cfuncs.jl")
include("api.jl")
include("adapt.jl")

precompile(NGrid,(UnivariateGrid,Vector{Int}))
precompile(interp,(Array{Float64},NGrid,Vector{Float64}))

export BasisFunction,
	   LinearBF,
       QuadraticBF,
	   UnivariateGrid,
	   CCGrid,
	   MaxGrid,
	   CC,
	   Max,
       NGrid,
       getW,
       interp,
	   shrink!,
	   grow!
end

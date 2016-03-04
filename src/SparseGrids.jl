__precompile__()
module SparseGrids
const GridSpecs = [(:CCGrid,:LinearBF,    :cc_bf_l, :cc_M, :cc_dM, :cc_dg)
			       (:CCGrid,:QuadraticBF, :cc_bf_q, :cc_M, :cc_dM, :cc_dg)
				   (:CC2Grid,:LinearBF,    :cc_bf_l, :cc2_M, :cc2_dM, :cc2_dg)]

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
	   GridType,
	   CCGrid,
	   CC2Grid,
	   MaxGrid,
	   CC,
	   CC2,
	   Max,
       NGrid,
       getW,
       interp,
	   shrink!,
	   grow!
end

__precompile__()
module SparseGrids
const GridSpecs = [(:CCGrid,:Linear,    :cc_bf_l, :cc_M, :cc_dM, :cc_dg)
			       (:CCGrid,:Quadratic, :cc_bf_q, :cc_M, :cc_dM, :cc_dg)]

using Combinatorics,Base.Threads
include("utils.jl")
include("univariategrids.jl")
include("grids.jl")
include("jlfuncs.jl")
include("cfuncs.jl")
include("api.jl")

precompile(NGrid,(Vector{Int},))



export BasisFunction,
	   Linear,
       Quadratic,
	   UnivariateGrid,
	   GridType,
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

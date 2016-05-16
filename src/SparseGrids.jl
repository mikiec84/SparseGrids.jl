__precompile__()
module SparseGrids


using Combinatorics,Base.Threads
include("utils.jl")
include("univariategrids.jl")
include("grids.jl")
include("jlfuncs.jl")
# include("cfuncs.jl")
include("api.jl")
include("analyse.jl")
include("adapt.jl")
include("newinterp.jl")

precompile(NGrid,(Vector{Int},))



export BasisFunction,
	   Linear,
       Quadratic,
	   UnivariateGrid,
	   GridType,
	   CCGrid,
	   CC,
	   Max,
       NGrid,
       getW,
       interp,
	   level,
	   shrink!,
	   grow!
end

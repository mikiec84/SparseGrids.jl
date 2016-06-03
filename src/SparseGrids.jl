__precompile__()
module SparseGrids


using Combinatorics,Base.Threads
include("utils.jl")
include("univariategrids.jl")
include("grids.jl")
include("jlfuncs.jl")
include("adapt.jl")
include("api.jl")


precompile(NGrid,(Vector{Int},))



export Linear,
       Quadratic,
       getW,
       interp,
	   level,
	   shrink!,
	   grow!,
	   @threads,
       NGrid
end

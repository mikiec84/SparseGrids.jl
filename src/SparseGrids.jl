module SparseGrids
import Base.show
export CurtisClenshaw,
	   NoBoundary,
	   Maximum,
	   interp,
	   getW,
	   getWinv,
	   getQ,
	   grow!,
	   shrink!

include("curtisclenshaw.jl")
include("noboundary.jl")
include("maximum.jl")

getWinvC(G::CurtisClenshaw.Grid) = CurtisClenshaw.getWinvC(G)
getWinvC(G::NoBoundary.Grid) = NoBoundary.getWinvC(G)
getWinvC(G::Maximum.Grid) = Maximum.getWinvC(G)
getQ(xi::Array{Float64,2},G::CurtisClenshaw.Grid,) = CurtisClenshaw.getQ(xi,G)
getQ(xi::Array{Float64,2},G::NoBoundary.Grid) = NoBoundary.getQ(xi,G)
getQ(xi::Array{Float64,2},G::Maximum.Grid) = Maximum.getQ(xi,G)

for g in (CurtisClenshaw,Maximum,NoBoundary)
	@eval interp(xi::Array{Float64,2},G::$g.Grid,A::Array{Float64,1}) 	= $g.interp(xi,G,A)
	@eval getW(G::$g.Grid,A::Array{Float64}) 							= $g.getW(G,A)
	@eval getWinv(G::$g.Grid) 											= $g.getWinv(G)
	@eval grow!(G::$g.Grid,id,bounds::Vector{Int}=15*ones(Int,G.d))		= $g.grow!(G,id,bounds)
	@eval shrink!(G::$g.Grid,id::Vector{Bool}) 							= $g.shrink!(G,id)
	@eval show(io::IO,G::$g.Grid) 										= println(io,typeof(G)," {",G.d,",",G.q,"}[",G.n,"]")
end


end

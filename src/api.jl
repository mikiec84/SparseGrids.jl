getW(G::NGrid,A::Array{Float64}) = length(G.L)>2 ? jl_getWbig(G,A) : c_getW(G,A)
(G::NGrid)(A::Array{Float64},x::Array{Float64,2}) = length(G.L)>2 ? c_interpbig(x,G,A) : c_interp(x,G,A)

# getW(G::NGrid,A::Array{Float64}) = c_getW(G,A)
# (G::NGrid)(A::Array{Float64},x::Array{Float64,2}) = c_interp(x,G,A)

# interp{T<:GridType}(xi::Array{Float64},G::NGrid{T,GaussianRadialBF},A::Array{Float64}) = jl_interp(xi,G,A)

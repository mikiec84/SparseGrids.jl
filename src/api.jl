getW(G::NGrid,A::Array{Float64}) = length(G.L)>2 ? getW1(G,A) : c_getW(G,A)


# interp{T<:GridType,BT<:BasisFunction}(xi::Array{Float64},G::NGrid{T,BT},A::Vector{Float64}) = length(G)>250 && length(G.L)>1? c_interpbig(xi,G,A) : c_interp(xi,G,A)
interp{T<:GridType,BT<:BasisFunction}(xi::Array{Float64},G::NGrid{T,BT},A::Array{Float64}) = length(G)>250 && length(G.L)>1 ? c_interpbig(xi,G,A) : c_interp(xi,G,A)


getW(G::NGrid{CC2Grid,LinearBF},A::Vector{Float64}) = jl_getW(G,A)
interp(xi::Array{Float64},G::NGrid{CC2Grid,LinearBF},A::Array{Float64,1}) = jl_interp(xi,G,A)
interp(xi::Array{Float64},G::NGrid{CC2Grid,LinearBF},A::Array{Float64,2}) = jl_interp(xi,G,A)

interp{T<:GridType}(xi::Array{Float64},G::NGrid{T,GaussianRadialBF},A::Array{Float64}) = jl_interp(xi,G,A)



(G::NGrid)(A::Array{Float64},x::Array{Float64,2}) = interp(x,G,A)

interp(xi::Array{Float64},G::NGrid,A::Vector{Float64}) = length(G)<1500 ? c_interp(xi,G,A) : c_interpbig(xi,G,A)

interp(xi::Array{Float64},G::NGrid,A::Array{Float64,2}) = length(G)<1500 ? c_interp(xi,G,A) : c_interpbig(xi,G,A)

getW(G::NGrid,A::Vector{Float64}) = c_getW(G,A)
getW(G::NGrid,A::Array{Float64,2}) = c_getW(G,A)

getW(G::NGrid,A::Array{Float64}) = jl_getWbig(G,A)
(G::NGrid)(A::Array{Float64},x::Array{Float64,2}) = jl_interp(G,A,x)

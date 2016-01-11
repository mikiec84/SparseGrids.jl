import SparseGrids:cc_dM,cc_M,cc_g,comb

cc_trapw(i::Int) = i==1 ? [1.0] : [.5;ones(cc_M(i)-2);.5]/cc_M(i)
function cc_dtrapw(i::Int)
    if i==1
        return [1.0]
    elseif i==2
        return [.5;.5]/3
    else
        return ones(cc_dM(i))/cc_M(i)
    end
end
trap(i::Int) = cc_g(i),cc_trapw(i)


function SmolyakWeights(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    g = Vector{Float64}[ug.g(l) for l = 1:maximum(L)+D-1]
    w = Vector{Float64}[cc_trapw(l) for l = 1:maximum(L)+D-1]
    G = Array(Array{Float64},0)
    W = Array(Array{Float64},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(g[covering]))
                push!(W, vec(prod(kron(w[covering]),2)*(-1)^(maximum(mL)-l)*binomial(D-1,D+l-maximum(mL)-1)) )
            end
        end
    end
    G = vcat(G...)
    W = vcat(W...)
    uG = unique(G,1)
    uW = zeros(size(uG,1))
    for i = 1:size(uG,1)
        ids = vec(all(G.==uG[i,:]',2))
        uW[i]=sum(W[ids])
    end
    return uG,uW
end


function SmolyakWeights1(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    g = Vector{Float64}[ug.dg(l) for l = 1:maximum(L)+D-1]
    w = Vector{Float64}[cc_dtrapw(l) for l = 1:maximum(L)+D-1]
    G = Array(Array{Float64},0)
    W = Array(Array{Float64},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(g[covering]))
                push!(W, vec(prod(kron(w[covering]),2)*(-1)^(maximum(mL)-l)*binomial(D-1,D+l-maximum(mL)-1)) )
            end
        end
    end
    G = vcat(G...)
    W = vcat(W...)
    return G,W
end



for i = 1:10
x,w=SmolyakWeights(CC,[i,i])
println(dot(f(x),w))
end


# cc_trapw(i::Int) = i==1 ? [1.0] : [.5;ones(cc_M(i)-2);.5]/cc_M(i)
# function cc_dtrapw(i::Int)
#     if i==1
#         return [1.0]
#     elseif i==2
#         return [.5;.5]/3
#     else
#         return ones(cc_dM(i))/cc_M(i)
#     end
# end



function cc_simpsonsw(i::Int,j::Int)
    @assert 1≤j≤CC.M(i)
    if i==1
        return 1.0
    end
    w = 1/(CC.M(i)-1)/3
    if j==1 || j==CC.M(i)
        return w
    else
        w*= mod(j,2)==0 ? 4 : 2
    end
    return w
end

cc_simpsonsw(i::Int) = Float64[cc_simpsonsw(i,j) for j = 1:CC.M(i)]

function cc_dsimpsonsw(i::Int)
    if i==1
        return cc_simpsonsw(i)
    elseif i==2
        return cc_simpsonsw(i)-[0;1.0;0]
    else
        return cc_simpsonsw(i)-vec(Float64[cc_simpsonsw(i-1)';zeros(1,CC.M(i-1))])[1:2CC.M(i-1)-1]
    end
end

function SmolyakWeights(ug::UnivariateGrid,L::Vector{Int},mL::UnitRange{Int}=0:maximum(L))
    D = length(L)
    g = Vector{Float64}[ug.g(l) for l = 1:maximum(L)+D]
    dw = Vector{Float64}[cc_dsimpsonsw(l) for l = 1:maximum(L)+D]
    G = Array(Array{Float64},0)
    W = Array(Array{Float64},0)
    for l = mL
        for covering in comb(D,D+l)
            if all(covering.≤L+1)
                push!(G,kron(g[covering]))
                push!(W, vec(prod(kron(dw[covering]),2)))
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

type point
    x::Vector{Float64}
    l::Vector{Int}
    j::Vector{Int}
    Δj::Vector{Int}
end

function position(x::Float64)
    l=0
    if x==0.5
        l =  1
    elseif x==0.0 || x==1.0
        l = 2
    else
        for l = 3:12
            mod(x,1/2^(l-1))==0.0 && break
        end
    end
    j = 0
    dj = 0
    if l==1
        j = 1
        dj = 1
    elseif l==2
        if x==0.0
            j = 1
            dj = 1
        elseif x==1.0
            j = 3
            dj = 2
        end
    else
        j = Int(div(x,1/2^(l-1))+1)
        dj = div(j-1,2)+1
    end
    return l,j,dj
end

function getstruct{D}(G::NGrid{D})
    out = Array(point,length(G))
    for i = 1:length(G)
        out[i] = point([],[],[],[])
        for j = 1:D
            x = G.grid[i,j]
            l,j,dj = position(x)
            push!(out[i].x,x)
            push!(out[i].l,l)
            push!(out[i].j,j)
            push!(out[i].Δj,dj)
        end
    end
    return out
end



level(n::point) = sum(n.l)-length(n.l)

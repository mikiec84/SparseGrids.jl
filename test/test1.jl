function interpsafe(xi::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64})
    w = getW(G,A)
    nx = size(xi,1)
    y = zeros(nx)
	Q = []

    for i = 1:nx
        for ii = 1:G.n
            temp2 = 1.0
            for d = 1:G.d
                temp2 *= CurtisClenshaw.basis_func(xi[i,d],G.grid[ii,d],G.lvl_s[ii,d])
				temp2==0.0 ? break : nothing
            end
			i==nx && temp2>0.0 ? push!(Q,[ii;temp2]) : nothing
            y[i]+= temp2*w[ii]
        end
    end
    y,Q
end

function interpunsafe(xi::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64})
    w = getW(G,A)
    nx = size(xi,1)
    y = zeros(nx)
	Q = []

    for i = 1:nx
		ii=1
        while ii<G.n
            temp2 = 1.0
            for d = 1:G.d
                temp2 *= CurtisClenshaw.basis_func(xi[i,d],G.grid[ii,d],G.lvl_s[ii,d])
				temp2==0.0 ? break : nothing
            end
			temp2>0 ? push!(Q,[ii;temp2]) : nothing
            y[i]+= temp2*w[ii]
            println(ii," ",temp2==0 ? ii+1 : G.nextid[ii])
			temp2 == 0.0 ? ii+=1 : ii=G.nextid[ii]
        end
    end
    y,Q
end
import SparseGrids.CurtisClenshaw:nodes,Index,getnextid

id = 8
bounds = [10,10]




# Gb =CurtisClenshaw.Grid(2,2)
G =CurtisClenshaw.Grid(2,2)
for i = 1:10
    A = rand(G.n)
    id = reverse(sortperm(abs(getW(G,A))-1e100*!G.active))
    CurtisClenshaw.grow!(G,id[1],[10,10])
end


A = rand(G.n)
nothing



xi = rand(1,2)
w = getW(G,A)

Q = zeros(G.n)
Q2 = zeros(G.n)
ids =Int[]

y = 0
cnt =0
for ii = 1:G.n
    cnt +=1
    temp2 = 1.0
    for d = 1:G.d
        temp2 *= CurtisClenshaw.basis_func(xi[1,d],G.grid[ii,d],G.lvl_s[ii,d])
        temp2==0.0 ? break : nothing
    end
    temp2>0 ? push!(ids,ii) : nothing
    # println(ii,"->",temp2)
    y+= temp2*w[ii]
    Q[ii] = y
end


y2 = 0
cnt2=0
ii=1
while ii<=G.n
    cnt2+=1
    temp2 = 1.0
    for d = 1:G.d
        temp2 *= CurtisClenshaw.basis_func(xi[1,d],G.grid[ii,d],G.lvl_s[ii,d])
        temp2==0.0 ? break : nothing
    end
    # println(ii,"->",temp2)
    y2+= temp2*w[ii]
    Q2[ii] = y2
    temp2 == 0.0 ? ii+=1 : ii=G.nextid[ii]
    cnt2>G.n ? break : nothsing
end


[Q Q2][ids,:]
@assert y==y2
println(cnt,"  ",cnt2)

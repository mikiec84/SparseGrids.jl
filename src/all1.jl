include("combinatorics.jl")
include("conversion.jl")

# function grid_kron_size(dim::Int64,seq::Int64)
# 	return dMi(seq)
# end


function grid_kron_size(dim::Int64,seq::Array{Int64,1})
	out = 1
	for i = 1:dim
		out *= dMi(seq[i])
	end
	return out
end

function grid_kron(seq::Vector{Int})
	dim = length(seq)
	vs = [dXi(seq[i]) for i = 1:dim]
    sz = map(dMi,seq)
    out = Array(Float64,prod(sz),dim)
    s = 1
    for d=1:dim
        snext = s*sz[d]
        for j = 1:prod(sz)
            out[j,d] = vs[d][div(rem(j-1, snext), s)+1]
        end
        s = snext
    end
    return out
end

function node_size(dim::Int,q::Int)
	q = q+dim+1
	cnt = 0
	for q1=dim:q-1
		inds = comb(dim,q1)
		Ninds = length(inds)
		for j = 1:Ninds
			tseq=inds[j]
			cnt = cnt+grid_kron_size(dim,tseq)
		end
	end
	return cnt
end

function nodes(dim::Int64,Q::Union(Vector{Int},Int) )
	if typeof(Q)==Int
		Q = fill!(Array(Int,dim),Q)
	end
	q = maximum(Q)+dim
	nG = node_size(dim,maximum(Q))
	GRID = zeros(nG,dim)
	cnt = 0
	next_ind = 0
	for q1 = dim:q
		inds = comb(dim,q1)
		bQ = Array(Bool,length(inds))
		for i = 1:length(inds)
			bQ[i]=all(inds[i].<=Q+1)
		end
		inds = inds[bQ]
		Ninds = length(inds)
		for j = 1:Ninds
			tseq = inds[j]
			new_grid_size = grid_kron_size(dim,tseq)
			GRID[next_ind+1:next_ind+new_grid_size,:] = grid_kron(tseq)
			next_ind += new_grid_size
		end
	end
	return GRID[1:next_ind,:]
end


function addnodes(dim::Int64,Q::Union(Vector{Int},Int) )
    if typeof(Q)==Int
        Q = fill!(Array(Int,dim),Q)
    end
    q = maximum(Q)+dim
    nG = node_size(dim,maximum(Q))
    GRID = zeros(nG,dim)
    cnt = 0
    next_ind = 0

    inds = comb(dim,q)
    bQ = Array(Bool,length(inds))
    for i = 1:length(inds)
        bQ[i]=all(inds[i].<=Q+1)
    end
    inds = inds[bQ]
    Ninds = length(inds)
    for j = 1:Ninds
        tseq = inds[j]
        new_grid_size = grid_kron_size(dim,tseq)
        GRID[next_ind+1:next_ind+new_grid_size,:] = grid_kron(tseq)
        next_ind += new_grid_size
    end

    return GRID[1:next_ind,:]
end




function addgrid(d::Int,q::Union(Int,Vector{Int}),bounds=0)
    if bounds==0
        bounds =[zeros(1,d);ones(1,d)];
    elseif (size(bounds,1)!=2 || size(bounds,2)!=d)
        bounds =[zeros(1,d);ones(1,d)];
    end

    grid = addnodes(d,q)
    n = size(grid,1)
    index = getind(grid,maximum(q))
    level = getlevel(index)

    lvl_l=[[findfirst(level.==i) for i = 1:maximum(q)];n+1]
    lvl_s = convert(Array{Float64},map(Mi,index))

    grid = nUtoX(grid,bounds)
    G = Grid(d,maximum(q),n,grid,index,level,lvl_l,lvl_s,bounds,ones(Bool,n))
    return G
end

function getlevel(index::Array{Int64})
	nG = size(index,1)
	D  = size(index,2)
	Gq = fill!(Array(Int64,nG),0)
	for i = 1:nG
		tsum =0
		for d=1:D
			tsum+=index[i,d]
		end
		Gq[i] = tsum-D
	end
	return Gq
end


type Grid
	d::Int64
	q::Int64
	n::Int64
	grid::Array{Float64,2}
	index::Array{Int64,2}
	level::Vector{Int64}
	lvl_l::Vector{Int32}
	lvl_s::Array{Int32,2}
	bounds::Array{Float64,2}
	active::Vector{Bool}
end

function Grid(d::Int,q::Union(Int,Vector{Int}),bounds=0)
	if bounds==0
		bounds =[zeros(1,d);ones(1,d)];
	elseif (size(bounds,1)!=2 || size(bounds,2)!=d)
		bounds =[zeros(1,d);ones(1,d)];
	end

	grid = nodes(d,q)
	n = size(grid,1)
	index = getind(grid,maximum(q))
	level = getlevel(index)

	lvl_l=[[findfirst(level.==i) for i = 1:maximum(q)];n+1]
	lvl_s = convert(Array{Float64},map(Mi,index))

	grid = nUtoX(grid,bounds)
	G = Grid(d,maximum(q),n,grid,index,level,lvl_l,lvl_s,bounds,ones(Bool,n))
	return G
end



function interp(xi1::Array{Float64},G::Grid,A::Array{Float64,1})
	xi = nXtoU(xi1,G.bounds)
	xi[xi.>1]=1.0
	xi[xi.<0]=0.0
	nx = size(xi,1)

	Aold = zeros(G.n)
	dA = zeros(G.n)
	w = zeros(G.n)

	grid = nXtoU(G.grid,G.bounds)
    ccall((libwget, "libsparse.so"),
		Void,
		(Ptr{Float64},Int32,Int32,Ptr{Int64},Ptr{Int32},Ptr{Float64},Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		pointer(grid),G.n,G.d,pointer(G.lvl_s),pointer(G.lvl_l),pointer(A),G.q,pointer(Aold),pointer(dA),pointer(w))
	xold = zeros(nx)
	dx = zeros(nx)
	ccall((libinterp, "libsparse.so"),
		Void,
		(Ptr{Float64},Int32,Ptr{Float64},Int32,Int32,
		Ptr{Float64},Ptr{Float64},Ptr{Float64},Int32,Ptr{Float64},
		Ptr{Float64},Ptr{Float64}),
		pointer(xi),nx,pointer(grid),G.n,G.d,
		pointer(G.lvl_s),pointer(G.lvl_l),pointer(A),G.q,pointer(w),
		pointer(xold),pointer(dx))
	yi = xold*1.0
	return yi
end



function getW(G::Grid,A)
	Aold = zeros(G.n)
	dA = zeros(G.n)
	w = zeros(G.n)
	grid = nXtoU(G.grid,G.bounds)
	ccall((libwget, "libsparse.so"),
		Void,
		(Ptr{Float64},Int32,Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64},Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		pointer(grid),G.n,G.d,pointer(G.lvl_s),pointer(G.lvl_l),pointer(A),G.q,pointer(Aold),pointer(dA),pointer(w))
	return w
end



function getWinv(G::Grid)
	drange=collect(1:G.d)
	grid = nXtoU(G.grid,G.bounds)
	qA= eye(G.n)
	qAold = zeros(G.n,G.n)
	qw= zeros(G.n,G.n)
	qtemp = zeros(G.n)

	for i = 1:G.lvl_l[1]-1
		@simd for j = 1:G.n
		    @inbounds qw[i,j] = qA[i,j] - qAold[i,j]
		end
    end
    for i in 1:G.n
    	fill!(qtemp,0)
        for ii in 1:G.lvl_l[1]-1
            temp2=1.0
            for d in drange
                @inbounds temp2 *= basis_func(grid[i,d],grid[ii,d],G.lvl_s[ii,d])
            end
            @simd for j = 1:G.n
	            @inbounds qtemp[j] += temp2*qw[ii,j]
	        end
        end

	    @simd for j = 1:G.n
		    @inbounds qAold[i,j] += qtemp[j]
		end
    end

	for q in 1:G.q
	    for i = G.lvl_l[q]:G.lvl_l[q+1]-1
	    	@simd for j = 1:G.n
		        @inbounds qw[i,j] = qA[i,j] - qAold[i,j]
		    end
	    end
	    for i in G.lvl_l[q]:G.n
	        fill!(qtemp,0)
	        for ii in G.lvl_l[q]:G.lvl_l[q+1]-1
	            temp2=1.0
	            for d in drange
	                @inbounds temp2 *= basis_func(grid[i,d],grid[ii,d],G.lvl_s[ii,d])
	            end
	            @simd for j = 1:G.n
		            @inbounds qtemp[j] += temp2*qw[ii,j]
		        end
	        end
	        @simd for j = 1:G.n
			    @inbounds qAold[i,j] += qtemp[j]
			end
	    end
	end
   return sparse(qw)
end


function getWinvC(G::Grid)
	Aold = zeros(G.n,G.n)
	dA = zeros(G.n,G.n)
	w = zeros(G.n,G.n)
	A = eye(G.n,G.n)
	grid = deepcopy(G.grid)
	ccall((libwgetinv, "libsparse.so"),
		Void,
		(Ptr{Float64},Int32,Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64},Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		pointer(grid),G.n,G.d,pointer(G.lvl_s),pointer(G.lvl_l),pointer(A),G.q,pointer(Aold),pointer(dA),pointer(w))
	return sparse(w)
end



function getQ(xi1::Array{Float64},G::Grid)
    xi = nXtoU(xi1,G.bounds)
    xi[xi.>1]=1.0
    xi[xi.<0]=0.0
    nx = size(xi,1)
    grid = nXtoU(G.grid,G.bounds)

    lvl_l = [1;G.lvl_l]

    Q = spzeros(nx,G.n)
    drange=collect(1:G.d)
    for i in 1:nx
        for q in 0:G.q
            for ii in lvl_l[q+1]:lvl_l[q+2]-1
                temp2=1.0
                for d in drange
                    @inbounds temp2 *= basis_func(xi[i,d],grid[ii,d],G.lvl_s[ii,d])

                end
                Q[i,ii]+=temp2
            end
        end
    end
   return Q
end






function adapt!(G::Grid,id::Int;bounds=[100])
	@assert id≤G.n
    G.active[id]=false
    targ = G.grid[id,:]
    nX = addgrid(G.d,min(G.level[id]+1,bounds),G.bounds)
    nearest=Any[]

    n = G.d*2
    for d=1:G.d
        if targ[d]==G.bounds[1,d] || targ[d]==G.bounds[2,d]
            n-=1
        end
    end

	# dst = Float64[norm(nX.grid[i,:]-targ)*(1000-999*(prod(nX.grid[i,:]-targ)==0)) for i = 1:nX.n]
    dst = Float64[norm(nX.grid[i,:]-targ) for i = 1:nX.n]

    grid = nX.grid[sortperm(dst)[1:n],:]

    id1 = ones(Bool,size(grid,1))
    for i = 1:size(grid,1)
        for j = 1:G.n
            if grid[i,:]==G.grid[j,:]
                id1[i] = false
            end
        end
    end
    grid = [G.grid;grid[id1,:]]
    active = [G.active;ones(Bool,length(id1))]
    n = size(grid,1)
    index = getind(nXtoU(grid,G.bounds),maximum(min(G.level[id]+1,bounds)))
    level = getlevel(index)
    id2 = sortperm(level)
    grid = grid[id2,:]
    active = active[id2]

    index = getind(grid,maximum(min(G.level[id]+1,bounds)))
    level = getlevel(index)

    lvl_l=[[findfirst(level.==i) for i = 1:maximum(level)];n+1]
    lvl_s = convert(Array{Float64},map(Mi,index))

    # grid = nUtoX(grid,G.bounds)
    return Grid(G.d,maximum(level),n,grid,index,level,lvl_l,lvl_s,G.bounds,active)
end

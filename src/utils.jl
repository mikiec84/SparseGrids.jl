
function XtoX(X::Array{Float64,1},bnds1::Array{Float64,1},bnds2::Array{Float64,1})
	lb1 = minimum(bnds1)
	ub1 = maximum(bnds1)
	lb2 = minimum(bnds2)
	ub2 = maximum(bnds2)
	out = Array(Float64,length(X))
	for i = 1:length(X)
		out[i] = (X[i]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
	end
	return out
end

function nXtoX(X::Array{Float64,2},bnds1::Array{Float64,2},bnds2::Array{Float64,2})
	N,R = size(X)
	out = Array(Float64,size(X))
	for r = 1:R
		lb1 = minimum(bnds1[:,r])
		ub1 = maximum(bnds1[:,r])
		lb2 = minimum(bnds2[:,r])
		ub2 = maximum(bnds2[:,r])
		for i = 1:N
			out[i,r] = (X[i,r]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
		end
	end
	return out
end
nXtoX(X::Array{Float64,1},bnds1::Array{Float64,2},bnds2::Array{Float64,2})=XtoX(X,bnds1[:],bnds2[:])

UtoX(U::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(U,[0.0,1.0],bnds)
XtoU(X::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(X,bnds,[0.0,1.0])

nUtoX(U::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(U,[zeros(1,size(U,2));ones(1,size(U,2))],bnds)
# nXtoU(X::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(X,bnds,[zeros(1,size(X,2));ones(1,size(X,2))])
nXtoU(X::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(X,bnds[:],[0.0,1.0])
nUtoX(U::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(U,[0.0,1.0],bnds[:])


function nXtoU(X::Array{Float64,2},bnds1::Array{Float64,2})
	N,R = size(X)
	out = Array(Float64,size(X))
	for r = 1:R
		lb1 = minimum(bnds1[:,r])
		ub1 = maximum(bnds1[:,r])
		@inbounds for i = 1:N
			out[i,r] = (X[i,r]-lb1)/(ub1-lb1)
		end
	end
	return out
end

function hsh(x::Vector{Int})
	h = x[1]*17
	for i = 2:length(x)
		h+=x[i]
		h*=17
	end
	h
end

function hsh(x::Array{Int,2})
	D = size(x,2)
	H = zeros(Int,size(x,1))
	for i = 1:size(x,1)
		H[i] = x[i,1]*17
		for j = 2:D
			H[i] += x[i,j]
			H[i] *= 17
		end
	end
	return H
end


function subs!(x::Expr,s::Pair)
    for i = 1:length(x.args)
        if x.args[i]==s.first
            x.args[i] = s.second
        elseif isa(x.args[i],Expr)
            subs!(x.args[i],s)
        end
    end
end

function subs!(x::Expr,list::Dict)
    for i = 1:length(x.args)
        if in(x.args[i],keys(list))
            x.args[i] = list[x.args[i]]
        elseif isa(x.args[i],Expr)
            subs!(x.args[i],list)
        end
    end
end


function _threadsforfixed(fixed,iter,lbody)
    fun = gensym("_threadsfor")
    lidx = iter.args[1]         # index
    range = iter.args[2]

    svars = map(gensym,fixed.args)
    newvar = quote end
    for i = 1:length(fixed.args)
        push!(newvar.args,:($(svars[i])=copy($(fixed.args[i]))))
    end
    subs!(lbody,Dict(zip(fixed.args,svars)))

    quote
        function $fun()
            tid = threadid()
            r = $(esc(range))
            # divide loop iterations among threads
            len, rem = divrem(length(r), nthreads())
            # not enough iterations for all the threads?
            if len == 0
                if tid > rem
                    return
                end
                len, rem = 1, 0
            end
            # compute this thread's iterations
            f = 1 + ((tid-1) * len)
            l = f + len - 1
            # distribute remaining iterations evenly
            if rem > 0
                if tid <= rem
                    f = f + (tid-1)
                    l = l + tid
                else
                    f = f + rem
                    l = l + rem
                end
            end
            # run this thread's iterations

            $(esc(quote
            $newvar
            end))
            for i = f:l
                local $(esc(lidx)) = Base.unsafe_getindex(r,i)
                $(esc(lbody))
            end
        end
        ccall(:jl_threading_run, Void, (Any,), Core.svec($fun))
    end
end


import Base.Threads.@threads
macro threads(args...)
    na = length(args)
    if na == 1
        ex = args[1]
        if !isa(ex, Expr)
            throw(ArgumentError("need an expression argument to @threads"))
        end
        if is(ex.head, :for)
            return Base.Threads._threadsfor(ex.args[1],ex.args[2])
        else
            throw(ArgumentError("unrecognized argument to @threads"))
        end
    elseif na==2
        fixed = args[1]
        ex = args[2]
        if is(ex.head, :for)
            return _threadsforfixed(fixed,ex.args[1],ex.args[2])
        else
            throw(ArgumentError("unrecognized argument to @threads"))
        end
    else
        throw(ArgumentError("wrong number of arguments (1 or 2)"))
    end
end


function comb(D::Int, Q::Int)
	D==Q ? (return Vector{Int}[ones(Int,D)]) : nothing
	L = Q - D + 1
	out = Array(Vector{Int},binomial(Q-1, D-1))
	tL1 = ones(Int, D)
	tL2 = L*ones(Int, D)
	p = 1
	cnt = 0

	while tL1[D] < L
		tL1[p] += 1
		if tL1[p] > tL2[p]
			tL1[p] = 1
			p += 1
		else
			for i = 1:p-1
				tL2[i] = tL2[p] - tL1[p] + 1
			end
			p = 1
			tL1[1] = tL2[1]
			cnt += 1
			out[cnt] = copy(tL1)
		end
	end
	return out
end


ndgrid(v::AbstractVector) = (copy(v),)

function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end

function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end

function ndgrid{T}(vs::AbstractVector{T}...)
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array(T, sz),n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

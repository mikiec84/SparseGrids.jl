import SparseGrids.CurtisClenshaw:M,dM,Xi,dXi,xi,dxi,Index,Node


function hsh(x::Vector{Int})
	h = x[1]*17
	for i = 2:length(x)
		h+=x[i]
		h*=17
	end
	h
end

function indexshift(i::Int,j::Int,q::Int)
	if i ==1
		return div(M(q)-1,2)+1
	elseif i==2
		if j==1
			return 1
		else
			return M(q)
		end
	else
		return div((M(q)-1)*(j*2-1),M(i)-1)+1
	end
end

function interp1(x::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64},T::Dict{Vector{Int},Int})
	w 		= getW(G,A)
	y 		= zeros(size(x,1))
	levels	= [UnitRange(G.lvl_l[i-1]+0,G.lvl_l[i]-1) for i = 2:length(G.lvl_l)]
	combs 	= vcat(G.index[1:1,:],[unique(G.index[i,:],1) for i in levels]...)
	nc 		= size(combs,1)

	for i = 1:length(y)
		J 		= ones(Int,G.q+1,G.d)
		D 		= zeros(G.q+1,G.d)
		id  	= zeros(Int,G.d)
		for d = 1:G.d
			for q = 1:G.q+1
				j 	= clamp(round(Int,x[i,d]*(dM(q))+1/2),1,dM(q))
				D[q,d] 	= abs(1.0-abs(x[i,d]-dxi(q,j))*(M(q)-1.0))
				J[q,d]  = indexshift(q,j,G.q+1)
			end
		end

		for ii = 1:nc
			b = 1.0
			for d = 1:G.d
				b*=D[combs[ii,d],d]
				id[d]=J[combs[ii,d],d]
			end
			y[i]+=b*w[T[id]]
		end
	end
	y
end




function interp3(x::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64},combs::Array{Int,2},combl::Vector{Tuple{Int,Int}},hGji::Array{UInt})
	w 		= getW(G,A)
	y 		= zeros(size(x,1))
	nc 		= size(combs,1)

	for i = 1:length(y)
		J 		= ones(Int,G.q+1,G.d)
		D 		= zeros(G.q+1,G.d)
		id  	= zeros(Int,G.d)
		for d = 1:G.d
			for q = 1:G.q+1
				j 	= clamp(round(Int,x[i,d]*(dM(q))+1/2),1,dM(q))
				D[q,d] 	= abs(1.0-abs(x[i,d]-dxi(q,j))*(M(q)-1.0))
				J[q,d]  = indexshift(q,j,G.q+1)
			end
		end

		for ii = 1:nc
			b = 1.0
			for d = 1:G.d
				b*=D[combs[ii,d],d]
				id[d]=J[combs[ii,d],d]
			end
			hid = hash(id)
			id1 = 0
			for iii = combl[ii][1]:combl[ii][2]
				hGji[iii]==hid ? (id1=iii ; break) : nothing
			end
			y[i]+=b*w[id1]
		end
	end
	y
end


function interp31(x::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64},combs::Array{Int,2},combl::Vector{Tuple{Int,Int}},hGji::Array{Int})
	w 		= getW(G,A)
	y 		= zeros(size(x,1))
	nc 		= size(combs,1)

	for i = 1:length(y)
		J 		= ones(Int,G.q+1,G.d)
		D 		= zeros(G.q+1,G.d)
		id  	= zeros(Int,G.d)
		for d = 1:G.d
			for q = 1:G.q+1
				j 	= clamp(round(Int,x[i,d]*(dM(q))+1/2),1,dM(q))
				D[q,d] 	= abs(1.0-abs(x[i,d]-dxi(q,j))*(M(q)-1.0))
				J[q,d]  = indexshift(q,j,G.q+1)
			end
		end

		for ii = 1:nc
			b = 1.0
			hid = 0
			for d = 1:G.d
				b*=D[combs[ii,d],d]
				hid +=J[combs[ii,d],d]
				hid *=17
			end
			id1 = 0
			for iii = combl[ii][1]:combl[ii][2]
				hGji[iii]==hid ? (id1=iii;break) : nothing
			end
			println(b)
			y[i]+=b*w[id1]
		end
	end
	y
end


Q = CurtisClenshaw.Index(5,5)
G = CurtisClenshaw.Grid(length(Q),Q.x);
A = rand(G.n)
levels	= [UnitRange(G.lvl_l[i-1]+0,G.lvl_l[i]-1) for i = 2:length(G.lvl_l)]
combs 	= vcat(G.index[1:1,:],[unique(G.index[i,:],1) for i in levels]...)
combl 	= [UnitRange(findfirst(all(combs[i:i,:].==G.index,2)),findlast(all(combs[i:i,:].==G.index,2))) for i = 1:size(combs,1)]
combl 	= Tuple{Int,Int}[(findfirst(all(combs[i:i,:].==G.index,2)),findlast(all(combs[i:i,:].==G.index,2))) for i = 1:size(combs,1)]
nc 		= size(combs,1)
T = Dict{Vector{Int},Int}()
Gji = zeros(Int,size(G.grid))
for i ∈ 1:G.n
   push!(T,[indexshift(G.index[i,d],clamp(round(Int,G.grid[i,d]*(dM(G.index[i,d]))+1/2),1,dM(G.index[i,d])),maximum(G.index[:,d])) for d = 1:G.d]=>i)
   Gji[i,:] = [indexshift(G.index[i,d],clamp(round(Int,G.grid[i,d]*(dM(G.index[i,d]))+1/2),1,dM(G.index[i,d])),maximum(G.index[:,d])) for d = 1:G.d]
end
hGji = vcat([hash(Gji[i,:]) for i = 1:G.n]...)
hGji1 = vcat([hsh(Gji[i,:]) for i = 1:G.n]...)


x = rand(1,G.d)
# @timecnt interp(x,G,A) 2
# @timecnt interp1(x,G,A,T) 2
# @timecnt interp3(x,G,A,combs,combl,hGji) 2
# @timecnt interp31(x,G,A,combs,combl,hGji1) 2




y=interp(x,G,A)
y3=interp3(x,G,A,combs,combl,hGji)
y31=interp31(x,G,A,combs,combl,hGji1)

# @timecnt ccall((:_Z6interplllllPlPd,"libsparse1.so"),Void,(Int32,Int32,Int32,Int32,Int32,Ptr{Int},Ptr{Float64}),G.d,G.q,G.n,size(x,1),nc,combs,x) 2




combsl = Int[(combl[i][1]) for i = 1:nc]
combsu = Int[(combl[i][2]) for i = 1:nc]
y = zeros(size(x,1))
w = getW(G,A)
ccall((:_Z6interplllllPdS_S_PlS0_S0_S0_,"libsparse1.so"),Void,
(Int32,Int32,Int32,Int32,Int32,
Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int},Ptr{Int}),
G.d,G.q,G.n,size(x,1),nc,pointer(x),pointer(w),pointer(y),pointer(combs),pointer(combsl),pointer(combsu),pointer(hGji1))



# combsl = Int[(combl[i][1]) for i = 1:nc]
# combsu = Int[(combl[i][2]) for i = 1:nc]
# ccall((:_Z5test3PlS_,"libsparse1.so"),Void,(Ptr{Int32},Ptr{Int32}),pointer(combsl),pointer(combsu))

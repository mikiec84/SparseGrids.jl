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

function interp31(x::Array{Float64},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64},combs::Array{Int,2},combl::Vector{Tuple{Int,Int}},hGji::Array{Int})
	w 		= SparseGrids.CurtisClenshaw.getW(G,A)
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
			y[i]+=b*w[id1]
		end
	end
	y
end

function interp31CC(x::Array{Float64,2},G::SparseGrids.CurtisClenshaw.Grid,A::Vector{Float64},nc::Int,combs::Array{Int},combsl::Vector{Int},combsu::Vector{Int},hGji::Vector{Int})
	y = zeros(size(x,1))
	w = getW(G,A)
	ccall((:_Z6interplllllPdS_S_PlS0_S0_S0_,"libsparse1.so"),Void,
	(Int32,Int32,Int32,Int32,Int32,
	Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Int},Ptr{Int}),
	G.d,G.q,G.n,size(x,1),nc,pointer(x),pointer(w),pointer(y),pointer(combs),pointer(combsl),pointer(combsu),pointer(hGji))
	return y
end

# TR = zeros(8,6,3)
# for Gd ∈ 1:8, Gq ∈ 1:6,nT ∈ 1:3

n = 1000
TR = zeros(n)
Td = zeros(Int,n)
Tq = zeros(Int,n)
Tn = zeros(Int,n)
nT = 10
for it = 1:n
	Gd = rand(2:7)
	Gq = rand(1:4,Gd)
	Td[it] = Gd
	Tq[it] = maximum(Gq)
	G = SparseGrids.CurtisClenshaw.Grid(Gd,Gq);
	Tn[it]=G.n
	A = rand(G.n)
	levels	= [UnitRange(G.lvl_l[i-1]+0,G.lvl_l[i]-1) for i = 2:length(G.lvl_l)]
	combs 	= vcat(G.index[1:1,:],[unique(G.index[i,:],1) for i in levels]...)
	nc 		= size(combs,1)
	combl 	= Tuple{Int,Int}[(findfirst(all(combs[i:i,:].==G.index,2)),findlast(all(combs[i:i,:].==G.index,2))) for i = 1:size(combs,1)]
	combsl = Int[(combl[i][1]) for i = 1:nc]
	combsu = Int[(combl[i][2]) for i = 1:nc]

	Gji = zeros(Int,size(G.grid))
	for i ∈ 1:G.n
	   Gji[i,:] = [indexshift(G.index[i,d],clamp(round(Int,G.grid[i,d]*(dM(G.index[i,d]))+1/2),1,dM(G.index[i,d])),maximum(G.index)) for d = 1:G.d]
	end
	hGji = vcat([hsh(Gji[i,:]) for i = 1:G.n]...)

	x = rand(G.n*nT,G.d)


	y	= interp(x,G,A)
	# y31 = interp31(x,G,A,combs,combl,hGji)
	y31cc = interp31CC(x,G,A,nc,combs,combsl,combsu,hGji)

	@assert y31cc==y

	# tic();interp(x,G,A) ;T1=toq()
	# tic(); interp31CC(x,G,A,nc,combs,combsl,combsu,hGji) ;T2=toq()
	T1 = @timecnt interp(x,G,A) 0.25
 	T2 = @timecnt interp31CC(x,G,A,nc,combs,combsl,combsu,hGji) 0.25
	println("$Gd,$Gq  ",round(100(T1/T2-1)),"%")
	# TR[Gd,Gq,nT] = T1/T2
	TR[it] = T1/T2
end


# histogram(vec(TR[:,:,1]))

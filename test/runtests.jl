using SparseGrids, Base.Test

@time for d = 1:8,q=1:3
	tic()
	G = CurtisClenshaw.Grid(d,q)
	A = rand(G.n)
	W = getW(G,A)
	iW = getWinv(G)

	@assert maximum(abs(A - interp(G.grid,G,A)))<1e-12
	@assert maximum(abs(iW*A - W))<1e-12

	t = toq()
	println("{$d,$q}[$(G.n)]   $t")
end


@time for d = 1:4,q=1:4
	tic()
	G = NoBoundary.Grid(d,q)
	A = rand(G.n)
	W = getW(G,A)
	iW = getWinv(G)
	@assert maximum(abs(A - interp(G.grid,G,A)))<1e-12
	@assert maximum(abs(iW*A - W))<1e-12
	t = toq()
	println("{$d,$q}[$(G.n)]   $t")
end


@time for d = 1:4,q=1:2
	tic()
	G = Maximum.Grid(d,q)
	A = rand(G.n)
	W = getW(G,A)
	iW = getWinv(G)
	@assert maximum(abs(A - interp(G.grid,G,A)))<1e-12
	@assert maximum(abs(iW*A - W))<1e-12

	t = toq()
	println("{$d,$q}[$(G.n)]   $t")
end

for seq in partitions(8,3)
	G = CurtisClenshaw.Grid(3,seq,[rand(1,3);rand(1,3)+1])
	A = rand(G.n)
	W = getW(G,A)
	iW = getWinv(G)
	@assert maximum(abs(A - interp(values(G),G,A)))<1e-12
	@assert maximum(abs(iW*A - W))<1e-12
end

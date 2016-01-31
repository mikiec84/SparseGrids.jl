G = NGrid(CC,[8])
A = rand(length(G),2)
@assert maximum(abs(interp(G.grid,G,A[:,1])-A[:,1]))<1e-15
@assert maximum(abs(interp(G.grid,G,A)-A))<1e-15
G = NGrid(CC,[8],B=QuadraticBF)
@assert maximum(abs(interp(G.grid,G,A[:,1])-A[:,1]))<1e-15
@assert maximum(abs(interp(G.grid,G,A)-A))<1e-15

for d = 2:4 
    G = NGrid(CC,4*ones(Int,d))
    A = rand(length(G),2)
    @assert maximum(abs(SparseGrids.c_interp(G.grid,G,A[:,1])-A[:,1]))<1e-13
    @assert maximum(abs(SparseGrids.c_interp(G.grid,G,A)-A))<1e-13
    @assert maximum(abs(SparseGrids.c_interpbig(G.grid,G,A[:,1])-A[:,1]))<1e-13
    @assert maximum(abs(SparseGrids.c_interpbig(G.grid,G,A)-A))<1e-13

    G = NGrid(CC,4*ones(Int,d),B=QuadraticBF)
    A = rand(length(G),2)
    @assert maximum(abs(SparseGrids.c_interp(G.grid,G,A[:,1])-A[:,1]))<1e-14
    @assert maximum(abs(SparseGrids.c_interp(G.grid,G,A)-A))<1e-14
    @assert maximum(abs(SparseGrids.c_interpbig(G.grid,G,A[:,1])-A[:,1]))<1e-14
    @assert maximum(abs(SparseGrids.c_interpbig(G.grid,G,A)-A))<1e-14
end

using SparseGrids,Base.Test
# Linear Basis Function


println("1→2 Dimensions")

for d ∈ 2:6, l ∈ 2:4, B ∈ [Linear,Quadratic]
    G = NGrid(ones(Int,l)*d,B=B)
    A = rand(length(G))
    @test maximum(abs(G(A,G.grid)-A))<1e-14
end


d,l,B=4,4,Quadratic
G = NGrid(ones(Int,d)*l,B=B)
for i = 1:1000
    x = rand(1000,d)
    A = rand(length(G),d)
    @assert all(SparseGrids.c_interp(G.grid,G,A)-A.<1e14) ||
    @assert all(SparseGrids.c_interpbig(G.grid,G,A)-A.<1e14)
    maximum(abs(SparseGrids.c_interp(x,G,A)-SparseGrids.c_interpbig(x,G,A)))<1e-14
end

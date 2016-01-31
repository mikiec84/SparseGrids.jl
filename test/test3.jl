using SparseGrids
G = NGrid(CC,[5,5,5,5])
A = rand(length(G))
x = zeros(10000,4)


for i = 1:200;interp(x,G,A) ;end


10000*length(G.L)*(maximum(G.L)+1)*200*3

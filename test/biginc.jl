function c_interpbig_s(xi::Array{Float64},G::NGrid{CCGrid,Linear},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z17interp_cc_l_big_slllllPdS_S_PsPiS1_Pl,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int32},Ptr{Int32},Ptr{Int}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]),pointer(G.hashG))
    return y
end

function c_interpbig1_s(xi::Array{Float64},G::NGrid{CCGrid,Linear},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z18interp_cc_l_big1_slllllPdS_S_PsS0_PiS1_Pl,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32},Ptr{Int}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]),pointer(G.hashG))
    return y
end

function c_interpbig1_p(xi::Array{Float64},G::NGrid{CCGrid,Linear},A::Vector{Float64})
    x 		= SparseGrids.nXtoU(xi,G.bounds)
    y 		= zeros(size(xi,1))
    w 		= getW(G,A)
    ccall((:_Z18interp_cc_l_big1_plllllPdS_S_PsS0_PiS1_Pl,lsparse),
        Void,
        (Int32,Int32,Int32,Int32,Int32,
        Ptr{Float64},Ptr{Float64},Ptr{Float64},
        Ptr{Int16},Ptr{Int16},Ptr{Int32},Ptr{Int32},Ptr{Int}),
        length(G.L),maximum(G.L),size(G.grid,1),size(x,1),size(G.coverings,1),
        pointer(x),pointer(w),pointer(y),
        pointer(G.coverings),pointer(G.coverings_dM),pointer(G.coveringsloc[1]),pointer(G.coveringsloc[2]),pointer(G.hashG))
    return y
end

res = Array(Any,6,6)
for d = 2:7,q= 1:6
G = NGrid(CC,q*ones(Int,d))
A = rand(length(G))
x = rand(10000,length(G.L))

y = interp(x,G,A)
yb = SparseGrids.c_interpbig(x,G,A)
y1 = c_interpbig1_p(x,G,A)
# ys = c_interpbig_s(x,G,A)
# y1s = c_interpbig1_s(x,G,A)

b1=@benchmark SparseGrids.c_interp(x,G,A)
bbig=@benchmark SparseGrids.c_interpbig(x,G,A)
bbig1p=@benchmark c_interpbig1_p(x,G,A)
# bbigs=@benchmark c_interpbig_s(x,G,A)
# bbig1s=@benchmark c_interpbig1_s(x,G,A)
res[d-1,q] = (b1,bbig,bbig1p)
end

sizes=Float64[SparseGrids.SmolyakSize(CC,q*ones(Int,d)) for d = 2:7,q= 1:6]
rspeed=Float64[mean(res[i,j][k].samples.elapsed_times)/mean(res[i,j][1].samples.elapsed_times)  for i = 1:6,j=1:6,k = 1:3]

Plot([scatter(x=vec(sizes),y=1./vec(rspeed[:,:,2]),mode="markers"),scatter(x=vec(sizes),y=1./vec(rspeed[:,:,3]),mode="markers")],Layout(xaxis=Dict("type"=>"log","autorange"=>true)))


Gs = NGrid[NGrid(CC,rand(1:5,5)) for i = 1:50]

res = zeros(50,5)
for i = 1:50
    A = rand(length(Gs[i]))
    x = rand(5000,length(Gs[i].L))
    res[i,1] = length(Gs[i])

    tic();for j = 1:5; SparseGrids.c_interp(x,Gs[i],A);end;res[i,2]=toq()
    tic();for j = 1:5; SparseGrids.c_interpbig(x,Gs[i],A);end;res[i,3]=toq()
    tic();for j = 1:5; SparseGrids,c_interpbig1_p(x,Gs[i],A);end;res[i,4]=toq()
end

Plot([scatter(x=res[:,1],y=res[:,2],mode="markers"),
      scatter(x=res[:,1],y=res[:,3],mode="markers"),
      scatter(x=res[:,1],y=res[:,4],mode="markers")],Layout(yaxis=Dict("type"=>"log","autorange"=>true)))

Plot(scatter(x=res[:,1],y=res[:,2]./res[:,4],mode="markers"),Layout(yaxis=Dict("type"=>"log","autorange"=>true)))

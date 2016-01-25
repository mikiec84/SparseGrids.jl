import SparseGrids:cc_M,cc_dM,cc_g,cc_dg,cc_itoi
# using Benchmarks,ProfileView

function f1(xi::Array{Float64},G::NGrid,A::Vector{Float64})
   w = getW(G,A)
   x = SparseGrids.nXtoU(xi,G.bounds)
   y = zeros(size(x,1))
   nc = size(G.coverings,1)
   mL      = maximum(G.L)
   D       = length(G.L)
   for i = 1:size(x,1)
       J = ones(Int,mL+1,D)
       B = zeros(mL+1,D)
       id      = zeros(Int,D)
       for d = 1:D
           for l = 1:mL+1
               j = round(Int,x[i,d]*(SparseGrids.cc_dM(l))+1/2)
               B[l,d] = SparseGrids.cc_bf_l(x[i,d],SparseGrids.cc_dg(l,j),Int16(SparseGrids.cc_M(l)))
               J[l,d]  =SparseGrids.cc_itoi(l,j,mL+1)
           end
       end
       for ii = 1:nc
           b = 1.0
           hid = 0
           for d = 1:D
               b*=B[G.coverings[ii,d],d]
               hid +=J[G.coverings[ii,d],d]
               hid *=17
           end
           id1 = Int32(0)
           for iii = G.coveringsloc[1][ii]:G.coveringsloc[2][ii]
               G.hashG[iii]==hid ? (id1=iii;break) : nothing
           end
           y[i]+=b*w[id1]
       end
   end
   y
end

function f2(x::Array{Float64},G::NGrid,A::Vector{Float64})
    w = getW(G,A)
    y = zeros(size(x,1))
    nc = size(G.coverings,1)
    mL      = maximum(G.L)
    D       = length(G.L)

    for i = 1:size(x,1)
        J = ones(Int,mL+1,D)
        B = zeros(mL+1,D)
        for d = 1:D
            for l = 1:mL+1
            j = round(Int,x[i,d]*(SparseGrids.cc_dM(l))+1/2)
                B[l,d] = SparseGrids.cc_bf_l(x[i,d],SparseGrids.cc_dg(l,j),Int16(SparseGrids.cc_M(l)))
                J[l,d] = j
            end
        end
        for ii = 1:nc
            b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
            id2 = J[G.coverings[ii,D],D]-1
            for d = D-1:-1:2
                b*=B[G.coverings[ii,d],d]
                id2 = id2*G.coverings_M[ii,d]+(J[G.coverings[ii,d],d]-1)
            end
            id2=(J[G.coverings[ii,1],1]-1)+G.coverings_M[ii,1]*id2+1+G.coveringsloc[1][ii]-1
           y[i]+=b*w[id2]
       end
   end
   y
end




G = NGrid(CC,[5,5,5])
A = rand(length(G))
x = rand(1,length(G.L))


@time for i = 1:10;y1=f1(x,G,A);end
@time for i = 1:10;y2=f2(x,G,A);end
@time for i = 1:10;y3=f3(x,G,A);end
# b1 = @benchmark f1(x,G,A)
# b2 = @benchmark f2(x,G,A)
# @assert y1==y2
#
# Profile.clear()
# @profile f1(x,G,A)
# ProfileView.view()
#
# Profile.clear()
# @profile f2(x,G,A)
# ProfileView.view()

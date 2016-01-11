import SparseGrids: jl_interpslow,jl_interp,jl_interpbig,c_interp,c_interpbig,nXtoU
import SparseGrids:GridType,UnivariateGrid,BasisFunction,SmolyakGrid,level,level_loc,nextid,hsh,nXtoU

function interpbig1(xi::Array{Float64},G::NGrid,A::Vector{Float64})
   w = getW(G,A)
   x = nXtoU(xi,G.bounds)
   y = zeros(size(x,1))
   nc = size(G.coverings,1)
   mL      = maximum(G.L)
   D       = length(G.L)

   for i = 1:length(y)
       J = ones(Int,mL+1,D)
       B = zeros(mL+1,D)
       id      = zeros(Int,D)
       for d = 1:D
           for l = 1:mL+1
               j = clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
               B[l,d] = cc_bf_l(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
               J[l,d]  = cc_itoj(l,j,mL+1)
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
           if (G.coveringsloc[2][ii]-G.coveringsloc[1][ii])> 20
               for iii = G.coveringsloc[1][ii]:G.coveringsloc[2][ii]
                   G.hashG[iii]==hid ? (id1=iii;break) : nothing
               end
            else
               id0 = Int32(G.coveringsloc[1][ii])
               id2 = Int32(G.coveringsloc[2][ii])
               while true
                   G.hashG[id0]==hid ? (id1=id0;break) : nothing
                   G.hashG[id2]==hid ? (id1=id2;break) : nothing
                   id1=div(id0+id2,Int32(2))
                   G.hashG[id1]>hid ? (id2=id1) : (id0=id1)
               end
           end
           y[i]+=b*w[id1]
       end
   end
   y
end


L = [8,8,5]
G = NGrid(CC,L)
x=  rand(length(G),length(L))

# @timecnt y=SparseGrids.c_interp(x,G,f(G.grid)) 1
# @timecnt yb=SparseGrids.c_interpbig(x,G,f(G.grid)) 1

@timecnt SparseGrids.jl_interpbig(x,G,f(G.grid)) 1
@timecnt interpbig1(x,G,f(G.grid)) 1

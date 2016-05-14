using SparseGrids
import SparseGrids: cc_bf_l,cc_bf_q,cc_dM,cc_dg,cc_M

function interp1{D,BF}(G::NGrid{D,BF},Tp::Vector{UInt32},Ts::Array{Int},A::Vector{Float64},x::Array{Float64})
    w = getW(G,A);
    y = zeros(size(x,1))
    mL      = maximum(G.L)
    lev = map(level,G[:,:]);nc = size(G.covers,1)
    T0 = zeros(UInt32,length(G))
    T1 = zeros(UInt32,length(G))
    J         = ones(Int,mL+1,D)
    B         = ones(mL+1,D)

    topc = [findlast((G.covers[:,d].==maximum(G.covers[:,d]))) for d=1:D]
    cTs = cumsum(prod(Ts,2))

    id = zeros(Int,D)

    for i = 1:size(x,1)
        for l = 2:mL+1
            dm = (l<3) ? 2^(l-1) : 2^(l-2)
            m = 2^(l-1)+1
            for d = 1:D
                j 	= clamp(round(Int,x[i,d]*dm+1/2),1,dm)
                xij =  l==2 ? Float64(j-1) : (2j-1)/(m-1.0)
                dx = (1.0-(m-1.0)*abs(x[i,d]-xij))
                B[l,d] 	= (dx>0.0) ? dx : 0.0
                J[l,d]  = j
            end
        end
        cnt=0
        for ii in topc
            cnt+=1
            b  = B[G.covers[ii,D],D]*B[G.covers[ii,1],1]
            id1 = J[G.covers[ii,D],D]-1
            for d = D-1:-1:2
                b*=B[G.covers[ii,d],d]
                id1 = id1*G.covers_dM[ii,d]+(J[G.covers[ii,d],d]-1)
            end
            id[cnt] = (J[G.covers[ii,1],1]-1)+G.covers_dM[ii,1]*id1+1
        end
        t1 =         (id[1]-1)*Ts[1]+1 : (id[1])*Ts[1]
        t2 = cTs[1]+((id[2]-1)*Ts[2]+1 : (id[2])*Ts[2])
        nT = c_intersect3(sub(Tp,t1),
                          sub(Tp,t2),
                          T0,length(t1),length(t2))
        Todd = false
        for d = 2:D-1
            dr = cTs[d]+((id[d+1]-1)*Ts[d+1]+1:(id[d+1])*Ts[d+1])
            if !Todd
                nT = c_intersect3(sub(Tp,dr),T0,T1,length(dr),Int(nT))
            else
                nT = c_intersect3(sub(Tp,dr),T1,T0,length(dr),Int(nT))
            end
            Todd= !Todd
        end
        T3 = Todd ? T1 : T0
        for ii in 1:nc
            b = 1.0
            for d = 1:D
                b*=B[G.covers[ii,d],d]
            end
            y[i]+=b*w[T3[ii]]
        end
    end
    y
end









G = NGrid([5,5,5,5])
T,Ts,Tp=getparents1(G);

A = sum(G[:,:],2)[:]
x = rand(10000,length(G.L))
yc = cinterp(G,Tp,Ts,A,x)
# y = G(A,x);
# y0 = SparseGrids.jl_interpbig(x,G,A);
# y1 = interp1(G,T,A,x);
# y2 = interp2(G,Tp,Ts,A,x);
# y3 = interp3(G,Tp,Ts,A,x);
# @assert maximum(abs(y-y0))<1e-14
# @assert maximum(abs(y-y1))<1e-14
# @assert maximum(abs(y-y2))<1e-14
# @assert maximum(abs(y-y3))<1e-14

x
# @time y = G(A,x);
# @time y0 = SparseGrids.jl_interpbig(x,G,A);
# @time y1 = interp1(G,T,A,x);
# @time y2 = interp2(G,Tp,Ts,A,x);
# @time y3 = interp3(G,Tp,Ts,A,x);
# @time yc = cinterp(G,Tp,Ts,A,x);


# b = @benchmark G(A,x);
# b0 =@benchmark SparseGrids.jl_interpbig(x,G,A);
# b1 =@benchmark interp1(G,T,A,x);
# b2 = @benchmark interp2(G,Tp,Ts,A,x);
# b3 = @benchmark interp3(G,Tp,Ts,A,x);
# bc = @benchmark cinterp(G,Tp,Ts,A,x)

# Profile.clear_malloc_data()
# interp1(G,uT,A,x);
# interp2(G,uT,A,x);
nothing

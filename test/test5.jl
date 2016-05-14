using Benchmarks,DevTools,MyDisplay
res=Array(Benchmarks.Results,5,2)

# for d = 2:6
#     G = NGrid(ones(Int,d)*6)
#     A = rand(length(G))
#     x = rand(length(G)*10,length(G.L))
#
#     G(A,x)
#     SparseGrids.jl_interp(G,A,x)
#     res[d-1,1] = @benchmark G(A,x)
#     sleep(15)
#     res[d-1,2] = @benchmark SparseGrids.jl_interp(G,A,x)
#     sleep(15)
# end

import SparseGrids:nXtoU,@threads,subs!
for b in [(Linear,Lbj),(Quadratic,Qbj)]
    for D = 2:12
        coverloop = :(@inbounds for ii in nc
            b   = B[G.covers[ii,$D],$D]*B[G.covers[ii,1],1]
            id1 = J[G.covers[ii,$D],$D]
        end)
        for d in D-1:-1:2
            push!(coverloop.args[2].args[2].args,:(b  *= B[G.covers[ii,$d],$d]))
            push!(coverloop.args[2].args[2].args,:(id1 = id1*G.covers_dM[ii,$d]+J[G.covers[ii,$d],$d]))
        end

        push!(coverloop.args[2].args[2].args,:(id1 = id1*G.covers_dM[ii,1]+G.covers_loc[ii]+J[G.covers[ii,1],1];yi += b*w[id1]))

        f=:(function jl_interp(G::NGrid{$D,$(b[1])},A::Vector{Float64},xi::Array{Float64,2})
            w         = getW(G,A)
            x         = nXtoU(xi,G.bounds)
            nx        = size(x,1)
            y         = zeros(Float64,nx)
            nc         = 1:size(G.covers,1)
            dr1 = 1:$D
            mL      = maximum(G.L)+1
            J         = zeros(Int,mL,$D)
            B         = ones(mL,$D)

            @threads [J,B] for i = 1:nx
                $(b[2])
                yi = 0.0
                coverloop
                y[i] = yi
            end
            y
        end)

        subs!(f,:coverloop=>coverloop)
        eval(current_module(),f)
    end
end

for b in [(Linear,Lbj),(Quadratic,Qbj)]
    for D = 2:12
        coverloop = :(@inbounds for ii in nc
            b   = B[G.covers[ii,$D],$D]*B[G.covers[ii,1],1]
            id1 = J[G.covers[ii,$D],$D]
        end)
        for d in D-1:-1:2
            push!(coverloop.args[2].args[2].args,:(b  *= B[G.covers[ii,$d],$d]))
            push!(coverloop.args[2].args[2].args,:(id1 = id1*G.covers_dM[ii,$d]+J[G.covers[ii,$d],$d]))
        end

        push!(coverloop.args[2].args[2].args,:(id1 = id1*G.covers_dM[ii,1]+G.covers_loc[ii]+J[G.covers[ii,1],1];for d=1:nA;y[i,d] += b*w[id1,d];end))

        f=:(function jl_interp(G::NGrid{$D,$(b[1])},A::Array{Float64,2},xi::Array{Float64,2})
            w         = getW(G,A)
            x         = nXtoU(xi,G.bounds)
            nx        = size(x,1)
            nA        = size(A,2)
            y         = zeros(Float64,nx,nA)
            nc         = 1:size(G.covers,1)
            dr1 = 1:$D
            mL      = maximum(G.L)+1
            J         = zeros(Int,mL,$D)
            B         = ones(mL,$D)

            for i = 1:nx
                $(b[2])
                coverloop
            end
            y
        end)

        subs!(f,:coverloop=>coverloop)
        eval(current_module(),f)
    end
end



G = NGrid([6,6,3,4,4])
A = rand(length(G),3)
x = rand(length(G)*10,length(G.L))
# bc = @benchmark G(A[:,1],x)
# bj = @benchmark SparseGrids.jl_interp(G,A[:,1],x)
b0=(@benchmark SparseGrids.jl_interp(G,A,x))
b1= (@benchmark jl_interp1(G,A,x))
 b2=(@benchmark jl_interp2(G,A,x))

Profile.init(delay=0.00000001)
Profile.clear()
@profile jl_interp1(G,A,x)
profiler()
Profile.clear()
@profile jl_interp2(G,A,x)
profiler()


#
#
# b0 =@benchmark G(A,x)
# sleep(20)
# # b1 =@benchmark jl_interp(G,A,x)
# # sleep(20)
# b5 =@benchmark jl_interp5(G,A,x)
# sleep(20)
# b6 =@benchmark jl_interp6(G,A,x)
# sleep(20)
# b7 =@benchmark SparseGrids.jl_interp(G,A,x)
#
# @time for i = 1:50; G(A,x) ; end
# @time for i = 1:50; SparseGrids.jl_interp(G,A,x) ; end
# plot([b0,b7])
# #
# Y = G(A,x)
# @assert SparseGrids.jl_interp(G,A,x)==Y
#
# # sleep(10)
# # b0 = @benchmark G(A,x)
# # sleep(20)
# # b1 = @benchmark SparseGrids.jl_interp(G,A,x)
# # plot([b0,b1])
# plot([(@benchmark jl_interp(G,A,x)),(@benchmark jl_interpfm(G,A,x))])
#
# Profile.clear()
# Profile.init(delay=0.000001)
# @profile jl_interp(G,A,x)
#
# s=pview();

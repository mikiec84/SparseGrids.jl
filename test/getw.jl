
function getW2(G::SparseGrids.CurtisClenshaw.Grid,A)
	Aold = zeros(G.n)
	dA = zeros(G.n)
	w = zeros(G.n)
	ccall((:_Z8w_get_l2PdiiPiS0_S_iS_S_S_, SparseGrids.CurtisClenshaw.libsparse),
		Void,
		(Ptr{Float64},Int32,Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64},Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		pointer(G.grid),G.n,G.d,pointer(G.lvl_s),pointer(G.lvl_l),pointer(A),G.q,pointer(Aold),pointer(dA),pointer(w))
	return w
end

G = CurtisClenshaw.Grid(4,4)
A =rand(G.n)
wi = getWinv(G)


coverings = unique(G.index,1)
[sum(abs(wi[G.level.==level,:]).>0,2) for level = 0:G.q]

[wi[all(G.index.==coverings[i,:]',2),:] for i = 1:70]

xi = rand(100,G.d)
q = getQ(xi,G)

sum(abs(q).>0,2)

x= rand(50)
f(x)= (x[:] = 0.0;nothing)
f1(x) = (fill!(x,0.0);nothing)




function getWinv2(G)
	drange 	= collect(1:G.d)
	qA		= speye(G.n)
	qAold 	= spzeros(G.n,G.n)
	qw		= spzeros(G.n,G.n)
	qtemp 	= spzeros(G.n)

	for i = 1:G.lvl_l[1]-1
		for j = 1:G.n
		    qw[i,j] = qA[i,j] - qAold[i,j]
		end
    end
    for i in 1:G.n
    	fill!(qtemp,0)
        for ii in 1:G.lvl_l[1]-1
            temp2=1.0
            for d in drange
                temp2 *= basis_func(G.grid[i,d],G.grid[ii,d],G.lvl_s[ii,d])
            end
            for j = 1:G.n
	            qtemp[j] += temp2*qw[ii,j]
	        end
        end

	    for j = 1:G.n
		    qAold[i,j] += qtemp[j]
		end
    end

	for q in 1:G.q
	    for i = G.lvl_l[q]:G.lvl_l[q+1]-1
	    	for j = 1:G.n
		        qw[i,j] = qA[i,j] - qAold[i,j]
		    end
	    end
	    for i in G.lvl_l[q]:G.n
	        fill!(qtemp,0)
	        for ii in G.lvl_l[q]:G.lvl_l[q+1]-1
	            temp2=1.0
	            for d in drange
	                @inbounds temp2 *= basis_func(G.grid[i,d],G.grid[ii,d],G.lvl_s[ii,d])
	            end
	            @simd for j = 1:G.n
		            @inbounds qtemp[j] += temp2*qw[ii,j]
		        end
	        end
	        @simd for j = 1:G.n
			    @inbounds qAold[i,j] += qtemp[j]
			end
	    end
	end
   return sparse(qw)
end

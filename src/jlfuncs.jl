 for i = 1:length(GridSpecs)
    eval(:(function buildW(G::NGrid{$(GridSpecs[i][2])})
        nc 		= size(G.coverings,1)
        mL      = maximum(G.L)
        D       = length(G.L)
        J 		= ones(Int,mL+1,D)
        B 		= zeros(mL+1,D)
        for i = 1:length(G)
            for d = 1:D
                for l = 1:G.level[i]+1
                    j 	= clamp(round(Int,G.grid[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
                    B[l,d] 	= $(GridSpecs[i][3])(G.grid[i,d],cc_dg(l,j),Int16(cc_M(l)))
                    J[l,d]  = j
                end
            end
            for ii = 1:G.coveringsid[i]
                b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
                id1 = J[G.coverings[ii,D],D]-1
                for d = D-1:-1:2
                    b*=B[G.coverings[ii,d],d]
                    id1 = id1*G.coverings_dM[ii,d]+(J[G.coverings[ii,d],d]-1)
                end
                id1=(J[G.coverings[ii,1],1]-1)+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
                if b>0
                    push!(G.IDs[i],id1)
                    push!(G.Bs[i],id1==i ? b : -b)
                end
            end
        end
        return
    end))



	eval(:(function jl_interpbig(xi::Array{Float64},G::NGrid{$(GridSpecs[i][2])},A::Vector{Float64})
		w 		= getW(G,A)
		x 		= nXtoU(xi,G.bounds)
		y 		= zeros(size(x,1))
		nc 		= size(G.coverings,1)
	    mL      = maximum(G.L)
	    D       = length(G.L)

		for i = 1:size(x,1)
			J 		= ones(Int,mL+1,D)
			B 		= zeros(mL+1,D)
			id  	= zeros(Int,D)
			for d = 1:D
				for l = 1:mL+1
					j 	= clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
					B[l,d] 	= $(GridSpecs[i][3])(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
					J[l,d]  = j
				end
			end

            for ii = 1:nc
                b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
                id1 = J[G.coverings[ii,D],D]-1
                for d = D-1:-1:2
                    b*=B[G.coverings[ii,d],d]
                    id1 = id1*G.coverings_dM[ii,d]+(J[G.coverings[ii,d],d]-1)
                end
                id1=(J[G.coverings[ii,1],1]-1)+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
                y[i]+=b*w[id1]
            end
		end
		y
	end))

    eval(:(function jl_interpbig(xi::Array{Float64},G::NGrid{$(GridSpecs[i][2])},A::Array{Float64,2})
        x 		= nXtoU(xi,G.bounds)
        nc 		= size(G.coverings,1)
        nA      = size(A,2)
		w 		= getW(G,A)
		y 		= zeros(size(x,1),nA)
	    mL      = maximum(G.L)
	    D       = length(G.L)

		for i = 1:size(x,1)
			J 		= ones(Int,mL+1,D)
			B 		= zeros(mL+1,D)
			id  	= zeros(Int,D)
            for d = 1:D
				for l = 1:mL+1
					j 	= clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
					B[l,d] 	= $(GridSpecs[i][3])(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
					J[l,d]  = j
				end
			end

            for ii = 1:nc
                b  = B[G.coverings[ii,D],D]*B[G.coverings[ii,1],1]
                id1 = J[G.coverings[ii,D],D]-1
                for d = D-1:-1:2
                    b*=B[G.coverings[ii,d],d]
                    id1 = id1*G.coverings_dM[ii,d]+(J[G.coverings[ii,d],d]-1)
                end
                id1=(J[G.coverings[ii,1],1]-1)+G.coverings_dM[ii,1]*id1+1+G.coveringsloc[1][ii]-1
                for iii = 1:nA
                    y[i,iii]+=b*w[id1,iii]
                end
            end
		end
		y
	end))
end


function jl_getWbig(G::NGrid,A::Vector{Float64})
    w = copy(A)
    for i = 1:length(G)
        for l ∈ 1:length(G.Bs[i])-1
            w[i]+=  G.Bs[i][l]*w[G.IDs[i][l]]
        end
    end
    return w
end

function jl_getWbig(G::NGrid,A::Array{Float64,2})
    nA = size(A,2)
    w = copy(A)
    for i = 1:length(G)
        for l ∈ 1:length(G.Bs[i])-1
            for a = 1:nA
                w[i,a]+=  G.Bs[i][l]*w[G.IDs[i][l],a]
            end
        end
    end
    return w
end

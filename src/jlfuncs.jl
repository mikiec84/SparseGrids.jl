for i = 1:length(GridSpecs)
    eval(:(function jl_getW(G::NGrid{$(GridSpecs[i][1]),$(GridSpecs[i][2])},A::Vector{Float64})
        Aold = zeros(size(G.grid,1))
        w = zeros(size(G.grid,1))
        N = size(G.grid,1)
        D = size(G.grid,2)
        for l = 0:maximum(G.L)
            for i=G.level_loc[l+1]:G.level_loc[l+2]-1
                w[i] = A[i]-Aold[i]
            end
            for i = G.level_loc[l+2]:N
                for ii=G.level_loc[l+1]:G.level_loc[l+2]-1
                    temp2=1.0
                    for d=1:D
                        temp2*=$(GridSpecs[i][3])(G.grid[i,d],G.grid[ii,d],G.level_M[ii,d])
                        temp2==0 && break
                    end
                    Aold[i] += temp2*w[ii]
                end
            end
        end
        return w
    end))

    eval(:(function jl_interpslow(xi::Array{Float64},G::NGrid{$(GridSpecs[i][1]),$(GridSpecs[i][2])},A::Vector{Float64})
        w = getW(G,A)
        x = nXtoU(xi,G.bounds)
        nx = size(x,1)
        y = zeros(nx)
        N = size(G.grid,1)
        D = size(G.grid,2)
        for i = 1:nx
            for ii = 1:N
                temp2 = 1.0
                for d = 1:D
                    temp2 *= $(GridSpecs[i][3])(x[i,d],G.grid[ii,d],G.level_M[ii,d])
                    temp2==0.0 ? break : nothing
                end
                y[i]+= temp2*w[ii]
            end
        end
        return y
    end))

    eval(:(function jl_interp(xi::Array{Float64},G::NGrid{$(GridSpecs[i][1]),$(GridSpecs[i][2])},A::Vector{Float64})
        w = getW(G,A)
        x = nXtoU(xi,G.bounds)
        nx = size(x,1)
        y = zeros(nx)
        N = size(G.grid,1)
        D = size(G.grid,2)
        for i = 1:nx
            ii = Int32(1)
            while ii<=N
                temp2 = 1.0
                for d = 1:D
                    temp2 *= $(GridSpecs[i][3])(x[i,d],G.grid[ii,d],G.level_M[ii,d])
                    temp2==0.0 ? break : nothing
                end
                y[i]+= temp2*w[ii]
                temp2==0.0 ? ii+=Int32(1) : ii=G.nextid[ii]
            end
        end
        return y
    end))

	eval(:(function jl_interpbig(xi::Array{Float64},G::NGrid{$(GridSpecs[i][1]),$(GridSpecs[i][2])},A::Vector{Float64})
		w 		= getW(G,A)
		x 		= nXtoU(xi,G.bounds)
		y 		= zeros(size(x,1))
		nc 		= size(G.coverings,1)
	    mL      = maximum(G.L)
	    D       = length(G.L)

		for i = 1:length(y)
			J 		= ones(Int,mL+1,D)
			B 		= zeros(mL+1,D)
			id  	= zeros(Int,D)
			for d = 1:D
				for l = 1:mL+1
					j 	= clamp(round(Int,x[i,d]*($(GridSpecs[i][5])(l))+1/2),1,$(GridSpecs[i][5])(l))
					B[l,d] 	= $(GridSpecs[i][3])(x[i,d],$(GridSpecs[i][7])(l,j),Int16($(GridSpecs[i][4])(l)))
					J[l,d]  = $(GridSpecs[i][6])(l,j,mL+1)
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
	end))

end

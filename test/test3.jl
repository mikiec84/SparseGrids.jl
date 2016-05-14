function getB(x)
    i=1
    D = length(x)
    mL = 5
    J         = ones(Int,mL+1,D)
    B         = ones(mL+1,D)
    for d = 1:D
        for l = 1:mL+1
            j 	= clamp(round(Int,x[i,d]*(cc_dM(l))+1/2),1,cc_dM(l))
            B[l,d] 	= bf(x[i,d],cc_dg(l,j),Int16(cc_M(l)))
            J[l,d]  = j
        end
    end
    B,J
end

function getB1(x)
    i=1
    D = length(x)
    mL = 5
    J         = ones(Int,mL+1,D)
    B         = ones(mL+1,D)
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
    B,J
end

x = rand(1,5)

B,J = getB(x)
B1,J1 = getB1(x)
@time  for i = 1:10000; getB(x);  end
@time  for i = 1:10000; getB1(x); end

@assert B==B1 && J==J1

15.33/13.31

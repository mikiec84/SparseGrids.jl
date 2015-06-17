
function XtoX(X::Array{Float64,1},bnds1::Array{Float64,1},bnds2::Array{Float64,1})
	lb1 = minimum(bnds1)
	ub1 = maximum(bnds1)
	lb2 = minimum(bnds2)
	ub2 = maximum(bnds2)
	out = Array(Float64,length(X))
	for i = 1:length(X)
		out[i] = (X[i]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
	end
	return out
end

function nXtoX(X::Array{Float64,2},bnds1::Array{Float64,2},bnds2::Array{Float64,2})
	N,R = size(X)
	out = Array(Float64,size(X))
	for r = 1:R
		lb1 = minimum(bnds1[:,r])
		ub1 = maximum(bnds1[:,r])
		lb2 = minimum(bnds2[:,r])
		ub2 = maximum(bnds2[:,r])
		for i = 1:N
			out[i,r] = (X[i,r]-lb1)/(ub1-lb1)*(ub2-lb2)+lb2
		end
	end
	return out
end
nXtoX(X::Array{Float64,1},bnds1::Array{Float64,2},bnds2::Array{Float64,2})=XtoX(X,bnds1[:],bnds2[:])

UtoX(U::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(U,[0.0,1.0],bnds)
XtoU(X::Array{Float64,1},bnds::Array{Float64,1}) = XtoX(X,bnds,[0.0,1.0])

nUtoX(U::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(U,[zeros(1,size(U,2));ones(1,size(U,2))],bnds)
nXtoU(X::Array{Float64,2},bnds::Array{Float64,2}) = nXtoX(X,bnds,[zeros(1,size(X,2));ones(1,size(X,2))])
nXtoU(X::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(X,bnds[:],[0.0,1.0])
nUtoX(U::Array{Float64,1},bnds::Array{Float64,2}) = XtoX(U,[0.0,1.0],bnds[:])

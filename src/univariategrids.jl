cc_M(i::Int)            = (i==1) ? 1 : 2^(i-1)+1
cc_iM(i::Int)           = (i==1) ? 1 : Int(log2(i-1)+1)
cc_dM(i::Int)           = (i<3) ? 2^(i-1) : 2^(i-2)
cc_g(i::Int,j::Int)     = (i==1) ?  0.5  : (j-1)/(cc_M(i)-1.0)
cc_g(i::Int)            = Float64[cc_g(i,j) for j = 1:cc_M(i)]
cc_dg(i::Int,j::Int)    = i==2 ? [0.0,1.0][j] : cc_g(i,2j)
cc_dg(i::Int)           = Float64[cc_dg(i,j) for j = 1:cc_dM(i)]
function cc_itoj(i,j,q)
    if i ==1
        return div(cc_M(q)-1,2)+1
    elseif i==2
        return j==1 ? 1 : cc_M(q)
    else
        return div((cc_M(q)-1)*(j*2-1),cc_M(i)-1)+1
    end
end
function cc_bf_l(x::Float64,xij::Float64,mi::Int16)
    if (mi==1)
        return 1.0
    end
    dx = (1.0-(mi-1.0)*abs(x-xij))
    return (dx>0.0) ? dx : 0.0
end
function cc_bf_q(x::Float64,xij::Float64,mi::Int16)
	if (mi==1)
		return 1.0
	end
	dx = 1.0-((mi-1.0)*(x-xij))^2
	return (dx>0.0) ?  dx : 0.0
end


m_M(i::Int) =  2^i+1
m_iM(i::Int) = Int(log2(i-1))
m_dM(i::Int) = i==1 ? 3 : 2^(i-1)
m_g(i::Int,j::Int)    = (j-1)/(m_M(i)-1)
m_g(i::Int)           = Float64[m_g(i,j) for j = 1:m_M(i)]
m_dg(i::Int,j::Int) = i==1 ? m_g(i)[j] : m_g(i,2j)
m_dg(i::Int)           = Float64[m_dg(i,j) for j = 1:m_dM(i)]
m_itoj(i,j,q) =	div((j-1)*(m_M(q)-1),m_M(i)-1)+1
function m_bf_l(x::Float64, xij::Float64, mi::Int16)
    if (abs(x-xij)<1/(mi-1))
        return 1.0-(mi-1)*abs(x-xij)
    else
        return 0.0
    end
end
function m_bf_q(x::Float64, xij::Float64, mi::Int16)
	dx = 1.0-((mi-1.0)*(x-xij))^2
	return dx>0.0 ? dx : 0.0
end

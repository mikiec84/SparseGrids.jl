abstract GridType

type UnivariateGrid{T<:GridType}
    M::Function
    iM::Function
    dM::Function
    g::Function
    dg::Function
	itoi::Function
end


type CCGrid <:GridType end

cc_M{T}(i::T)            = (i==1) ? T(1) : T(2)^(i-T(1))+T(1)
cc_iM(i::Int)           = (i==1) ? 1 : Int(log2(i-1)+1)
cc_dM(i::Int)           = (i<3) ? 2^(i-1) : 2^(i-2)
cc_g(i::Int,j::Int)     = (i==1) ?  0.5  : (j-1)/(cc_M(i)-1.0)
cc_g(i::Int)            = Float64[cc_g(i,j) for j = 1:cc_M(i)]
cc_dg(i::Int,j::Int)    = i==2 ? [0.0,1.0][j] : cc_g(i,2j)
cc_dg(i::Int)           = Float64[cc_dg(i,j) for j = 1:cc_dM(i)]
function cc_itoi(i1,j,i2)
    if i1 ==1
        return div(cc_M(i2)-1,2)+1
    elseif i1==2
        return j==1 ? 1 : cc_M(i2)
    else
        return div((cc_M(i2)-1)*(j*2-1),cc_M(i1)-1)+1
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

function cc_simpsonsw(i::Int,j::Int)
    @assert 1≤j≤CC.M(i)
    if i==1
        return 1.0
    end
    w = 1/(CC.M(i)-1)/3
    if j==1 || j==CC.M(i)
        return w
    else
        w*= mod(j,2)==0 ? 4 : 2
    end
    return w
end

cc_simpsonsw(i::Int) = Float64[cc_simpsonsw(i,j) for j = 1:CC.M(i)]

function cc_dsimpsonsw(i::Int)
    if i==1
        return cc_simpsonsw(i)
    elseif i==2
        return cc_simpsonsw(i)-[0;1.0;0]
    else
        return cc_simpsonsw(i)-vec(Float64[cc_simpsonsw(i-1)';zeros(1,CC.M(i-1))])[1:2CC.M(i-1)-1]
    end
end

const CC = UnivariateGrid{CCGrid}(cc_M,cc_iM,cc_dM,cc_g,cc_dg,cc_itoi)



type MaxGrid <:GridType end
m_M(i::Int) =  2^i+1
m_iM(i::Int) = Int(log2(i-1))
m_dM(i::Int) = i==1 ? 3 : 2^(i-1)
m_g(i::Int,j::Int)    = (j-1)/(m_M(i)-1)
m_g(i::Int)           = Float64[m_g(i,j) for j = 1:m_M(i)]
m_dg(i::Int,j::Int) = i==1 ? m_g(i)[j] : m_g(i,2j)
m_dg(i::Int)           = Float64[m_dg(i,j) for j = 1:m_dM(i)]
m_itoi(i,j,q) =	div((j-1)*(m_M(q)-1),m_M(i)-1)+1
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

const Max = UnivariateGrid{MaxGrid}(m_M,m_iM,m_dM,m_g,m_dg,m_itoi)






type CC2Grid <: GridType end


function cc2_M(i::Int)
    if i==1
        return 2
    elseif i==2
        return 3
    else
        return 2^(i-1)+1
    end
end

function cc2_dM(i::Int)
    if i==1
        return 2
    else
        return 2^(i-2)
    end
end

function cc2_iM(i::Int)
    if i==2
        return 1
    elseif i== 3
        return 2
    else
        return Int(log2(i-1)+1)
    end
end

cc2_g(i::Int,j::Int)    = (j-1)/(cc2_M(i)-1.0)
cc2_g(i::Int)            = Float64[cc2_g(i,j) for j = 1:cc2_M(i)]

function cc2_dg(i::Int,j::Int)
    if i==1
        return [0.0,1.0][j]
    else
        return cc2_g(i,2j)
    end
end
cc2_dg(i::Int)           = Float64[cc2_dg(i,j) for j = 1:cc2_dM(i)]

function cc2_itoi(i,j,q)
    if i ==1
        return div(cc2_M(q)-1,2)+1
    elseif i==2
        return j==1 ? 1 : cc2_M(q)
    else
        return div((cc2_M(q)-1)*(j*2-1),cc2_M(i)-1)+1
    end
end


const CC2 = UnivariateGrid{CC2Grid}(cc2_M,cc2_iM,cc2_dM,cc2_g,cc2_dg,cc2_itoi)

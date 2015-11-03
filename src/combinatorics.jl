
function swap!(x::Vector{Int},i::Int,j::Int)
	ta = x[i]
	x[i]=x[j]
	x[j]=ta
	return
end

function next_permute!(seq::Vector{Int})
	n = length(seq)
	j = n-1
	while seq[j]>=seq[j+1]
		j-=1
		if j==0
			return true
		end
	end
	l = n
	while seq[j]>=seq[l]
		l-=1
	end
	swap!(seq,j,l)

	seq[j+1:n]=reverse(seq[j+1:n])
	return false
end

function get_perms(seq::Vector{Int})
	a = sort(seq)
	out = Array(Vector{Int},1)
	out[1] = [a;]
	if length(seq)==1
		return out
	end
	done = false
	while !done
		done=next_permute!(a)
		push!(out,[a;])
	end
	if out[end]==out[end-1]
		out = out[1:end-1]
	end
	return out::Vector{Vector{Int}}
end

function comb(D::Int,Q::Int)
	if D ==1
		return [Index([Q])]
	end
	Pq = partitions(Q)
	out = Array(Vector{Int},0)

	for p in Pq
		if length(p)==D
			pPq=get_perms(p)
			for j = 1:length(pPq)
				push!(out,pPq[j])
			end
		end
	end

	return [Index(x) for x in out]
end

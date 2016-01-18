function test{T<:UnivariateGrid}(ug::T)
    @assert ug.M(1) == ug.dM(1)
    for i = 2:10
        @assert ug.M(i)-ug.M(i-1)==ug.dM(i)
    end
    @assert ug.g(1) == ug.dg(1)
    for i = 2:10
        @assert setdiff(Set(ug.g(i)),Set(ug.g(i-1)))==Set(ug.dg(i))
    end
end

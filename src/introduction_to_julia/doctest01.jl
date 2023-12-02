"""
    bps(a::Number) 

어떤 기능을 하는지 설명한다. 여기서는 행렬에 2를 곱한다.


# Example

```julia-repl
julia> bps([1 2;3 4])

2×2 Matrix{Int64}:
 2  4
 6  8
```

see also [`sin`](@ref).

"""
function bps(a::Number) 
    return a*2
end
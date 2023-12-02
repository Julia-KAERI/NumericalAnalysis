using LinearAlgebra

export SimplePolynomial, derivative, integrate, monic, definite_integrate

"""
    SimplePolynomial{T}

간단한 polynomial 자료형.
"""
struct SimplePolynomial{T}
    coeffs :: Vector{T}
    
    function SimplePolynomial(a::AbstractVector{P}) where P <: Number
        if length(a) == 0 
            return new{P}(zeros(T, 1))
        else 
            last_nz = findlast(!iszero, a)
            a_last = max(1, isnothing(last_nz) ? 0 : last_nz)
            return new{P}(a[1:a_last])
        end
    end

    function SimplePolynomial{T}(a::AbstractVector{P}) where {T <: Number, P<:Number}
        if length(a) == 0 
            return new{T}(zeros(T, 1))
        else 
            last_nz = findlast(!iszero, a)
            a_last = max(1, isnothing(last_nz) ? 0 : last_nz)
            return new{T}(a[1:a_last])
        end
    end

end

function (p::SimplePolynomial)(x::T) where T <: Number
    return evalpoly(x, p.coeffs)
end


function (p::SimplePolynomial)(x::T) where T <: Matrix{N} where N<:Number
    r = UniformScaling(p.coeffs[1])
    @assert size(x)[1] == size(x)[2] # 정사각 행렬에 대해서만 가능하다.
    for i in 2:length(p.coeffs)
        @inbounds r +=  p.coeffs[i]*x^(i-1)
    end
    return r
end

function Base.show(io::IO, p::SimplePolynomial{T}) where T<:Number
    result = ""
    n = length(p)
    if n == 1 && iszero(p.coeffs[1])
        result = "0"
    else 
        for (i, v) in enumerate(p.coeffs[end:-1:1])
            vp = string(abs(v))
            if v > zero(T) && i>1
                result *= " + " * vp * " x^$(n-i)"
            elseif v< zero(T) && i > 1
                result *= " - " * vp * " x^$(n-i)"
            elseif v > zero(T)
                result *= vp * " x^$(n-i)"
            elseif v<zero(T)
                result *= " -"*vp * " x^$(n-i)"
            end

            if i == n
                if !iszero(v)
                    result = result[1:end-4]
                else 
                    result = result
                end
            end
        end
    end
    println(io, "SimplePolynomial($(result[1:end]))")
end


Base.length(p::SimplePolynomial) = length(p.coeffs)

order(p::SimplePolynomial) = length(p)-1
degree(p::SimplePolynomial) = order(p)


function Base.zero(a::P) where P<:SimplePolynomial
    return SimplePolynomial([zero(eltype(a.coeffs)), ])
end

function Base.one(a::P) where P<:SimplePolynomial
    return SimplePolynomial([one(eltype(a.coeffs)), ])
end

function Base.:-(b::P) where {P<: SimplePolynomial}
    coeffs = -b.coeffs
    return SimplePolynomial(coeffs)
end

function Base.:+(a::T, b::SimplePolynomial{P}) where {T <: Number, P <: Number} 
    rT = promote_type(T, P)
    coeffs = rT.(b.coeffs)
    coeffs[1] += a
    return SimplePolynomial(coeffs)
end

function Base.:+(b::SimplePolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return a+b
end

function Base.:+(a::SimplePolynomial{P1}, b::SimplePolynomial{P2}) where {P1 <: Number, P2 <: Number} 
    rT = promote_type(P1, P2)
    if length(b) > length(a)
        coeffs = zeros(rT, length(b))
        coeffs[1:length(a)] = a.coeffs[:]
        coeffs += b.coeffs
    else 
        coeffs = zeros(rT, length(a))
        coeffs[1:length(b)] = b.coeffs[:]
        coeffs += a.coeffs
    end
    return SimplePolynomial(coeffs)
end

function Base.:-(a::SimplePolynomial{P1}, b::SimplePolynomial{P2}) where {P1 <: Number, P2 <: Number} 
    return a + (-b)
end

function Base.:-(b::SimplePolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return b+(-a)
end

function Base.:-(a::T, b::SimplePolynomial{P}) where {T <: Number, P <: Number} 
    return a+(-b)
end

function Base.:*(a::T, b::SimplePolynomial{P}) where {T <: Number, P <: Number} 
    return SimplePolynomial(b.coeffs*a)
end

function Base.:*(b::SimplePolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return SimplePolynomial(b.coeffs*a)
end

function Base.:*(a::SimplePolynomial{P1}, b::SimplePolynomial{P2}) where {P1 <: Number, P2 <:Number} 
    rT = promote_type(P1, P2)
    ord1, ord2 = degree(a), degree(b)
    ord = ord1*ord2
    coef = zeros(rT, ord+2)
    
    for i in 0:ord1, j in 0:ord2
        @inbounds coef[i+j+1] += a.coeffs[i+1]*b.coeffs[j+1]
    end
    return SimplePolynomial(coef)
end



function Base.:/(b::SimplePolynomial{P}, a::T) where {P <: Number, T <: Number} 
    return SimplePolynomial(b.coeffs/a)
end

"""
    monic(p::P) where {P<:SimplePolynomial}

return monic polynomial of which highest coefficient is 1
"""
function monic(p::P) where P<:SimplePolynomial
    return p/p.coeffs[end]
end

function derivative(p::SimplePolynomial)
    if length(p) < 2 
        return SimplePolynomial([one(eltype(p.coeffs)), ])
    else
        coeffs = p.coeffs[2:end] .* (1:(length(p)-1))
        return SimplePolynomial(coeffs)
    end
end


"""
    integrate(p::P, C::N) where {P<:SimplePolynomial, N<:Number}

integrate polynomial p with constant C.
"""
function integrate(p::SimplePolynomial, a::Union{Nothing, Number}=nothing, b::Union{Nothing, Number}=nothing) 

    if eltype(p.coeffs) <: Integer
        coeffs = zeros(Float64, length(p)+1)
    else 
        coeffs = zeros(eltype(p.coeffs), length(p)+1)
    end
    
    for i in 1:length(p.coeffs)
        coeffs[i+1] = p.coeffs[i]/(i)
    end
    
    
    if a === nothing && b === nothing # 상수항이 0 인 부정적분
        coeffs[1] = zero(eltype(coeffs))
        return SimplePolynomial(coeffs)
    elseif a === nothing || b === nothing # 상수항이 a 혹은 b 로 주어진 부정적분
        coeffs[1] = a
        return SimplePolynomial(coeffs)
    else # a 에서 b 구간 까지의 정적분
        return evalpoly(b, coeffs) - evalpoly(a, coeffs)
    end
end


"""
    polynomial_from_roots(xp::Vector{T}) where T<:Number

get monic polynomial having roots of xp[1],..., xp[end], i.e (x-xp[1]) cdots (x-xp[end])
"""
function polynomial_from_roots(xp::AbstractVector{T}) where T<:Number 
    return prod([SimplePolynomial([-x0, 1]) for x0 in xp])
end


function valdermond_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    
    N = length(xp)
    @assert length(xp) == length(yp)
    V = [x^(j-1) for x in xp, j in 1:length(xp)]
    return SimplePolynomial(V\yp)
end

function lagrange_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}

    N = length(xp)
    @assert length(xp) == length(yp)
    
    r = SimplePolynomial([zero(T2), ])
    for i in 1:N
        coef = yp[i]
        rt = one(T2)
        for j in 1:N
            if i ≠ j
                @inbounds coef = coef/(xp[i]-xp[j])
                @inbounds rt = rt*SimplePolynomial([-xp[j], 1.0])
            end
        end
        r += rt*coef
    end
    return r
end

function newton_polynomial(
    xp::AbstractVector{T1}, 
    yp::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    n = length(xp)    
    @assert n == length(yp)
    T = promote_type(T1, T2)
    N = LowerTriangular(ones(T, n, n))
    for j in 2:n, i in j:n
        @inbounds N[i, j] = N[i, j-1]*(xp[i] - xp[j-1]) 
    end
    a = N\yp
    r = SimplePolynomial([a[1], ])
    for i in 2:(n)
        @inbounds r += a[i] * polynomial_from_roots(xp[1:i-1])
    end
    return r

end

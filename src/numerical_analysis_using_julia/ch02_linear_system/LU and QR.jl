using LinearAlgebra, BenchmarkTools 

"""
    II(n)

n \times n 단위행렬을 리턴한다.
"""
function II(n, T::Type=Float64) 
    return Matrix{T}(I, n, n)
end

"""
    Ls(A, b)

하삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Ls(L::Matrix{T}, b::Vector{T})::Vector{T} where T<:AbstractFloat
    n = size(L)[1]
    x = zeros(T, n)
    x[1] = b[1]/L[1, 1]
    for i in 2:1:n
        x[i] = b[i]
        for j in 1:1:(i-1)
            x[i] -= L[i, j]*x[j]
        end
        x[i] = x[i]/L[i, i]
    end
    return x
end


"""
    Us(A, b)

상삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Us(U::Matrix{T}, b::Vector{T})::Vector{T} where T<:AbstractFloat
    n = size(U)[1]
    x = zeros(T, n)
    x[n] = b[n]/U[n, n]
    for i in (n-1):-1:1
        x[i] = b[i]
        for j in (i+1):1:n
            x[i] -= U[i, j] * x[j]
        end
        x[i] = x[i]/U[i, i]
    end
    return x
end

"""
    L1(n, i, j)

n \tims m 행렬의 i번째 행과 j 번째 행을 교환하는 행 기본 연산
"""
function L1(n::Integer, i::Integer, j::Integer, T::Type = Float64)
    E = II(n, T)
    E[i, i] = zero(T)
    E[j, j] = zero(T)
    E[i, j] = one(T)
    E[j, i] = one(T)
    return E
end

"""
    iL1(n, i, j)

E1(n, i, j) 의 역행렬.. 이지만 E1(n, i, j) 와 같다.
"""
function iL1(n::Integer, i::Integer, j::Integer, T::Type = Float64)
    return L1(n, i, j, T)
end


"""
    L2(n, i, c)

n \times m 행렬의 i 번째 행에 상수 c 를 곱하는 행 기본 연산. 
"""
function L2(n::Integer, i::Integer, c::Number, T::Type = Float64)
    E = II(n, T)
    E[i, i] = convert(c, T)
    return E
end


"""
    iL2(n, i, c)

E2(n, i, c) 의 역행렬로 L2(n, i, 1.0/c) 와 같다.
"""
function iL2(n::Integer, i::Integer, c::Number, T::Type = Float64)
    return L2(n, i, 1.0/c, T)
end


"""
    L3(n, i, j, c)

n \times m 행렬의  i-번째 행을 i 번째행 + c \times j-번째 행으로 변환하는 행 기본 연산.
"""
function L3(n::Integer, i::Integer, j::Integer, c::Number, T::Type = Float64)
    E = II(n, T)
    E[i, j] = convert(T, c)
    return E
end


"""
    iL3(n, i, j, c)

L3(n, i, j, c) 의 역행렬로 L3(n, i, j, -c) 와 같다.
"""
function L3(n::Integer, i::Integer, j::Integer, c::Number, T::Type = Float64)
    return L3(n, i, j, -c)
end


function LU(M::Array{T, 2}) where T<:AbstractFloat
    U = copy(M)
    m, n = size(U)[1:2]
    @assert m == n
    L = II(n, T)
    @inbounds for i in 2:n
        L[i, 1] = M[i, 1]/M[1, 1]
        U[i, 1] = 0.0 
    end
    @inbounds for i = 2:n 
        for j = 2:n 
            fc = zero(T)
            for k = 1:i-1
                fc += L[i, k]*U[k, j]
            end
            if i>j
                L[i, j] = (M[i, j] - fc)/U[j, j]
                U[i, j] = zero(T)
            else 
                U[i, j] = (M[i, j] - fc) 
            end
        end
    end
    return L, U
end



function PLU(A::Matrix{T}, b::Vector{T}, pivoting=true) where T<:AbstractFloat
    U = copy(A)
    c = copy(b)
    m, n = size(U)[1:2]
    @assert m == n
    
    if pivoting
        P = II(n, T)
        for i in 1:n-1
            maxindex = argmax(U[i:n, i], dims=1)[1] +i -1
            if maxindex ≠ i
                for k in 1:1:n
                    U[i, k], U[maxindex, k] = U[maxindex, k], U[i, k]
                    P[i, k], P[maxindex, k] = P[maxindex, k], P[i, k]
                end
                c[i], c[maxindex] = c[maxindex], c[i]
            end
        end
    end
    L, U1 = LU(U)
    return P, L, U1, c
end

function solve_lu(A::Matrix{T}, b::Vector{T}, pivoting = true)::Vector{T} where T<:AbstractFloat
    P, L, U, c = PLU(A, b)
    y = Ls(L,  c)
    x = Us(U, y)
    return x
end


Ac = Float32.([3 4 2 5; 5 3 4 -2; -2 -3 1 0; -2 -3 1 -2])
b = Float32.([3, 2, 1, -2])
l1, u1 = LU(Ac)
p4, l4, u4, b4 = PLU(Ac, b)
x = solve_lu(Ac, b, true)
x2 = solve_lu(Ac, b, false)

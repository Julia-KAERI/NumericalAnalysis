using LinearAlgebra, BenchmarkTools 

"""
    II(n::Integer, T::Type = Float64)

n \times n 단위행렬을 리턴한다.
"""
function II(n::Integer, T::Type = Float64) 
    return Matrix{T}(I, n, n)
end

"""
    Ls(A, b)

하삼각행렬 L 에 대해 Lx=b 의 해 x 를 구한˜다.
"""
function Ls(L::Matrix{T}, b::Vector{T}) where T<:Number
    M, N = size(L)    
    @assert M == N
    x = zeros(T, N)
    x[1] = b[1]/L[1, 1]
    for i in 2:N
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

상삼각행렬 U 에 대해 Ux=b 의 해 x 를 구한다.
"""
function Us(U::Matrix{T}, b::Vector{T}) where T<:Number
    M, N = size(U)    
    @assert M == N
    x = zeros(T, N)
    x[N] = b[N]/U[N, N]
    for i in (N-1):-1:1
        x[i] = b[i]
        for j in (i+1):1:N
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
function L1(n, i, j)
    E = II(n)
    E[i, i] = 0.0
    E[j, j] = 0.0
    E[i, j] = 1.0
    E[j, i] = 1.0
    return E
end

"""
    iL1(n, i, j)

E1(n, i, j) 의 역행렬.. 이지만 E1(n, i, j) 와 같다.
"""
function iL1(n, i, j)
    return L1(n, i, j)
end


"""
    L2(n, i, c)

n \times m 행렬의 i 번째 행에 상수 c 를 곱하는 행 기본 연산. 
"""
function L2(n, i, c)
    E = II(n)
    E[i, i] = c
    return E    
end


"""
    iL2(n, i, c)

E2(n, i, c) 의 역행렬로 L2(n, i, 1.0/c) 와 같다.
"""
function iL2(n, i, c)
    return L2(n, i, 1.0/c)
end


"""
    L3(n, i, j, c)

n \times m 행렬의  i-번째 행을 i 번째행 + c \times j-번째 행으로 변환하는 행 기본 연산.
"""
function L3(n, i, j, c)
    E = II(n)
    E[i, j] =c
    return E
end


"""
    iL3(n, i, j, c)

L3(n, i, j, c) 의 역행렬로 L3(n, i, j, -c) 와 같다.
"""
function iL3(n, i, j, c)
    return L3(n, i, j, -c)
end

"""
    GE1(A)

Basic Gauss elimination without pivoting.
"""
function GE1(A::Matrix{T}) where T<:Number
    B = copy(M)
    m, n = size(A)[1:2]
    M = min(m, n)
    for i in 1:1:M, j in (i+1):1:m
        A[j, :] = A[j, :] .- (A[i, :] .* (A[j , i]/A[i, i]))
    end
    return A
end

A1 = Float64.([1 3 2; 3 -4 2; 3 3 1])
GE1(A1)

A2 = Float64.([1 2 3 ; 2 3 4 ; 3 4 5; 1 2 2])
GE1(A2)

A3 = Float64.([1.1 2.3 3 ; 2 3 4.0 ; 3 4 5.0; 1 2 2])
GE1(A3)

function find_nonzero_rows_after(M::Matrix{T}, i::Integer) where T<:Number
    h, w = size(M)
    k = 0
    for j in i+1:1:h
        if M[j, i] ≠ 0.0
            k = j
            break
        end
    end
    return k
end

function GE2(M::Array{T, 2})::Array{T, 2} where T<:AbstractFloat 
    A = copy(M)
    m, n = size(A)[1:2]
    M = min(m, n)
    for i in 1:1:M, j in (i+1):1:m
        if A[i, i] == 0.0
            k = find_nonzero_rows_after(A, i)
            if k ≠ nothing
                A[i, :], A[k, :] = A[k, :], A[i, :]
            else
                continue
            end
        end
        A[j, :] = A[j, :] .- (A[i, :] .* (A[j , i]/A[i, i]))
        A[j, i] = 0.0
    end
    return A
end


A1 = Float64.([1.1 3 2; 3 -4 2; 3 3 1])
GE1(A1)
A2 = Float64.([1 2 3 ; 2 3 4 ; 3 4 5; 1 2 2])
GE1(A2)
A3 = Float64.([1.1 2.3 3 ; 2 3 4.0 ; 3 4 5.0; 1 2 2])
GE1(A3)
GE2(A1)
A4 = Float64.([1.1 2.3 3 ; 2 3 4.0 ; 0 0 0 ;3 4 5.0; 1 2 2])
GE2(A4)

function GJE(M::Array{T, 2})::Array{T, 2} where T<:AbstractFloat
    A = GE2(M)
    m, n = size(A)[1:2]
    M = min(m, n)
    for i in 1:1:M
        if A[i, i] ≠ 0.0
            A[i, :] = A[i, :]./A[i, i]
        end
    end
    return A
end

GJE(A4)







A1 = Float64.([0 2 3 4;1 0 3 4; 1 3 0 1; 2 4 2 2]);
A2 = Float64.([1. 3. 1.; 0. 0. -0.; 3. 11. 5.]);
A3 = [1.357  0.124  0.628 ;  0.34  1.65  0.124 ; 0.148  0.822  1.112]
GJE(A3)





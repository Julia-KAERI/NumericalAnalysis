using LinearAlgebra, BenchmarkTools 

"""
    II(n)

n \times n 단위행렬을 리턴한다.
"""
# function II(n::T, mtype::N=Float64) where {T<:Integer, N<: Type}
#     if isabstracttype(mtype) || ~(mtype <: AbstractFloat) 
#         error("mtype should be concrete type of floating number")
#     end
    
#     return Matrix{mtype}(I, n, n)
# end

function II(n::Integer, T::Type = Float64) 
    return Matrix{T}(I, n, n)
end

"""
    Ls(A, b)

하삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Ls(L::AbstractMatrix{T}, b::Vector{T})::Vector{T} where T<:AbstractFloat
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
function Us(U::AbstractMatrix{T}, b::Vector{T})::Vector{T} where T<:AbstractFloat
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


function decomposition_lu(M::AbstractMatrix{T}) where T<:AbstractFloat
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



function decomposition_plu(A::AbstractMatrix{T}, b::AbstractVector{T}, pivoting=true) where T<:AbstractFloat
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
    L, U1 = decomposition_lu(U)
    return P, L, U1, c
end


"""
    solve_lu(A::AbstractMatrix{T}, b::Vector{T}, pivoting = true)::Vector{T} where T<:AbstractFloat

solve Ax = b then get x by LU decomposition of A. If pivoting == false, pivoting 
in LU decomposition is avoided
"""
function solve_lu(A::AbstractMatrix{T}, b::AbstractVector{T}, pivoting = true)::Vector{T} where T<:AbstractFloat
    P, L, U, c = decomposition_plu(A, b)
    y = Ls(L,  c)
    x = Us(U, y)
    return x
end

function gram_schmidt(A::T) where T<:AbstractMatrix{N} where N<:AbstractFloat
    m, n = size(A)
    U = zeros(N, (m, n))
    U[:,1] = A[:,1]/norm(A[:,1])
    for i = 2:1:n
        u = A[:, i]
        for k in 1:1:i-1
            uk = U[:, k]
            ck = (uk'*A[:,i])/(uk'*uk)
            u -= ck .* uk
        end
        U[:,i] = u/norm(u)
    end
    return U
end

function decomposition_qr(A::T) where T<:AbstractMatrix{N} where N<:AbstractFloat
    Q = gram_schmidt(A)
    R = zeros(N, size(A))
    m, n = size(A)
    for j in 1:1:n
        for i in 1:1:j
            for k = 1:1:m
                R[i, j] += (Q[k, i]*A[k, j])
            end
        end
    end
    return Q, R
end



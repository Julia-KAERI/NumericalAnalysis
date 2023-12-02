

"""
    Ls(A, b)

하삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Us(U::Matrix{T}, b::Vector{T}) where T<:Number
    M, N = size(U)    
    @assert M == N == length(b)
    
    if zero(T) in diag(U)
        error("Input matrix has zero diagonal element")
    end
    
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
    Us(A, b)

상삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Ls(L::Matrix{T}, b::Vector{T}) where T<:Number
    M, N = size(L)    
    @assert M == N == length(b)
       
    if zero(T) in diag(L)
        error("Input matrix has zero diagonal element")
    end
    
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
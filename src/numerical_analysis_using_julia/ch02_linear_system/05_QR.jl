using LinearAlgebra, BenchmarkTools

function Gram_Schmidt(A::Matrix{T})::Matrix{T} where T<:AbstractFloat
    m, n = size(A)
    @assert m == n
    U = zeros(T, (m, n))
    U[:,1] = A[:,1]/norm(A[:,1])
    for i = 2:n
        @views vi = A[:, i]
        vv = zeros(T, n)
        for k in 1:i-1
            @views uk = U[:, k]
            ck = dot(uk, vi)#/dot(uk, uk)
            vv .+= ck .* uk
        end
        un = vi .- vv
        U[:,i] = un/norm(un)
    end
    return U
end

function QR(A::Matrix{T}) where T<:AbstractFloat
    Q = Gram_Schmidt(A)
    R = zeros(T, size(A))
    m, n = size(A)
    for j in 1:1:n
        for i in 1:1:j
            @views R[i, j] = dot(Q[:, i], A[:, j])
        end
    end
    return Q, R
end

U = Float64.([1 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 -1 1])
V = Gram_Schmidt(U)
v1, v2, v3, v4 = V[:,1], V[:,2], V[:,3], V[:,4]
Q, R = QR(U)
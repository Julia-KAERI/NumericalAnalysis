using LinearAlgebra, BenchmarkTools 

"""
    II(n)

n \times n 단위행렬을 리턴한다.
"""
function II(n)
    return Matrix{Float64}(I, n, n)
end

function Gram_Schmidt(A::AbstractMatrix)
    m, n = size(A)
    U = zeros(Float64, (m, n))
    U[:,1]=A[:,1]/norm(A[:,1])
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

function QR(A::AbstractMatrix)
    Q = Gram_Schmidt(A)
    R = zeros(Float64, size(A))
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


Ac = Float64.([3 4 2 5; 5 3 4 -2; -2 -3 1 0; -2 -3 1 -2])
G = Gram_Schmidt(Ac)
q, r = QR(Ac)
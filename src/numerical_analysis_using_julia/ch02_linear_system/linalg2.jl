using LinearAlgebra, BenchmarkTools 

"""
    II(n)

n \times n 단위행렬을 리턴한다.
"""
function II(n)
    return Matrix{Float64}(I, n, n)
end

"""
    Ls(A, b)

하삼각행렬 A 에 대해 Ax=b 의 해 x 를 구한다.
"""
function Ls(L::AbstractMatrix, b::AbstractVector)
    n = size(L)[1]
    x = zeros(Float64, n)
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

function Ls2(L::Array{T, 2}, b::Vector{T}) where T<:Number
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
function Us(U::AbstractMatrix, b::AbstractVector)
    n = size(U)[1]
    x = zeros(Float64, n)
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


function LU(A::AbstractMatrix)
    U = copy(A)
    m, n = size(U)[1:2]
    @assert m == n
    L = II(n)
    for i = 2:n 
        L[i, 1] = A[i, 1]/A[1, 1]
        U[i, 1] = 0.0
        for j = 2:n
            fc = 0.0
            for k = 1:min(i, j)-1
                fc += L[i, k]*U[k, j]
            end
            if i>j
                L[i, j] = (A[i, j] - fc)/U[j, j]
                U[i, j] = 0.0
            else
                U[i, j] = A[i, j] - fc
            end
        end
    end
    return L, U
end


function LU1(A::AbstractMatrix)
    U = copy(A)
    m, n = size(U)[1:2]
    @assert m == n
    L = II(n)
    for k = 1:n
        L[k+1:n,k] = U[k+1:n,k] / U[k,k]
        U[k+1:n,:] = U[k+1:n,:] - L[k+1:n,k]*U[k:k,:]  
    end
    return L, U

end

function LU2(A::AbstractMatrix)
    U = copy(A)
    m, n = size(U)[1:2]
    @assert m == n
    L = II(n)
    L[2:n, 1] = A[2:n, 1]/A[1, 1]
    for j = 2:n, i = 2:n
        if i>j
            L[i, j] = (A[i, j] - L[i, 1:i-1]'*U[1:i-1, j])/U[j, j]
        else 
            U[i, j] = (A[i, j] - L[i, 1:i-1]'*U[1:i-1, j])
            U[i, 1:i-1] .= 0.0
        end
    end
    return L, U
end

function LU3(A::AbstractMatrix)
    U = copy(A)
    m, n = size(U)[1:2]
    @assert m == n
    L = II(n)
    for i in 2:n
        L[i, 1] = A[i, 1]/A[1, 1]
        U[i, 1] = 0.0 
    end
    for i = 2:n 
        for j = 2:n 
            fc = 0.0
            for k = 1:i-1
                fc += L[i, k]*U[k, j]
            end
            
            if i>j
                L[i, j] = (A[i, j] - fc)/U[j, j]
                U[i, j] = 0.0
            else 
                U[i, j] = (A[i, j] - fc) 
            end
        end
    end
    return L, U
end




function PLU(A::AbstractMatrix, b::AbstractVector, pivoting=true)
    U = copy(A)
    c = copy(b)
    m, n = size(U)[1:2]
    @assert m == n
    
    if pivoting
        P = II(n)
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


function solve_lu(A::AbstractMatrix, b::AbstractVector, pivoting = true)
    P, L, U, c = PLU(A, b)
    y = Ls(L,  c)
    x = Us(U, y)
    return x
end


Ac = Float64.([3 4 2 5; 5 3 4 -2; -2 -3 1 0; -2 -3 1 -2])
b = [3, 2, 1, -2]
l1, u1 = LU1(Ac)
l2, u2 =  LU2(Ac)
l3, u3 = LU3(Ac)
p4, l4, u4, b4 = PLU(Ac, b)
x = solve_lu(Ac, b, true)
x2 = solve_lu(Ac, b, false)


function det_PLU(A::AbstractMatrix)
    U = copy(A)
    m, n = size(U)[1:2]
    @assert m == n
    p = 1    
    for i in 1:n-1
        maxindex = argmax(U[i:n, i], dims=1)[1] +i -1
        if maxindex ≠ i
            for j = 1:1:n
                U[i, j], U[maxindex, j] = U[maxindex, j], U[i, j]
            end
            p *= -1
        end
    end
    L, U1 = LU(U)
    return p * prod(diag(L))*prod(diag(U1))
end
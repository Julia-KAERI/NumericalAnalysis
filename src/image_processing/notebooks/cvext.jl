using Distributions
using OpenCV

const cv = OpenCV

function arr2mat(arr::Matrix{T}) where T
    cv.Mat(permutedims(stack([arr, ]), [3,2,1]))
end

function arr2mat(arr::Array{T, 3}) where T
    cv.Mat(arr)
end

function img2arr(img)
    T = typeof(img[1, 1].val.i)
    broadcast(q->T(q.val.i), img)
end

function img2mat(img) 
    T = typeof(img[1, 1].val.i)
    tm = broadcast(q->T(q.val.i),img)
    cv.Mat(permutedims(stack([tm, ]), [3,2,1]))
end

function mat2arr(mat::OpenCV.Mat)
    return permutedims(mat.data, [3,2,1])
end

function histogram1d(mat::OpenCV.Mat{T}) where T<:Union{UInt8, UInt16}
    tm = Int32(typemax(T))

    v = cv.calcHist(cv.InputArray[mat,], Int32[0], fill(UInt8(1), size(mat)), Int32[tm+1], Float32[0, tm+1])
    return (0:1:tm, Int64.(v[1,1,:]))    
end



function gaussian_noise(img::OpenCV.Mat{T}, μ, σ, N) where T<:Union{UInt8, UInt16}
    m, n = size(img)[2:3]
    MM = typemax(T)
    Y, X, V = rand(1:m, N), rand(1:n, N), round.(T, rand(truncated(Normal(μ, σ), 0, MM), N))
    img2 = copy(img)
    for (y, x, v) ∈ zip(Y, X, V)
        img2[1, y, x] = v
    end
    return img2
end

function salt_pepper_noise(img::OpenCV.Mat{T}, N) where T<:Union{UInt8, UInt16}
    m, n = size(img)[2:3]
    MM = typemax(T)
    Y, X, V = rand(1:m, N), rand(1:n, N), sample([0, MM], N)
    img2 = copy(img)
    for (y, x, v) ∈ zip(Y, X, V)
        img2[1, y, x] = v
    end
    return img2
end


function cvPoint(x::T1, y::T2) where {T1<:Real, T2<:Real}
    T = promote_type(T1, T2)
    return cv.Point{T}(T(x), T(y))
end

function cvPoint(x::T1, y::T2, z::T3) where {T1<:Real, T2<:Real, T3<:Real}
    T = promote_type(T1, T2, T3)
    return cv.Point3{T}(T(x), T(y), T(z))
end

function cvSize(w::T1, h::T2) where {T1<:Integer, T2<:Integer}
    return cv.Size(Int32(w), Int32(h))
end;


function cvConvertTo(mat::OpenCV.Mat, t::T) where T<:Type
    return cv.Mat(convert.(t, mat.data))
end
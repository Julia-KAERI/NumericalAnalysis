using Plots, BenchmarkTools, LinearAlgebra



function bisection(f::Function, a::T, b::T, xtol::T=1e-10, maxiter::Integer=10_000) where T<:Real
    Nitter = 0
    a, b = minmax(a, b)
    c = zero(T)
    hf = one(T)/2
    f(a)*f(b) <= 0 || error("f(a)*f(b) should be negative") 
    
    while ((b-a) > xtol) 
        Nitter +=1
        Nitter < maxiter || error("Interation 횟수가 정해진 최댓값에 도달했습니다")
        c = (a+b)*hf
        if f(c) == 0.0
            break
        elseif f(a)*f(c) < 0 
            a, b = a, c
        else 
            a, b = c, b 
        end
    end
    println("Number of Iterlation = $Nitter")
    return c    
end

function newton_method(f::Function, df::Function, xi::T, MaxIter::Int64=100_000, etol::T = 1.0e-10, dfmin::T = 1.0e-6)::T where T <:AbstractFloat
    Niter = 0
    for i in 1:MaxIter
        if abs(f(xi)) < etol
            break
        elseif abs(df(xi)) < dfmin 
            error("df ≈ 0.0")
        end
        xi = xi - f(xi)/df(xi)
        Niter += 1       
        
        if abs(f(xi)) < etol 
            println("$Niter 회 반복 후 답 : $xi ")
            return xi
        end
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함")
end

function newton_secant(f::Function, x0::T, x1::T, MaxIter::Int64 = 100_000, etol::T = 1.0e-10, dfmin::T = 1.0e-10) where T <: AbstractFloat
    Niter = 0
    xp, xc = x0, x1

    for i in 1:MaxIter
        gx = (f(xc)-f(xp))/(xc-xp)
        if abs(f(xc)) < etol 
            return xc
        elseif abs(gx)<dfmin
            error("df ≈ 0.0")
        end
        xp, xc = xc, xc - f(xc)/gx
        Niter += 1
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함.")
end

function interp1d(xi::Vector{T}, yi::Vector{T}, method = "nearest", is_regular_grid = false) where T<:AbstractFloat
    @assert sizeof(xi) == sizeof(yi)
    if method == "nearest"
        return x->interp1d_nn(x, xi, yi, is_regular_grid)
    elseif method == "linear"
        return x->interp1d_linear(x, xi, yi, is_regular_grid)
    elseif method == "cubic"
        return x->interp1d_cubic(x, xi, yi, is_regular_grid)
    else 
        error("$method is not supported interpolation method")
    end
end

@doc raw"""
    interp1d_nn(x, xarr, yarr, is_regular_grid)

1D interpolation using nearest neighbor approximation    
"""
function interp1d_nn(x::T, xarr::Vector{T}, yarr::Vector{T}, is_regular_grid=false)::T where T<:AbstractFloat
    y = zero(T) 
    
    # Regular grid 의 경우
    if is_regular_grid 
        dx = (xarr[end] - xarr[1])/(length(xarr)-1)
        ind = round(Int64, (x-xarr[1])/dx +1 )
        if 1 <= ind <= length(xarr)
            y = yarr[ind]
        end
    
    # Irregular grid 의 경우
    else
        @inbounds for i in 1:1:(length(xarr)-1)
            if (x >= xarr[i]) && (x < xarr[i+1])
                if (x-xarr[i]) <= (xarr[i+1]-x)
                    y = yarr[i]
                else 
                    y = yarr[i+1]
                end
                break
            end
        end

        if x == last(xarr)
            y = last(yarr)
        end
    end
    return y
end

@doc raw"""
    interp1d_li(x, xarr, yarr, is_regular_grid)

1D interpolation using linear approximation    
"""
function interp1d_linear(x::T, xarr::Vector{T}, yarr::Vector{T}, is_regular_grid = false) where T<:AbstractFloat
    y = zero(T)

    # Regular grid 의 경우
    if is_regular_grid
        dx = (xarr[end] - xarr[1])/(length(xarr)-1)
        ind::Int64 = fld( x-xarr[1], dx)+1
        if 1 <= ind < length(xarr)
            d = x-xarr[ind]
            y = (yarr[ind+1]-yarr[ind])/dx * d + yarr[ind]
        end
    
    # Irregular grid의 경우
    else   
        @inbounds for i in 1:1:(length(xarr)-1)
            if x == xarr[i]
                y = yarr[i]
            elseif (x-xarr[i])*(x  - xarr[i+1]) < 0.0
                y = (yarr[i+1]-yarr[i])/(xarr[i+1]-xarr[i]) * (x-xarr[i]) + yarr[i]
                break
            end
        end
    end

    if x == last(xarr)
        y = last(yarr)
    end
    
    return y
end

@doc raw"""
    interp1d_cubic(x, xarr, yarr, is_regular_grid)

1D qubic interpolation 
"""
function interp1d_cubic(x::T, xarr::Vector{T}, yarr::Vector{T}, is_regular_grid)::T where T<:AbstractFloat
    y = zero(T)
    
    ll = length(xarr)
    xl = Array{T, 1}(undef, 4)
    yl = Array{T, 1}(undef, 4)

    if (x-xarr[1]) * (x - xarr[end]) > 0.0
        return y
    end

    if is_regular_grid
        dx = (xarr[end] - xarr[1])/(length(xarr)-1)
        ind::Int64 = fld( x-xarr[1], dx)+1
        if 1 <= ind < 3
            xl = first(xarr, 4)
            yl = first(yarr, 4)
        elseif ll-2 <= ind <= ll
            xl = last(xarr, 4)
            yl = last(yarr, 4)
        elseif 2 < ind < ll-2
            xl = xarr[ind-1:ind+2]
            yl = yarr[ind-1:ind+2]
        else 
            return y
        end
    else
        
        if (x>=xarr[1]) && (x<xarr[3])
            xl = first(xarr, 4)
            yl = first(yarr, 4)
        elseif (x >= xarr[ll-2]) && (x <= xarr[ll])
            xl = last(xarr, 4)
            yl = last(yarr, 4)
        else
            @inbounds for i in 3:1:(ll-3)
                if (x>= xarr[i]) && (x<xarr[i+1])
                    xl = xarr[i-1:i+2]
                    yl = yarr[i-1:i+2]
                end
            end
        end
    end

    if length(xl) == 0
        return 0.0
    else 
        y = yl[1]*(x-xl[2])*(x-xl[3])*(x-xl[4])/(xl[1]-xl[2])/(xl[1]-xl[3])/(xl[1]-xl[4])
        y += yl[2]*(x-xl[1])*(x-xl[3])*(x-xl[4])/(xl[2]-xl[1])/(xl[2]-xl[3])/(xl[2]-xl[4])
        y += yl[3]*(x-xl[1])*(x-xl[2])*(x-xl[4])/(xl[3]-xl[1])/(xl[3]-xl[2])/(xl[3]-xl[4])
        y += yl[4]*(x-xl[1])*(x-xl[2])*(x-xl[3])/(xl[4]-xl[1])/(xl[4]-xl[2])/(xl[4]-xl[3])
        return y
    end
end


function integrate_simple(f::Function, a::T, b::T, N::Integer) where T<:AbstractFloat
    a, b = minmax(a, b)
    dx = (b-a)/N
    result = f(a+(b-a)/N/2)*dx
    for k in 1:N
        result += f(a+(b-a)/N/2 + k*(b-a)/N)*dx
    end
    return sign(b-a)*result 
end


function ode_euler(f::Function, x0::Real, y0::Real, Npoints::Integer, h = 1.0e-6)
    x0, y0, h = promote(x0, y0, h)
    xp = x0 .+ collect(0:1:(Npoints-1)) * h
    yp = zero(xp)
    yp[1] = y0
    for i in 2:Npoints
        yp[i] = yp[i-1]+ f(xp[i-1]) * h
    end
    return xp, yp
end


function ode_runge_kutta2(f::Function, x0::Real, y0::Real, Npoints::Integer, h = 1.0e-6) 
    x0, y0, h = promote(x0, y0, h)
    xp = x0 .+ collect(0:1:(Npoints-1)) * h
    yp = zero(xp)
    yp[1] = y0
    for i in 2:Npoints
        xn, yn = xp[i-1], yp[i-1]
        k1 = f(xn, yn)
        k2 = f(xn + h/2, yn + k1*h/2)
        yp[i] = yn + (k1+k2)*h/2
    end
    return xp, yp
end


function ode_runge_kutta4(f::Function, x0::Real, y0::Real, Npoints::Integer, h = 1.0e-6) where T<:AbstractFloat
    x0, y0, h = promote(x0, y0, h)
    xp = x0 .+ collect(0:1:(Npoints-1)) * h
    yp = zero(xp)
    yp[1] = y0
    for i in 2:Npoints
        xn, yn = xp[i-1], yp[i-1]
        k1 = f(xn, yn)
        k2 = f(xn + h/2, yn + h*k1/2)
        k3 = f(xn + h/2, yn + h* k2/2)
        k4 = f(xn + h, yn + h*k3)
        yp[i] = yn + (k1 + 2*k2 + 2*k3 + k4)*h/6
    end
    return xp, yp
end

function ode_runge_kutta4_2(f::Function, t0::T, x0::Vector{T}, Npoints::Integer, epsilon = 1.0e-6) where T<:AbstractFloat
    N = length(x0)
    eps = convert(T, epsilon)
    tp = t0 .+ collect(0:1:(Npoints-1)) * eps
    xp = zeros(T,  N, Npoints)
    xp[:,1] = x0[:]
    for i in 2:Npoints
        tn, xn = t[i-1], xp[:, i-1]
        k1 = f(tn, xn)
        k2 = f(tn + eps/2, xn .+ (eps/2) .* k1)
        k3 = f(tn + eps/2, xn .+ (eps/2) .* k2)
        k4 = f(tn + eps, xn .+ eps.*k3)
        xp[:, i] = xn .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4).*(eps/6)
    end
    return tp, xp
end
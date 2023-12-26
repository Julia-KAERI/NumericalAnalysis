using CairoMakie, LaTeXStrings

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


f2(x) = x*(x-1.)*(x-2.)+5
df2(x) = 3*x^2-6*x+2
newton_method(f2, df2, 0.0)
newton_secant(f2, 0.0, 0.1)

function fixedpoint(f::Function, xi::T, xtol::T=1.0e-10, MaxIter::Int64 = 100_000)::T where T<:AbstractFloat
    Niter = 0
    for i in 1:MaxIter
        c = f(xi)
        if abs(c-xi) < xtol
            return c
        else
            xi = c
            Niter += 1
        end
    end
    error("최대 반복 횟수 $MaxIter 에 도달하였으나 답을 찾지 못함.")
end

h(x) = (2x^2+1)^(1/3)
fixedpoint(h, 2.0)
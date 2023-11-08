# julia-dec
using Plots
using IterativeSolvers
using LinearAlgebra
function newton(x0, f, Df, tol, itmax)
    x = x0
    if typeof(x0)==float
        return newton1d(x0,f,tol,itmax) 
    else
        n = length(x0) end
    delta_x = zeros(n)
    for i in range(itmax)
        J = Df(x,f)
        if det(J) == 0
            print("Singularity encountered")
            return "Jacobian is singular. Choose different x0" end
        cg!(delta_x,J, -f(x)) # Solving LinSys

        x += delta_x

        if norm(delta_x) < tol
            return x end
    end  
end
function newton(f,x0,tol,df,its)
    for _ in 0:its
        fx_0 = f(x0)
        dir = df(x0,f)
        if dir == 0
            return "Newton's method failed: derivative close to zero" end
        x0 = x0 - fx_0 / dir
        fx_1 = f(x0)
        rel_error = abs(fx_1 / fx_0)
        if rel_error < tol
            return x0 end
    end
    return -1
end

function df(x,f)
    h = 0.0000001
    return (f(x+h)-f(x))/h
end

function f1(x)
    return sin(x)
end

print(round(newton(f1,3,0.00001,df,1000),digits=15))

    
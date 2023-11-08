using Plots
using LinearAlgebra

function explicit_euler(func::Function,t0::Int64,y0::Float64,N::Int64,tn::Int64)
	h::Float64 = (tn-t0)/N
	times::Array = range(t0,tn,N)
	res::Array = zeros(N)
	res[1] = y0
	for i in range(1,N-1)
		res[i+1] = res[i] + h*func(times[i],res[i])
	end
	return times,res
end

function heun(func::Function,t0::Int64,y0::Float64,N::Int64,tn::Int64)
	h::Float64 = (tn-t0)/N
	times::Array = range(t0,tn,N)
	res::Array = zeros(N)
	res[1] = y0
	for i in range(1,N-1)
		xp = res[i] + h*func(times[i],res[i])
        res[i+1] = res[i]+h*(func(times[i],res[i])+func(times[i+1],xp))/2
	end
	return times,res,error(times,res)
end 

function error(times,res)
	an = f_a.(times)
	return abs.(an-res)
end

function f(t,y)
	return 3*exp.(-4*t)-2*y
end

function f_a(t)
	return  0.5*(5*exp(-2*t)-3*exp(-4*t))
end
n = range(1,3000)
x,ys = explicit_euler(f,0,1.,3000,3)
err = error(x,ys)
p = plot(n,err)

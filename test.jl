using Plots
theme(:default)

function explicit_euler(func, t0, y0, N, tn)
    h = (tn - t0) / N
    times = LinRange(t0, tn, N)  # Use LinRange instead of range for better performance
    res = similar(times)        # Preallocate res array
    res[1] = y0
    for i in 1:N-1
        res[i+1] = res[i] + h * func(times[i], res[i])
    end
    return times, res
end

function heun(func, t0, y0, N, tn)
    h = (tn - t0) / N
    times = LinRange(t0, tn, N)  # Use LinRange instead of range for better performance
    res = similar(times)        # Preallocate res array
    res[1] = y0
    for i in 1:N-1
        xp = res[i] + h * func(times[i], res[i])
        res[i+1] = res[i] + h * (func(times[i], res[i]) + func(times[i+1], xp)) / 2
    end
    return times, res
end 

function f(t, y)
    return 3 * exp(-4 * t) - 2 * y
end

function f_a(t)
    return 0.5 * (5 * exp(-2 * t) - 3 * exp(-4 * t))
end

function erra(func, times, res)
    an = func.(times)
    return abs.(an - res)
end

function plot_conv()
	t0 = 0.0
	tn = 1.0
	n0 = 1
	n1 = 3
    N_vals = [10^i for i in n0:n1]
    mini_euler = [maximum(erra(f_a, explicit_euler(f, t0, tn, N, 1.0)...)) for N in N_vals]
    mini_heun = [maximum(erra(f_a, heun(f, t0, tn, N, 1.0)...)) for N in N_vals]
    plot(N_vals, mini_euler, color="blue", label="Euler", marker=:circle)
    plot!(N_vals, mini_heun, color="red", label="Heun", marker=:circle, xlabel="N", ylabel="Max Error")
    plot!(N_vals, [(((tn-t0)/N_vals[i]))^2 for i in n0:n1], color="pink", label="Order 2", marker=:circle)
	plot!(xaxis=:log, yaxis=:log)
	
end

function ploterr()
    N = 40000
    n = range(1, 1000)
    x, ys = explicit_euler(f, 0.03, 1.0, 1000, 1.0)
    err = erra(f_a, x, ys)

    p = plot(n, err, ylim=[1e-10, 1.0], xlabel="n", ylabel="erra")
    plot!(xaxis=:log, yaxis=:log)
end

plot_conv()
ploterr()

using Plots

h = 1e-3

function reduce_with_prev(f, a, init)
    s = collect(similar(a, Tuple{Any, Any}))
    isempty(s) && return(s)
    s[1] = init
    for i in firstindex(s)+1:lastindex(s)
        s[i] = f(s[i-1], a[i])
    end
    return s
end

function explicit_euler(a, b, n, f, init)
    h = (b - a)/n
    xs = range(a, b, n)
    reduce_with_prev(((y_, x) -> y_ .+ f(x, y_) .* h), xs, init)
end

function runge_kutta_3(a, b, n, f, init)
    h = (b - a)/n
    xs = range(a, b, n)
    reduce_with_prev((y_, x) -> begin
        ϕ0 = f(x, y_)
        ϕ1 = f(x + h/3, y_ .+ ϕ0.*(h/3))
        ϕ2 = f(x + 2h/3, y_ .+ ϕ1.*(2h/3))
        y_ .+ ((ϕ0 .+ ϕ2 .* 3) .* h)./4
    end, xs, init)
end

a = 1
b = 1.5
f = (x, y) -> begin
    u, v = y
    u_ = u^3 - v^2 - x
    v_ = v^3 + u^2 + x
    (u_, v_)
end

n = 1000
fee = explicit_euler(a, b, n, f, (0.0, 0.0))
frk = runge_kutta_3(a, b, n, f, (0.0, 0.0))

println(fee[end] .- frk[end])

xs = range(a, b, n)
plot(xs, fee .|> (x -> x[1]), label="u, Explicit Euler")
plot!(xs, fee .|> (x -> x[2]), label="v, Explicit Euler")
plot!(xs, frk .|> (x -> x[1]), label="u, Runge-Kutta 3")
plot!(xs, frk .|> (x -> x[2]), label="v, Runge-Kutta 3")
using Plots

function quad_euler(a, b, fs, n)
    h = (b - a) / n
    f_0 = (fs[2] - fs[1]) / h
    f_n = (fs[n] - fs[n-1]) / h
    k = (1 / 12) * h^2 * (f_0 - f_n)
    h * sum(2:n .|> i -> (fs[i] + fs[i-1]) / 2 + k)
end

function runge_rule(a, b, m4, ε)
    Int(round(((b - a)^2 * m4) / (720 * ε)))
end

f1(x) = cos(x)^2
f1_(x) = x / 2 + sin(2x) / 4

f1_(π) - f1_(0)

a1 = 0
b1 = π
ε = 1e-6
f1_4(x) = 8(cos(x)^2 - sin(x)^2)
n = get_n(a1, b1, maximum(f1_4.(range(a1, b1, 1000))), ε)
fe = quad_euler(a1, b1, f1.(range(a1, b1, n)), n)

function chebyshev_range(a, b, n)
    0:(n-1) .|> (k -> begin
        (a + b) / 2 + ((b - a) / 2 * cos((π * (2k + 1)) / 2n))
    end)
end

f2(x) = exp(x) / sqrt(1 - x^2)
f2hermit(n) = π / (n + 1) * sum(chebyshev_range(-1, 1, n) .|> exp)
f2avgrect(n) = 2 / n * sum(range(-1, 1, n)[1:end-1] .|> x -> f2(x + 1 / n))
f2avgrect(10000000)
f2hermit(10000000)

ns = 10:100:100000
f2hermit(10000) - f2avgrect(10000)
plot(ns, [f2hermit.(ns), f2avgrect.(ns)])
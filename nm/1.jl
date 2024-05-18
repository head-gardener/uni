using Plots
using LinearAlgebra

# jacobi method
function jacobi(a, f)
  g, b, x = canonize(a, f)
  converge_(x, x -> b * x + g)
end

# canonize Ax=f to x=Bx+g using jacobi method,
# returning g, b and an initial approximation - x₀
function canonize(a, f)
  # restriction lifted due to the specifics of the task.
  # @assert(has_sdd(a), "sdd required!")
  @assert(size(a, 1) == size(f, 1), "sizes don't match")

  g = [f[i] / a[i, i] for i in axes(a, 1)]
  b = [i == j ? 0 : -a[i, j] / a[i, i]
       for i in axes(a, 1), j in axes(a, 2)]
  x₀ = ones(size(a, 1))

  (g, b, x₀)
end

# apply f to x₀ until difference falls below ε,
# or number of iteration exceeds k.
# returns resulting x and number of iterations passed.
function converge(x₀, f, ε, k)
  i = 1
  x = x₀

  while i <= k
    x₀ = x
    x = f(x)
    i += 1
    if abs(norm(x) - norm(x₀)) < ε
      break
    end
  end

  return (x, i)
end

# converge with ε and k set to some defaults
converge_(x₀, f) = converge(x₀, f, 1e-6, 1000)

function f1(x)
  x * sin(x)
end

function f2(x)
  sin(abs(x)) + 1
end

function div_diff(fs, xs, a, b)
  if a == b
    fs[a]
  else
    (div_diff(fs, xs, a + 1, b) - div_diff(fs, xs, a, b - 1)) /
    (xs[b] - xs[a])
  end
end

# nth_newton(n, fs, xs) = (x-x_0)...(x-x_n-1)f(x_0;...;x_n)
function nth_newton(n, fs, xs)
  if n == 0
    _ -> fs[1]
  else
    x -> prod(xs[1:min(n, end)] .|> (x_ -> x - x_)) *
         div_diff(fs, xs, 1, n + 1)
  end
end

function newton_poly(fs, xs, n_, n)
  cs = n_:n .|> (x -> nth_newton(x, fs, xs))
  (x -> cs .|> (f -> f(x)) |> sum)
end

function extend_newton_poly(p, fs, xs, n_, n)
  n = newton_poly(fs, xs, n_, n)
  (x -> p(x) + n(x))
end

function chebyshev_range(a, b, n)
  0:(n-1) .|> (k -> begin
    (a + b) / 2 + ((b - a) / 2 * cos((π * (2k + 1)) / 2n))
  end)
end

function cubic_spline(fs, xs, n)
  M_ = Matrix{Float64}(I, n, n)
  Mf = zeros(n)
  for i in 2:n-1
    h = xs[i] - xs[i-1]
    h_ = xs[i+1] - xs[i]
    M_[i, i-1] = h / 6
    M_[i, i] = (h + h_) / 3
    M_[i, i+1] = h_ / 6
    Mf[i] = (fs[i+1] / fs[i]) / h_ - (fs[i] / fs[i-1]) / h
  end
  # notice how jacobi method requires sdd,
  # which always holds for a uniform grid.
  (M, _) = jacobi(M_, Mf)
  (x -> begin
    i = max(2, findfirst(x_ -> x_ >= x, xs))
    h = xs[i] - xs[i-1]
    M[i-1] * ((xs[i] - x)^3 / 6h) +
    M[i] * ((x - xs[i-1])^3 / 6h) +
    (fs[i-1] - M[i-1] * (h^2 / 6)) * ((xs[i] - x) / h) +
    (fs[i] - M[i] * (h^2 / 6)) * ((x - xs[i-1]) / h)
  end)
end

x = range(-2, 2, length=100)
n = 5
# x_ = chebyshev_range(-2, 2, n)
y1 = f1.(x)
y2 = f2.(x)
# y2 = cubic_spline(y_, x_, n).(x)

function make_newton_poly(r, f, a, b, n)
  x_ = r(a, b, n+1)
  y_ = f.(x_)
  newton_poly(y_, x_, 0, n)
end

function make_cubic_spline(f, a, b, n)
  x_ = range(a, b, n)
  y_ = f.(x_)
  cubic_spline(y_, x_, n)
end

function to_analytical_rep(r, f, a, b)
  xs = r(a, b, 3)
  fs = f.(xs)
  "$(fs[1]) + (x - ($(xs[1])))*($(div_diff(fs, xs, 1, 2)))" *
  "(x - ($(xs[1])))(x - ($(xs[2])))*($(div_diff(fs, xs, 1, 3)))"
end

function plot_(x, f1, f2, label1, label2)
  p = plot(x, f1.(x), label=label1)
  plot!(p, x, f2.(x), label=label2)
end

to_analytical_rep(range, f1, -2, 2)

plot_(x, f1, make_newton_poly(range, f1, -2, 2, 2), "f1(x)", "P2(x)")
plot_(x, f1, make_newton_poly(range, f1, -2, 2, 4), "f1(x)", "P4(x)")
plot_(x, f1, make_newton_poly(range, f1, -2, 2, 8), "f1(x)", "P8(x)")
plot_(x, f1, make_newton_poly(range, f1, -2, 2, 16), "f1(x)", "P16(x)")

to_analytical_rep(range, f2, -2, 2)

plot_(x, f2, make_newton_poly(range, f2, -2, 2, 2), "f2(x)", "P2(x)")
plot_(x, f2, make_newton_poly(range, f2, -2, 2, 4), "f2(x)", "P4(x)")
plot_(x, f2, make_newton_poly(range, f2, -2, 2, 8), "f2(x)", "P8(x)")
plot_(x, f2, make_newton_poly(range, f2, -2, 2, 16), "f2(x)", "P16(x)")

to_analytical_rep(chebyshev_range, f1, -2, 2)

plot_(x, f1, make_newton_poly(chebyshev_range, f1, -2, 2, 2), "f1(x)", "P2(x)")
plot_(x, f1, make_newton_poly(chebyshev_range, f1, -2, 2, 4), "f1(x)", "P4(x)")
plot_(x, f1, make_newton_poly(chebyshev_range, f1, -2, 2, 8), "f1(x)", "P8(x)")
plot_(x, f1, make_newton_poly(chebyshev_range, f1, -2, 2, 16), "f1(x)", "P16(x)")

to_analytical_rep(chebyshev_range, f2, -2, 2)

plot_(x, f2, make_newton_poly(chebyshev_range, f2, -2, 2, 2), "f2(x)", "P2(x)")
plot_(x, f2, make_newton_poly(chebyshev_range, f2, -2, 2, 4), "f2(x)", "P4(x)")
plot_(x, f2, make_newton_poly(chebyshev_range, f2, -2, 2, 8), "f2(x)", "P8(x)")
plot_(x, f2, make_newton_poly(chebyshev_range, f2, -2, 2, 16), "f2(x)", "P16(x)")

plot_(x, f1, make_cubic_spline(f1, -2, 2, 3), "f1(x)", "S3(x), n = 3")
plot_(x, f1, make_cubic_spline(f1, -2, 2, 5), "f1(x)", "S3(x), n = 5")
plot_(x, f1, make_cubic_spline(f1, -2, 2, 9), "f1(x)", "S3(x), n = 9")
plot_(x, f1, make_cubic_spline(f1, -2, 2, 17), "f1(x)", "S3(x), n = 17")

plot_(x, f2, make_cubic_spline(f2, -2, 2, 3), "f2(x)", "S3(x), n = 3")
plot_(x, f2, make_cubic_spline(f2, -2, 2, 5), "f2(x)", "S3(x), n = 5")
plot_(x, f2, make_cubic_spline(f2, -2, 2, 9), "f2(x)", "S3(x), n = 9")
plot_(x, f2, make_cubic_spline(f2, -2, 2, 17), "f2(x)", "S3(x), n = 17")
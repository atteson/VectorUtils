using VectorUtils
using Random
using StatsBase

Random.seed!(1)

n = 100
xstar = rand(n);
x0 = zeros(n);
x1 = ones(n);

g(x) = x - xstar

(x,y) = brent( g, x0, x1 )
@assert( maximum(abs.(x - xstar)) .< 1e-6 )

g(x) = (x - xstar).^3
(x,y) = brent( g, x0, x1 )
@assert( maximum(abs.(x - xstar)) .< 1e-5 )
@assert( maximum(abs.(y - g(xstar))) .< 1e-10 )

m = 1_000
n = 10_000
mu = 0.0002 .+ 0.0001.*randn(m);
sigma = exp.(log(0.01) .+ 0.5.*randn(m));
r = mu .+ sigma .* randn(m, n);

g(r) = w -> mean(log.(1 .+ w .* r), dims=2)
x0 = vec(-1 ./ maximum(r, dims=2));
x0 += eps.(x0);
x1 = vec(-1 ./ minimum(r, dims=2));
x1 -= eps.(x1);
@assert( all(x0 .* x1 .< 0) )
@assert(sum(isinf.(g(r)(x0)))==0)
@assert(sum(isinf.(g(r)(x1)))==0)

(x, y) = brent( g(r), x0, x1 );
@time (x, y) = brent( g(r), x0, x1 );

using CUDA
rc = CuArray(r);
x0c = CuArray(x0);
x1c = CuArray(x1);
(xc, yc) = brent( g(rc), x0c, x1c );
@time (xc, yc) = brent( g(rc), x0c, x1c );


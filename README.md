# Normal Hermite Splines

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IgorKohan.github.io/NormalHermiteSplines.jl/stable)
[![Build Status](https://travis-ci.com/IgorKohan/NormalHermiteSplines.jl.svg?branch=master)](https://travis-ci.com/github/IgorKohan/NormalHermiteSplines.jl)
[![codecov.io](https://codecov.io/github/IgorKohan/NormalHermiteSplines.jl/coverage.svg?branch=master)](https://codecov.io/github/IgorKohan/NormalHermiteSplines.jl?branch=master)


![Problem definition](/images/problem-definition.png)

## Instalation
Start Julia and run the following commands:
```
julia> using Pkg
julia> Pkg.add("NormalHermiteSplines")
```
## Example usage
To use this package, begin your code with
```
using NormalHermiteSplines
```
The first example is the function ```φ(x,y) = sin(4.0 * sqrt(x^2 + y^2))```

 <img src="/images/m_t_6.png" width="256"/>  <img src="/images/m_cf_6.png" width="256"/>  <img src="/images/m_grid_6,2,3.png" width="256"/> 

We'll construct an interpolating normal spline using this function and its gradient values sampled on 1000 Halton nodes ([1]) distributed in the [-1,1]x[-1,1] square.
```
using NormalHermiteSplines

nodes = get_2D_halton_nodes(1000)             # generates Halton data set in [0,1]x[0,1] 
n_1 = size(nodes, 2)
u = Vector{Float64}(undef, n_1)               # function values
grid = get_2D_grid2(100)                      # uniform Cartesian grid of size 101x101 
d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)  # directional derivative nodes 
es = Matrix{Float64}(undef, 2, 2 * n_1)       # derivative direction 
du = Vector{Float64}(undef, 2 * n_1)          # directional derivative values 

k = 0
for i = 1:n_1
    nodes[1, i] = nodes[1, i] * 2.0 - 1.0     # transforming Halton nodes to [-1,1]x[-1,1] 
    nodes[2, i] = nodes[2, i] * 2.0 - 1.0
    d = sqrt(nodes[1, i]^2 + nodes[2, i]^2)
    u[i] = sin(4.0 * d)
    k += 1
    grad = [0.0; 0.0]
    if d > 0.0
        grad[1] = 4.0 * nodes[1, i] * cos(4.0 * d) / d
        grad[2] = 4.0 * nodes[2, i] * cos(4.0 * d) / d
    end
    d_nodes[1,k] = nodes[1,i]
    d_nodes[2,k] = nodes[2,i]
    du[k] = grad[1]
    es[1,k] = 1.0
    es[2,k] = 0.0
    k += 1
    d_nodes[1,k] = nodes[1,i]
    d_nodes[2,k] = nodes[2,i]
    du[k] = grad[2]
    es[1,k] = 0.0
    es[2,k] = 1.0
end

# Hermite spline must be constructed with RK_H1 or RK_H2 kernel.
# Here value of the 'scaling parameter' ε is estimated in the interpolate procedure.
rk = RK_H1()               
#
spline = interpolate(nodes, u, d_nodes, es, du, rk)
σ = evaluate(spline, grid)
```
The spline surface and filled 2-D contour plots:

 <img src="/images/s_t_6,2,3,1,_0.0,_.png" width="256"/>  <img src="/images/s_cf_6,2,3,1,_0.0,_.png" width="256"/> 

Approximation error plots:

 <img src="/images/delta_s_6,2,3,1,_0.0,_.png" width="256"/>  <img src="/images/delta_cf_6,2,3,1,_0.0,_.png" width="256"/> 

Spline was evaluated on a uniform Cartesian grid of size 101x101. Accuracy of the interpolation was measured by calculating the Root Mean Square Error (RMSE) and the Maximum Absolute Error (MAE). For this case
```RMSE```: 1.5E-03, ```MAE```: 1.1E-01, estimated value of the scaling parameter ```ε``` is 5.7E+00, estimation of the Gram matrix condition number is 1.0E+12.

The second example is the geometric composition ```Ψ``` of a circle and a square. The function ```ψ``` on ```Ω=[0,1]x[0,1]``` is given as ```Ψ = C + S```, where ```C``` and ```S``` denote the characteristic functions of these circle and square.

<img src="/images/m_t_1.png" width="256"/>  <img src="/images/m_cf_1.png" width="256"/>  <img src="/images/m_grid_1,2,5.png" width="256"/> 

We'll construct an interpolating normal spline using function ```Ψ``` values sampled on the 5000 Halton nodes distributed in the [0,1]x[0,1] square.
```
using NormalHermiteSplines

nodes = get_2D_halton_nodes(5000)             # generates Halton data set in [0,1]x[0,1] 
n_1 = size(nodes, 2)
u = Vector{Float64}(undef, n_1)               # function values
grid = get_2D_grid(100)                       # uniform Cartesian grid of size 101x101
for i = 1:n_1
    val = 0.0
    x = nodes[1,i]
    y = nodes[2,i]
    d = (x - 0.75)^2 + (y - 0.75)^2
    r2 = 0.04
    if d <= r2
        val = 1.0
    end
    if x <= 0.25 && x >= 0.05 && y <= 0.25 && y >= 0.05
        val = 1.0
    end
    u[i] = val
end

# Here spline is being constructed with ```RK_H0``` kernel,
# value of the 'scaling parameter' ```ε``` is defined explicitly.
rk = RK_H0(0.001)
#
spline = interpolate(nodes, u, rk)
σ = evaluate(spline, grid)
```
The spline surface and filled 2-D contour plots:

 <img src="/images/s_t_1,2,5,0,_0.001,_.png" width="256"/>  <img src="/images/s_cf_1,2,5,0,_0.001,_.png" width="256"/> 

Approximation error plots:

 <img src="/images/delta_s_1,2,5,0,_0.001,_.png" width="256"/>  <img src="/images/delta_cf_1,2,5,0,_0.001,_.png" width="256"/> 

Spline was evaluated on a uniform Cartesian grid of size 101x101. For this case ```RMSE```: 7.3E-02, ```MAE```: 0.95, value of the scaling parameter ```ε``` is 0.001, estimation of the Gram matrix condition number is 1.0e+13.

Further examples are given in documentation.

## Documentation

For more information and explanation see [Documentation](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Hermite-Splines/)

## References:
1. [Halton sequence](https://en.wikipedia.org/wiki/Halton_sequence)









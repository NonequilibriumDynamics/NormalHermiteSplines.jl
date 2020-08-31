<div align="center">
  <img src="/images/logo.png" width="400" alt="Normal Splines">
</div>

# Normal Hermite Splines

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://IgorKohan.github.io/NormalHermiteSplines.jl/v0.3.0)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IgorKohan.github.io/NormalHermiteSplines.jl/stable)
[![Build Status](https://travis-ci.com/IgorKohan/NormalHermiteSplines.jl.svg?branch=master)](https://travis-ci.com/github/IgorKohan/NormalHermiteSplines.jl)
[![codecov.io](https://codecov.io/github/IgorKohan/NormalHermiteSplines.jl/coverage.svg?branch=master)](https://codecov.io/github/IgorKohan/NormalHermiteSplines.jl?branch=master)


![Problem definition](/images/problem-definition.gif)

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
The first example is the function ```φ(x,y)=sin(4.0*sqrt(x^2+y^2))``` defined on ```Ω=[-1,1]x[-1,1]```.

 <img src="/images/m_t_6.png" width="256"/>  <img src="/images/m_cf_6.png" width="256"/>  <img src="/images/m_grid_6,2,3.png" width="256"/> 

We'll construct an interpolating normal spline using this function and its gradient values sampled on 1000 Halton nodes ([1]) distributed in the [-1,1]x[-1,1] square.
```
    using NormalHermiteSplines

    nodes = get_2D_halton_nodes(1000)           # generates 1000 Halton nodes in [0,1]x[0,1] 
    u = Vector{Float64}(undef, n_1)             # function values
    d_nodes = Matrix{Float64}(undef, 2, 2*n_1)  # directional derivative nodes 
    es = Matrix{Float64}(undef, 2, 2*n_1)       # derivative directions 
    du = Vector{Float64}(undef, 2*n_1)          # directional derivative values 

    grid = get_2D_grid2(100)                    # creates uniform Cartesian grid of size 101x101 

    n_1 = size(nodes, 2)
    k = 0
    for i = 1:n_1
        nodes[1, i] = nodes[1, i]*2.0-1.0       # transforming Halton nodes to [-1,1]x[-1,1] 
        nodes[2, i] = nodes[2, i]*2.0-1.0
        d = sqrt(nodes[1, i]^2 + nodes[2, i]^2)
        u[i] = sin(4.0*d)
        k += 1
        grad = [0.0; 0.0]
        if d > 0.0
            grad[1] = 4.0*nodes[1, i]*cos(4.0*d)/d
            grad[2] = 4.0*nodes[2, i]*cos(4.0*d)/d
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

    σ1 = evaluate_one(spline, [0.5; 0.5])
    #  ≈ 0.308071

    g1 = evaluate_gradient(spline, [0.5; 0.5])
    #  ≈ -2.690
    #  ≈ -2.690
```
The spline surface and filled 2-D contour plots:

 <img src="/images/s_t_6,2,3,1,_0.0,_.png" width="256"/>  <img src="/images/s_cf_6,2,3,1,_0.0,_.png" width="256"/> 

Approximation error plots:

 <img src="/images/delta_s_6,2,3,1,_0.0,_.png" width="256"/>  <img src="/images/delta_cf_6,2,3,1,_0.0,_.png" width="256"/> 

Spline was evaluated on a uniform Cartesian grid of size 101x101. Accuracy of the interpolation was measured by calculating the Root Mean Square Error (RMSE) and the Maximum Absolute Error (MAE). For this case
```RMSE```: 1.6E-03, ```MAE```: 1.1E-01, estimated value of the scaling parameter ```ε``` is 8, estimation of the Gram matrix condition number is 1.0E+11.

The second example is the function ```Ψ(x,y,z)=cos(π*x)*cos(y-0.5)*sin(π*(z-0.5))``` defined on ```Ω=[0,1]x[0,1]x[0,1]```.

<img src="/images/m_grid_6,1,50,_.png" width="256"/>  <img src="/images/m_grid_6,2,50,_.png" width="256"/>  <img src="/images/m_nodes_6,2,4.png" width="256"/> 

We'll construct an interpolating normal spline using function ```Ψ``` values sampled on the 1000 non-uniform random nodes distributed in the [0,1]x[0,1]x[0,1] cube.
```
    using NormalHermiteSplines

    nodes = get_3D_random_grid(9)       # generates 1000 non-uniform random nodes
    n_1 = size(nodes, 2)
    u = Vector{Float64}(undef, n_1)     # function values
    grid = get_3D_grid(50)              # creates uniform Cartesian grid of size 51x51x51
                                        # in [0, 1]x[0, 1]x[0, 1]
    for i = 1:n_1
        x = nodes[1,i]
        y = nodes[2,i]
        z = nodes[3,i]
        u[i] = cos(π*x)*cos(y-0.5)*sin(π*(z-0.5))
    end

    # Here spline is being constructed with RK_H2 kernel,
    # the 'scaling parameter' ε is defined explicitly.
    rk = RK_H2(5.0)
    #
    spline = interpolate(nodes, u, rk)
    σ = evaluate(spline, grid)

    σ1 = evaluate_one(spline, [0.8; 0.6; 0.8])
    #  ≈ -0.6512
    
    g1 = evaluate_gradient(spline, [0.8; 0.6; 0.8])
    #  ≈ -1.486
    #  ≈  0.065
    #  ≈ -1.486
```
The spline plots:

 <img src="/images/s_grid_6,false,2,4,2,_5.0,1,50,_.png" width="256"/>  <img src="/images/s_grid_6,false,2,4,2,_5.0,2,50,_.png" width="256"/> 

Approximation error plots:

 <img src="/images/delta_t_6,false,2,4,2,_5.0,1,50,_.png" width="256"/>  <img src="/images/delta_t_6,false,2,4,2,_5.0,2,50,_.png" width="256"/> 

Spline was evaluated on a uniform Cartesian grid of size 51x51x51. For this case ```RMSE```: 2.9E-03, ```MAE```: 5.1E-02, value of the scaling parameter ```ε``` is 5.0, estimation of the Gram matrix condition number is 1.0e+13.

Further examples are given in documentation.  

## Documentation

For more information and explanation see [Documentation](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Hermite-Splines/)

The normal splines method for one-dimensional function interpolation and linear ordinary differential and integral equations was proposed in [2]. Multivariate generalization of the normal splines method was developed for two-dimensional problem of low-range computerized tomography in [3] and applied for solving a mathematical economics problem in [4]. Further results were reported at the seminars and conferences [5,6,7].

## References:
1. [Halton sequence](https://en.wikipedia.org/wiki/Halton_sequence)
2. V. Gorbunov, The method of normal spline collocation. [USSR Computational Mathematics and Mathematical Physics, Vol. 29, No. 1, 1989](https://www.researchgate.net/publication/265357408_Method_of_normal_spline-collocation)
3. I. Kohanovsky, Normal Splines in Computing Tomography (in Russian). [Avtometriya, No. 2, 1995](https://www.iae.nsk.su/images/stories/5_Autometria/5_Archives/1995/2/84-89.pdf)
4. V. Gorbunov, I. Kohanovsky, K. Makedonsky, Normal splines in reconstruction of multi-dimensional dependencies. [Papers of WSEAS International Conference on Applied Mathematics, Numerical Analysis Symposium, Corfu, 2004](http://www.wseas.us/e-library/conferences/corfu2004/papers/488-312.pdf)
5. I. Kohanovsky, Multidimensional Normal Splines and Problem of Physical Field Approximation, International Conference on Fourier Analysis and its Applications, Kuwait, 1998.
6. I. Kohanovsky, Inequality-Constrained Multivariate Normal Splines with Some Applications in Finance. [27th GAMM-Seminar Leipzig on Approximation of Multiparametric functions, 2011](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)
7. V. Gorbunov, I. Kohanovsky, Heterogeneous Parallel Method for the Construction of Multi-dimensional Smoothing Splines. [ESCO 2014 4th European Seminar on Computing, 2014](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)













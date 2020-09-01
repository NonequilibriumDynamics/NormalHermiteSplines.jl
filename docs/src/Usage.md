# Example Usage

## 1D interpolation case

Let's interpolate function ``f(x)``

```math
f(x) =
    \begin{cases}
       0 \, , &  1 \le x \lt 6 \\
       1 \, , &  6 \le x \le 10 \\
       -x/5 + 3 \, , &  6 \le x \le 15 \\
       0 \, , &  15 \le x \le 20 
    \end{cases}
```
by values of the function in nodes ``\{1, 2, 3, ..., 20\}`` (case A) and by values of the function and values of its first derivatives in the same nodes (case B).

```@meta
DocTestSetup = quote
    using Random
    using NormalHermiteSplines
end
```

A)

```@example A
    using NormalHermiteSplines

    x = collect(1.0:1.0:20)       # function nodes
    u = x.*0.0                    # function values in nodes
    for i in 6:10
        u[i] = 1.0
    end
    for i in 11:14
        u[i] = -0.2 * i + 3.0
    end

    # Build a differentiable spline by values of function in nodes
    # (a spline built with RK_H0 kernel is a continuous function,
    #  a spline built with RK_H1 kernel is a continuously differentiable function,
    #  a spline built with RK_H2 kernel is a twice continuously differentiable function).
    # Here value of the 'scaling parameter' ε is estimated in the interpolate procedure.
    spline = prepare(x, RK_H1())
    
    # An estimation of the Gram matrix condition number
    cond = get_cond(spline)
```

```@example A
    # A value of the 'scaling parameter' of Bessel Potential space
    # the normal spline was built in.
    ε = get_epsilon(spline)
```

```@example A
    spline = construct(spline, u)

    p = collect(1.0:0.2:20)        # evaluation points
    σ = evaluate(spline, p)
    σ = nothing                    # hide
```

![Example 1A](images/1d-usage/example-1a.svg)

```@example A
    dσ = similar(p)
    for i=1:length(p)
        dσ[i] = evaluate_derivative(spline, p[i])
    end
```
![Example 2A](images/1d-usage/example-1a-der.svg)


Evaluate the spline at some points:

```@example A
    p = [3.1, 8.1, 12.1, 18.1]
    σ = evaluate(spline, p)
```

Evaluate the spline derivatives at the same points:

```@example A
    dσ = similar(p)
    for i=1:length(p)
        dσ[i] = evaluate_derivative(spline, p[i])
    end
    dσ
```

Construct spline by different function values in nodes and evaluate new spline at the same evaluation points:

```@example A
    u2 = 2.0 .* u
    spline = construct(spline, u2)
    σ = evaluate(spline, p)
```

B)

```@example B
    using NormalHermiteSplines

    x = collect(1.0:1.0:20)       # function nodes
    u = x.*0.0                    # function values in nodes
    for i in 6:10
        u[i] = 1.0
    end
    for i in 11:14
        u[i] = -0.2 * i + 3.0
    end

    s = x                         # function first derivative nodes
    v = x.*0.0                    # function first derivative values
    for i in 11:14
        v[i] = -0.2
    end

    # Build a differentiable spline by values of function,
    # and values of its first derivatives in nodes
    # (a spline built with RK_H0 kernel is a continuous function,
    #  a spline built with RK_H1 kernel is a continuously differentiable function,
    #  a spline built with RK_H2 kernel is a twice continuously differentiable function).
    # Here value of the 'scaling parameter' ε is estimated in the interpolate procedure.
    spline = interpolate(x, u, s, v, RK_H1())

    p = collect(1.0:0.2:20)      # evaluation points
    σ = evaluate(spline, p)
    σ = nothing                    # hide
```

![Example 1B](images/1d-usage/example-1b.svg)

```@example B
    dσ = similar(p)
    for i=1:length(p)
        dσ[i] = evaluate_derivative(spline, p[i])
    end
```
![Example 2B](images/1d-usage/example-1b-der.svg)

Evaluate the spline at some points:

```@example B
    p = [3.1, 8.1, 12.1, 18.1]
    σ = evaluate(spline, p)
```

Evaluate the spline derivatives at the same points:

```@example B
    dσ = similar(p)
    for i=1:length(p)
        dσ[i] = evaluate_derivative(spline, p[i])
    end
    dσ
```

## 2D interpolation case

Let's interpolate function ``\phi (x,y)  = \frac{2}{3}cos(10x)sin(10y) + \frac{1}{3}sin(10xy)`` defined on domain ``\Omega = [0,1]^2``.

```@raw html
<img src="../images/2d-usage/m-cf-32.png" width="256"/>
``` ```@raw html
<img src="../images/2d-usage/m-t-32.png" width="256"/>
```  ```@raw html
<img src="../images/2d-usage/m-grid-32,3.png" width="197"/>
```
We built a spline using function ``\phi`` values sampled on 200 pseudo-random points uniformly distributed on ``\Omega`` (case A).

Spline plot                                                                        Approximation error plots

```@raw html
<img src="../images/2d-usage/s-cf-32,32,3,1,0.0,-.png" width="256"/>
``` ```@raw html
<img src="../images/2d-usage/delta-cf-32,32,3,1,0.0,-.png" width="256"/>
```  ```@raw html
<img src="../images/2d-usage/delta-s-32,32,3,1,0.0,-.png" width="256"/>
```
and 

using function ``\phi`` values sampled on 200 pseudo-random points uniformly distributed on ``\Omega`` and 80 values of function ``\phi`` gradient given at the border of ``\Omega`` (case B):

```@raw html
<img src="../images/2d-usage/m-grid-33,3.png" width="197"/> 
```
here function ``\phi`` gradient locations are marked with 'red'. 

Spline plot                                                                         Approximation error plots

```@raw html
<img src="../images/2d-usage/s-cf-33,33,3,1,0.0,-.png" width="256"/>
``` ```@raw html
<img src="../images/2d-usage/delta-cf-33,33,3,1,0.0,-.png" width="256"/>
```  ```@raw html
<img src="../images/2d-usage/delta-s-33,33,3,1,0.0,-.png" width="256"/>
```
ccc

ddd

```@example 2A
    using Random
    using NormalHermiteSplines

    # generating 200 uniform random nodes
    m = 200
    nodes = Matrix{Float64}(undef, 2, m)
    rng = MersenneTwister(0);
    rnd = rand(rng, Float64, (2, m))
    for i = 1:m
        nodes[1, i] = rnd[1, i]
        nodes[2, i] = rnd[2, i]
    end

    # creating the uniform Cartesian grid of size 51x51 on [0, 1]x[0, 1]
    t = 50
    x = collect(range(0.0, 1.0; step = 1.0/t))
    y = collect(range(0.0, 1.0; step = 1.0/t))
    t1 +=  1
    grid = Matrix{Float64}(undef, 2, t1^2)
    for i = 1:t1
        for j = 1:t1
            r = (i - 1) * t1 + j
            grid[1, r] = x[i]
            grid[2, r] = y[j]
        end
    end

    n_1 = size(nodes, 2)
    u = Vector{Float64}(undef, n_1)     # function values
    for i = 1:n_1
        x = nodes[1,i]
        y = nodes[2,i]
        u[i] = (2.0*cos(10.0*x)*sin(10.0*y) + sin(10.0*x*y))/3.0
    end

    # Here spline is being constructed with RK_H1 kernel,
    # the 'scaling parameter' ε is defined explicitly.
    rk = RK_H1()
    #
    spline = interpolate(nodes, u, rk)
    cond = get_cond(spline)
```

```@example 2A
    # A value of the 'scaling parameter' of Bessel Potential space
    # the normal spline was built in.
    ε = get_epsilon(spline)
```


```@example 2A

    σ = evaluate(spline, grid)
    σ = nothing
```


Value of function ``\phi`` in evaluation point ``p = [0.5; 0.5]``
```@example 2A
    p = [0.5; 0.5]
    f1 = (2.0*cos(10.0*p[1])*sin(10.0*p[2]) + sin(10.0*p[1]*p[2]))/3.0
```

Value of spline in that evaluation point:
```@example 2A
    σ1 = evaluate_one(spline, p)
```

```@example 2A
    g1 = evaluate_gradient(spline, p)
    #  ≈ -1.486
    #  ≈  0.065
    #  ≈ -1.486
```









## Choice of the scaling parameter
aa

## Q & A
ss

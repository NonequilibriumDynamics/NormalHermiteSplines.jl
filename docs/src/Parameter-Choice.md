# Choice of the scaling parameter

Approximating properties of the normal spline are getting better with the smaller value of the scaling parameter $\varepsilon$ (this parameter also is known as the "shape" parameter in RBF literature), and if the value of this parameter is small enough then normal spline become similar to Duchon's $D^m -$spline [1]. See also
[Relation to Polyharmonic Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Relation-to-Polyharmonic-Splines/).

However with decreasing value of the parameter $\varepsilon$ the condition number of the corresponding problem Gram matrix is increasing and the problem becomes numerically unstable. The Gram matrix of the interpolating problem even can lost its positive deiniteness property if $\varepsilon$ is small.

## TODO - rewrite it! 
 Also, it was pointed out [2] that RBF interpolation with small value of the scaling (shape) parameter $\varepsilon$ may cause the Runge's phenomenon, but very small values of $\varepsilon$ have never been used here and Runge-type oscillations have not been observed. 

Therefore, when choosing the value of the $\varepsilon$, a compromise is needed. In practice, it is necessary to choose such value of the scaling parameter that condition number of the problem Gram matrix is a small enough number. As a rule, the heuristic algorithm implemented within the interpolation procedure produces a good estimation of the scaling parameter value (this algorithm applies if the value of the scaling parameter was not provided explicitly in creation of the reproducing kernel object.)

The following API functions could be useful for selecting a suitable value of the scaling parameter $\varepsilon$:

- ```get_cond```
- ```get_epsilon```
- ```estimate_epsilon```  
- ```estimate_accuracy```

Let's consider an example with interpolation of function 

```math
\phi (x,y)  = \frac{2}{3}cos(10x)sin(10y) + \frac{1}{3}sin(10xy)
```
sampled on set of 200 pseudo-random nodes uniformly distributed on unit square ``\Omega = [0,1]^2``.

```
    using NormalHermiteSplines
    ....
    ....
    # Here spline is being constructed with RK_H1 kernel,
    # the value of the 'scaling parameter' ε is estimated
    # in the interpolate procedure.
    rk = RK_H1()
    #
    spline = interpolate(nodes, u, rk)
    σ = evaluate(spline, grid)
```
(the complete code example is here: [Example usage](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Usage/#D-interpolation-case-2/))

Let's get values of scaling parameter, estimation of the Gram matrix condition number (algorithm is taken from [3]) and assessed value of the interpolation quality (value of the maximum of relative residual error calculated
using data of the function value interpolation nodes).
```
    ε = get_epsilon(spline)
```

```
    cond = get_cond(spline)
```


**References**

[1] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Springer, Berlin, Vol. 571, 1977.

[2] B. Fornberg, J. Zuev, The Runge phenomenon and spatially variable shape parameters in RBF interpolation,
Comput. Math. Appl., Vol.54, No.3, 2007.

[3] C. Brás, W. Hager, J. Júdice, An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP 20, 2012.

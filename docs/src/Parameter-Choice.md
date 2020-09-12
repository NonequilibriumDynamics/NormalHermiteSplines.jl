# Choice of the scaling parameter

Approximating properties of the normal spline are getting better with the smaller value of the scaling parameter $\varepsilon$ (this parameter also is known as the "shape" parameter in RBF literature), and if the value of this parameter is small enough then normal spline become similar to Duchon's $D^m -$spline [2]. See also
[Relation to Polyharmonic Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Relation-to-Polyharmonic-Splines/).

However with decreasing value of the parameter $\varepsilon$ the condition number of the corresponding problem Gram matrix is increasing and the problem becomes numerically unstable (see "uncertainty principle" of Schaback [4]). The Gram matrix of the interpolating problem even can lost its positive deiniteness property if $\varepsilon$ is small. Also, as it was pointed out in [3] the RBF interpolation with small value of the "shape" (scaling) parameter may cause the Runge phenomenon i.e. undesirable interpolant oscillations which are most likely observed at the domain border. 

Therefore, when choosing the value of the $\varepsilon$, a compromise is needed. In practice, it is necessary to choose such value of the scaling parameter that condition number of the problem Gram matrix is a small enough number. As a rule, the heuristic algorithm implemented within the interpolation procedure produces a good estimation of the scaling parameter value (this algorithm applies if the value of the scaling parameter was not provided explicitly in creation of the reproducing kernel object.)

The following API functions could be useful for selecting a suitable value of the scaling parameter $\varepsilon$:

- ```get_cond```
- ```get_epsilon```
- ```estimate_accuracy```
- ```estimate_epsilon```  

Let's consider an example with interpolation of function 

```math
\phi (x,y)  = \frac{2}{3}cos(10x)sin(10y) + \frac{1}{3}sin(10xy)
```
sampled on set of 200 pseudo-random nodes uniformly distributed on unit square ``\Omega = [0,1]^2``.

```
    using NormalHermiteSplines
    ....
    spline = interpolate(nodes, u, rk)
    σ = evaluate(spline, grid)
    ε = get_epsilon(spline)
    cond = get_cond(spline)
    valid_digits = estimate_accuracy(spline)
    ....
```
(the complete code example can be found here: [Example usage](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Usage/#D-interpolation-case-2/))

The results of this function interpolation with reproducing kernel ```RK_H0``` are displayed in Table I, results of interpolation with reproducing kernel ```RK_H1``` are displayed in Table II and results of interpolation received with reproducing kernel ```RK_H2``` - in Table III.
Condition number (```κ```) of the interpolation problem Gram matrix were estimated by procedure described in [1]. The root mean squared error (```RMSE```), the maximum absolute error (```MAE```) and estimated accuracy (```EAC```) are defined as follows:







**References**

[1] C. Brás, W. Hager, J. Júdice, An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP 20, 2012.

[2] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Springer, Berlin, Vol. 571, 1977.

[3] B. Fornberg, J. Zuev, The Runge phenomenon and spatially variable shape parameters in RBF interpolation,
Comput. Math. Appl., Vol.54, No.3, 2007.

[4]] R. Schaback, Error estimates and condition numbers for radial basis functions interpolation, Adv. in Comput. Math. 3, 1995.
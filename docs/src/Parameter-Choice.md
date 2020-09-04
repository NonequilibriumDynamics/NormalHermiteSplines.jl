# Choice of the scaling parameter

Approximating properties of the normal spline are getting better with the smaller value of the scaling parameter $\varepsilon$, and if the value of this parameter is small enough then normal spline become similar to Duchon's $D^m -$spline [1]. See also
[Relation to Polyharmonic Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Relation-to-Polyharmonic-Splines/).

However with decreasing value of the parameter $\varepsilon$ the condition number of the corresponding problem Gram matrix is increasing and the problem becomes numerically unstable. The Gram matrix of the interpolating problem can lost its positive deiniteness property if $\varepsilon$ is small. Therefore, when choosing the value of the $\varepsilon$, a compromise is needed. In practice, it is necessary to choose such value of the $\varepsilon$ that condition number of Gram matrix is small enough. 

The following API functions could be useful for selecting a suitable value of the scaling parameter:

- ```assess_interpolation```
- ```get_cond```
- ```get_epsilon```
- ```estimate_epsilon```  

Let's consider an example with interpolating a test function 

```math
\phi (x,y)  = \frac{2}{3}cos(10x)sin(10y) + \frac{1}{3}sin(10xy)
```
 sampled on set of 200 pseudo-random nodes uniformly distributed on unit square ``\Omega = [0,1]^2``.




**References**

[1] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Springer, Berlin, Vol. 571, 1977.

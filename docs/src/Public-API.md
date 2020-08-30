# Public API

## API Summary

| Function                            | Description                                                                                        |
|:----------------------------------- |:-------------------------------------------------------------------------------------------------- |
|```prepare```                        |Prepare the spline by constructing and factoring a Gram matrix of the interpolation problem.        |
|```construct```                      |Construct the spline by calculating its coefficients.                                               |
|```interpolate```                    |Prepare and construct the spline.                                                                   |
|```evaluate```                       |Evaluate the spline value at the required locations                                                 |
|```evaluate_one```                   |Evaluate the spline value at the required location                                                  |
|```evaluate_gradient```              |Evaluate gradient of the spline at the required location.                                           |
|```evaluate_derivative```            |Evaluate the 1D spline derivative at the required location.                                         |
|```assess_interpolation```           |Assess interpolation result.                                                                        |
|```get_cond```                       |Get an estimation of the Gram matrix condition number.                                              |
|```get_epsilon```                    |Get the 'scaling parameter' of Bessel Potential space the normal spline was built in.               |
|```estimate_epsilon```               |Get an estimation of the 'scaling parameter' of Bessel Potential space the spline being built in.   |

## Functions
```@docs
prepare
construct
interpolate
evaluate
evaluate_one
evaluate_gradient
evaluate_derivative
assess_interpolation
get_cond
get_epsilon
estimate_epsilon
```

## Types

### Bessel Potential space Reproducing Kernels

```@docs
RK_H0
RK_H1
RK_H2
```

### NormalSpline structure

```@docs
NormalSpline
```

## Index
```@index
Order = [:function, :type]
```

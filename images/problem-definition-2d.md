This Julia package implements the normal splines method for solving following interpolation problem:

*Problem:* Given points ``\{(x_i, y_i), (x_i, y_i) \in R^2\}_{i=1}^{n}`` find a  function ``\psi`` such that

```math
\tag{1}
\begin{aligned}

 a_i \frac{\partial^2{\psi}}{\partial{x}^2 }(x_i, y_i) + b_i \frac{\partial^2{\psi}}{\partial{y}^2 }(x_i, y_i) + c_i \frac{\partial^2{\psi}}{\partial{x} \partial{y} }(x_i, y_i) + d_i \frac{\partial{\psi}}{\partial{x}}(x_i, y_i) + e_i \frac{\partial{\psi}}{\partial{y}}(x_i, y_i) + f_i \psi (x_i, y_i) = u_i \, ,
\\
 \qquad \qquad 
a_i, \, b_i, \, c_i, \, d_i, \, e_i, \, f_i, \, u_i \, \in R^1 \, , \ \ (x_i, y_i) \in R^2 \, , \ \ i = 1, 2, \dots, n \, .
\end{aligned}
``` 
We assume that function ``\psi`` is an element of the Bessel Potential space ``H^s_\varepsilon (R^2)`` which is defined as:

```math
   H^s_\varepsilon (R^2) = \left\{ \varphi | \varphi \in S' (R^2) , \ 
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^2) \right\} , \quad
  \varepsilon \gt 0 , \ \ s = r + \frac{3}{2} \, , \quad r = 2,3,\dots \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^2)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. The parameter ``\varepsilon`` can be considered as a "scaling parameter", it allows to control approximation properties of the normal spline which usually are getting better with smaller values of ``\varepsilon``, also it can be used to reduce the ill-conditioness of the related computational problem (in traditional theory ``\varepsilon = 1``).

The Bessel Potential space ``H^s_\varepsilon (R^2)`` is a  Reproducing kernel Hilbert space, an element ``\psi`` of that space can be treated as a ``r``-times continuously differentiable function.

The normal splines method consists in finding a solution of system (1) having minimal norm in Hilbert space ``H^s_\varepsilon (R^2)``, thus the interpolation normal spline ``\sigma`` is defined as follows:

```math
\tag{2}
   \sigma = {\rm arg\,min}\{  \| \psi - z \|^2 : (1), \forall \psi \in H^s_\varepsilon (R^2) \}  \, , \quad z \in H^s_\varepsilon (R^2) \, ,
```
where ``z`` is a 'prototype` element (a-priori guess of the solution).

Normal splines method is based on the following functional analysis results:

* Bessel Potential space embedding theorem
* The Riesz representation theorem for Hilbert spaces
* Reproducing kernel properties 

Using these results it is possible to reduce the task (2) to solving a system of linear equations with symmetric positive definite Gram matrix.

Detailed explanation is given in the package documentation.

```@meta
Author = "Igor Kohanovsky"
```

# NormalHermiteSplines.jl

*Multivariate Normal Hermite Splines in Julia*

`NormalHermiteSplines.jl` implements the normal splines method for solving following interpolation problem:

*Problem:*   Given points ``\{p_i, p_i \in R^n\}_{i=1}^{n_1}``, ``\{s_j, s_j \in R^n\}_{j=1}^{n_2}`` and a set of unit vectors ``\{e_j, e_j \in R^n\}_{j=1}^{n_2}`` find a function ``f`` such that

```math
\tag{1}
\begin{aligned}
& f(p_i) =  u_i \, , \quad  i = 1, 2, \dots, n_1 \, ,
\\  
& \frac{ \partial{f} }{ \partial{e_j} }(s_j) =  v_j \, , \quad  j = 1, 2, \dots, n_2 \, ,
\\
& n_1 \gt 0 \, ,  \ \  n_2 \ge 0 \, .
\end{aligned}
```
where ``\frac{ \partial{f} }{ \partial{e_j} }(s_j) = \nabla f(s_j) \cdot e_j = \sum _{k=1}^{n}  \frac{ \partial{f} }{ \partial{x_k} } (s_j) e_{jk}`` is a directional derivative of ``f`` at the point ``s_j`` in the direction of ``e_j``.

We assume that function ``f`` is an element of the Bessel potential space ``H^s_\varepsilon (R^n)`` which is defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \ \ s = n/2 + 1/2 + r \, , \quad r = 1,2,\dots \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. The parameter ``\varepsilon`` can be considered as a "scaling parameter", it allows to control approximation properties of the normal spline which usually are getting better with smaller values of ``\varepsilon``, also it can be used to reduce the ill-conditioness of the related computational problem (in traditional theory ``\varepsilon = 1``).

The Bessel potential space ``H^s_\varepsilon (R^n)`` is a Hilbert space, an element ``f`` from that space can be treated as a ``r``-times continuously differentiable function.

The normal splines method consists in finding a solution of system (1) having minimal norm in Hilbert space ``H^s_\varepsilon (R^n) ,`` thus the interpolation normal spline ``\sigma`` is defined as follows:

```math
\tag{2}
   \sigma = {\rm arg\,min}\{  \| f \|^2 : (1), \forall f \in H^s_\varepsilon (R^n) \} \, .
```

Normal splines method is based on the following functional analysis results:

* Bessel potential space embedding theorem
* The Riesz representation theorem for Hilbert spaces
* Reproducing kernel properties

Using these results it is possible to reduce the task (2) to solving a system of linear equations with symmetric positive definite Gram matrix.

Detailed explanation is given in [Multivariate Normal Hermite Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Hermite-Splines/)

The normal splines method for one-dimensional function interpolation and linear ordinary differential and integral equations was proposed in [1]. Multivariate generalization of the normal splines method was developed for two-dimensional problem of low-range computerized tomography in [2] and applied for solving a mathematical economics problem in [3]. Further results were reported at the seminars and conferences [4,5,6].

#### References:
1. V. Gorbunov, The method of normal spline collocation. [USSR Computational Mathematics and Mathematical Physics, Vol. 29, No. 1, 1989](https://www.researchgate.net/publication/265357408_Method_of_normal_spline-collocation)
2. I. Kohanovsky, Normal Splines in Computing Tomography (in Russian). [Avtometriya, No. 2, 1995](https://www.iae.nsk.su/images/stories/5_Autometria/5_Archives/1995/2/84-89.pdf)
3. V. Gorbunov, I. Kohanovsky, K. Makedonsky, Normal splines in reconstruction of multi-dimensional dependencies. [Papers of WSEAS International Conference on Applied Mathematics, Numerical Analysis Symposium, Corfu, 2004](http://www.wseas.us/e-library/conferences/corfu2004/papers/488-312.pdf)
4. I. Kohanovsky, Multidimensional Normal Splines and Problem of Physical Field Approximation, International Conference on Fourier Analysis and its Applications, Kuwait, 1998.
5. I. Kohanovsky, Inequality-Constrained Multivariate Normal Splines with Some Applications in Finance. [27th GAMM-Seminar Leipzig on Approximation of Multiparametric functions, 2011](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)
6. V. Gorbunov, I. Kohanovsky, Heterogeneous Parallel Method for the Construction of Multi-dimensional Smoothing Splines. [ESCO 2014 4th European Seminar on Computing, 2014](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)

## Contents

```@contents
Pages = [
      "index.md",
      "Public-API.md",
      "Usage.md",
      "Parameter-Choice.md",
      "Numerical-Tests.md",
      "Tests-with-real-data.md",
      "Interpolating-Normal-Splines.md",
      "Reproducing-Kernel-of-Bessel-Potential-space.md",
      "Relation-to-Polyharmonic-Splines.md",
      "Updating-Cholesky-Factorization.md"
]
Depth = 3
```

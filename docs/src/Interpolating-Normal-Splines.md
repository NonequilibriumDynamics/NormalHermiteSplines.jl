# Interpolating Normal Splines

Consider the following interpolation problem:

*Problem:* Given points ``\{p_i, p_i \in R^n\}_{i=1}^{n_1}``, ``\{s_j, s_j \in R^n\}_{j=1}^{n_2}`` and a set of unit vectors ``\{e_j, e_j \in R^n\}_{j=1}^{n_2}`` find a function ``f`` such that

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

We assume that function ``f`` is an element of the Bessel Potential space ``H^s_\varepsilon (R^n)`` which is defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \ \ s = n/2 + 1/2 + r \, , \quad r = 1,2,\dots \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. The parameter ``\varepsilon`` can be considered as a "scaling parameter", it allows to control approximation properties of the normal spline which usually are getting better with smaller values of ``\varepsilon``, also it can be used to reduce the ill-conditioness of the related computational problem (in traditional theory ``\varepsilon = 1``).

Theoretical properties of spaces ``H^s_\varepsilon`` at ``\varepsilon \gt 0`` are identical — they are Hilbert spaces with inner product

```math
\langle \varphi , \psi \rangle =
\int ( \varepsilon ^2  + | \xi |^2 )^s
\mathcal F [\varphi ] \overline{\mathcal F [\psi ] } \, d \xi
```
and  norm

```math
\| \varphi \| = \left( \langle \varphi , \varphi \rangle \right)^{1/2} =
\| (  \varepsilon ^2 + | \xi |^2 )^{s/2} {\mathcal F} [\varphi ] \|_{L_2} \ .
```
It is easy to see that all these norms are equivalent. It means that space ``H^s_\varepsilon (R^n)`` is equivalent to ``H^s (R^n) =  H^s_1 (R^n)`` ([3], [5]).

 Hilbert space ``H^s_\varepsilon (R^n)`` is continuously embedded in Hölder space ``C_b^r(R^n)`` ([3],[17],[20]) of functions continuous and bounded with their first ``r`` derivatives, it means function ``f`` can be treated as an element of function class ``C^r(R^n)`` of functions continuous with their first ``r`` derivatives. Therefore functionals ``F_i`` and ``F'_j``

```math
\begin{aligned}
  & F_i(\varphi) = \varphi (p_i) \, , \ \ F'_j(\varphi) = \frac{ \partial{\varphi} }{ \partial{e_j} }(s_j) \, , \ \ \forall \varphi \in H^s_\varepsilon (R^n) \, ,
   \quad
   p_i, \, s_j \in R^n \, ,
   \\
  & i = 1, 2, \dots, n_1 \, , \ \   j = 1, 2, \dots, n_2 \, ,
\end{aligned}
```
are linear continuous functionals in ``H^s_\varepsilon``.

We also assume that all points ``\{p_i\}`` are different and in a case when among points ``\{s_j\}`` there are coincident ones, we stipulate that the corresponding unit vectors defining the directions of the directional derivatives at such points are linearly independent. Note that some points ``\{p_i\}`` may coincide with some ``\{s_j\}``. Under these restrictions all functionals $F_i, \, F'_j$ are linearly independent.

 In accordance with Riesz representation theorem [1] these linear continuous functionals can be represented in the form of inner product of some elements ``h_i, h'_j \in H^s_\varepsilon`` and ``\varphi \in H^s_\varepsilon``, for any ``\varphi \in H^s_\varepsilon``:

```math
 F_i(\varphi) = {\langle h_i,  \varphi \rangle} \, ,  \quad  F'_j(\varphi) = {\langle h'_j,  \varphi \rangle} \, ,  \quad  \forall \varphi \in H^s_\varepsilon \, ,
\\  i = 1, 2, \dots, n_1 \, , \ \  j = 1, 2, \dots, n_2 \, .
```
Elements ``h_i`` and ``h'_j`` are continuously differentiable functions. Thereby the original system of constraints (1) can be written in form:

```math
\tag{2}
\begin{aligned}
& f(p_i) = F_i(f) = {\langle h_i,  f \rangle}  = u_i \, ,
\\
& \frac{ \partial{f} }{ \partial{e_j} }(s_j) = F'_j(f) = {\langle h'_j,  f \rangle} = v_j \, ,  
 \\
 &  h_i, \, h'_j, \, f \in H^s_\varepsilon \, ,
 \quad
  i = 1, 2, \dots, n_1 \, , \ \  j = 1, 2, \dots, n_2 \, .
 \end{aligned}
```
here all functions ``h_i`` and ``h'_j`` are linear independent and system of constrains (2) defines a nonempty convex and closed set (as an intersection of hyper-planes) in the Hilbert space ``H^s_\varepsilon``.

Problem of reconstruction of function ``f`` satisfying system of constraints (2) is undetermined. We reformulate it as a problem of finding solution of this system of constraints that has minimal norm:

```math
\tag{3}
   \sigma = {\rm arg\,min}\{  \| f - z \|^2 : (2), z \in H^s_\varepsilon , \forall f \in H^s_\varepsilon \} \, ,
```
where ``z \in H^s_\varepsilon`` is a "prototype" function. Solution of this problem exists and it is unique ([6], [16]) as a projection of element ``z`` on the nonempty convex closed set in Hilbert space ``H^s_\varepsilon``. Element ``\sigma`` is an interpolating normal spline.

In accordance with generalized Lagrange method ([13], [16]) solution of the problem (3) can be presented as:

```math
\tag{4}
\sigma =  z + \sum _{i=1}^{n_1} \mu_i  h_i  + \sum _{j=1}^{n_2} \mu'_j  h'_j  \, ,
```
where coefficients ``\mu_i`` and ``\mu'_j`` are defined by system of linear equations

```math
\tag{5}
    \sum _{l=1}^{n_1} g_{il} \mu_l + \sum _{j=1}^{n_2} g'_{ij} \mu'_j   =  u_i - {\langle h_i,  z \rangle} \, , \quad 1 \le i \le n_1 \, , \\      \sum _{i=1}^{n_1} g'_{ij} \mu_i + \sum _{m=1}^{n_2} g''_{jm} \mu'_m\ =  v_j - {\langle h'_j,  z \rangle}  \, , \quad 1 \le j \le n_2 \, ,
```
Matrix of system (5) is the positive definite symmetric Gram matrix of the set of linearly independent elements  ``\{h_i\}, \{h'_j\}`` and coefficients ``g_{il}, g'_{ij}, g''_{jm}`` are defined as follows:

```math
\tag{6}
   g_{il} = {\langle h_i,  h_l \rangle} \, , \ \ g'_{ij} = {\langle h_i,  h'_j \rangle} \, , \ \  g''_{jm} = {\langle h'_j,  h'_m \rangle} \, .
```

Space ``H^s_\varepsilon (R^n)`` is a reproducing kernel Hilbert space ([5],[18],[19],[20]). We denote its reproducing kernel as ``V(\eta, \xi)``.

Recall the definition of the reproducing kernel ([4], [7]). The reproducing kernel of space ``H^s_\varepsilon`` is a such function ``V(\eta, \xi)`` that

* for every ``\xi \in R^n``, ``\ V(\eta, \xi)`` as function of ``\eta`` belongs to ``H^s_\varepsilon``
* for every ``\xi \in R^n`` and every function ``\varphi \in H^s_\varepsilon``

```math
\tag{7}
\varphi(\xi) = {\langle V(\eta, \xi),  \varphi(\eta) \rangle}
```
Reproducing kernel is a symmetric function:

```math
V(\eta, \xi) = V(\xi, \eta) \, ,
```
also in the considered case (``s = n/2 + 1/2 + r, \, r \ge 1``) it is a continuously differentiable function. Differentiating the identity (7) allows to get the identities for derivatives:

```math
\tag{8}
\frac {\partial \varphi(\xi)}{\partial \xi_k} = {\left \langle \frac{\partial V(\cdot, \xi)} {\partial \xi_k}, \varphi \right \rangle}
```
which holds for any ``\varphi \in H^s_\varepsilon`` and ``\xi \in R^n``, it means that function ``\frac{\partial {V(\cdot , \xi)} }{\partial{\xi_k}}`` represents a point-wise functional defined as value of function ``\frac{ \partial {\varphi (\cdot)} }{\partial{\xi_k}}`` at the point ``\xi``.

Now it is possible to express functions ``h_i`` and ``h'_j`` via the reproducing kernel ``V``. Comparing (2) with (7) and (8) we receive:

```math
\tag{9}
\begin{aligned}
&  h_i (\eta) =  V(\eta, p_i) \, ,  \qquad \qquad \qquad \qquad \qquad \ \ i = 1, 2, \dots, n_1 \, \\  
&  h'_j (\eta) =  \frac{\partial V(\eta, s_j)}{\partial e_j} =  \sum_{k=1}^n  \frac{ \partial {V(\eta, s_j)} }{\partial{\xi_k}} e_{jk}  \, , \quad  j = 1, 2, \dots, n_2 \ .  
\end{aligned}
```
The coefficients (6) of the Gram matrix can be presented as ([7], [8], [10]):

```math
\tag{10}
\begin{aligned}
   & g_{il} = {\langle h_i,  h_l \rangle} = {\langle V(\cdot, p_i),  V(\cdot, p_l) \rangle} = V(p_i, p_l) \, ,
   \\
    & g'_{ij} = {\langle h_i,  h'_j \rangle} = {\left \langle V(\cdot, p_i), \frac{\partial V(\cdot, s_j)}{\partial e_j} \right \rangle} =  \frac{\partial V(p_i, s_j)}{\partial e_j} =
    \\
    & \qquad \qquad \qquad =  \sum_{k=1}^n  \frac{ \partial {V(p_i, s_j)} }{\partial{\xi_k}} e_{jk}  \ .
\end{aligned}
```
With the help of (7) and (10), we can also calculate ``g''_{jm}`` ([8], [10]):

```math
\tag{11}
\begin{aligned}
 g''_{jm} = {\langle h'_j,  h'_m \rangle} \, & = \, {\left \langle \frac{\partial V(\cdot, s_j)}{\partial e_j}, \frac{\partial V(\cdot, s_m)}{\partial e_m} \right \rangle} \,  = \, \frac {\partial^2 V(s_j, s_m)} {\partial e_j \partial e_m} =
 \\
 &  =  \sum_{r=1}^n \sum_{k=1}^n  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_k}} e_{jk} e_{mr} \ .
\end{aligned}
```
Further

```math
\tag{12}
\begin{aligned}
  & {\langle h_i,  z \rangle}  =  {\langle V(\cdot, p_i),  z \rangle} = z(p_i)  \, ,
  \\
  & {\langle h'_j,  z \rangle}  = \frac{\partial z(s_j)}{\partial e_j} = \sum_{k=1}^n  \frac{ \partial {z(s_j)} }{\partial{x_k}} e_{jk} \, .
\end{aligned}
```

Here normal hermite splines will be constructed in Bessel Potential spaces ``H^{s_1}_\varepsilon (R^n) , \, s_1 = n/2 + 3/2`` and ``H^{s_2}_\varepsilon (R^n) , \, s_2 = n/2 + 5/2``. Elements of space ``H^{s_1}`` can be treated as  continuously differentiable functions and elements of space ``H^{s_2}`` can be treated as twice continuously differentiable functions. Note, the spline is infinitely differentiable everywhere in ``R^n`` excepting the nodes ``p_i`` and ``s_j``.

Reproducing kernel of Bessel Potential space was presented in [5] and its simplified form was given in [14], [18], [19], [20]. For space ``H^{s}_\varepsilon (R^n), \, s = n/2 + 1/2 + r, \, r \ge 0`` it can be written as:

```math
 V(\eta , \xi) = \exp (-\varepsilon |\eta - \xi|) \,
     \sum_{k=0}^{r} \frac{(r+k)!}{2^k k! (r-k)!} (\varepsilon |\eta - \xi|)^{r-k} \ ,
```
``(``a constant multiplier is omitted here.``)``

This reproducing kernel is known as the Matérn kernel [23].

Therefore for space ``H^{s_1}_\varepsilon (R^n)`` with accuracy to constant multiplier we get:

```math
\tag{13}
V(\eta, \xi) =  \exp (-\varepsilon |\eta - \xi|)
             (1 + \varepsilon |\eta - \xi|) \, .
```
and for space ``H^{s_2}_\varepsilon (R^n)``:
```math
\tag{14}
V(\eta, \xi) =  \exp (-\varepsilon |\eta - \xi|)
             (3 + 3\varepsilon |\eta - \xi| + \varepsilon ^2 |\eta - \xi| ^2 ) \, .
```
Let's write down expressions of ``h_i, h'_j, g_{il}, g'_{ij}, g''_{jm}`` for space ``H^{s_1}_\varepsilon (R^n)``:

```math
\tag{15}
\begin{aligned}
& h_i (\eta) =  \exp (-\varepsilon |\eta - p_i |) (1 + \varepsilon |\eta - p_i|) \, , \qquad  \quad \  i = 1, 2, \dots, n_1 \, , \\
& h'_j (\eta) = \varepsilon^2 \exp (-\varepsilon | \eta - s_j | ) \sum _{k=1}^n (\eta_k - s_{jk}) e_{jk} \, , \quad  j = 1, 2, \dots, n_2 \, , \\
& g_{il}= \exp (-\varepsilon | p_i - p_l | )(1 + \varepsilon | p_i - p_l | ) \, ,    \quad  \ \  i = 1, 2, \dots, n_1 \, , \ \ l = 1, 2, \dots, n_1 \, , \\
& g'_{ij} = \varepsilon^2 \exp (-\varepsilon |p_i - s_j | ) \sum _{k=1}^n (p_{ik} - s_{jk}) e_{jk}  \, , \ \ i = 1, 2, \dots, n_1 \, , \ \ j = 1, 2, \dots, n_2 \, ,
\\
& g''_{jm} =  \sum_{r=1}^n \sum_{k=1}^n  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_k}} e_{jk} e_{mr}
\\
& \quad \qquad j \ne m \, , \quad j = 1, 2, \dots, n_2 \, , \ \ m = 1, 2, \dots, n_2 \, ,
\\
& \text{where}
\\
&  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_r}} = \varepsilon^2 \exp (-\varepsilon | s_j - s_m |) \left (1 - \varepsilon \frac {(s_{jr} - s_{mr})^2}{| s_j - s_m |} \right) \, ,
\\
&  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_k}} = -\varepsilon^3 \exp (-\varepsilon | s_j - s_m |) \frac {(s_{jr} - s_{mr})(s_{jk} - s_{mk})}{| s_j - s_m |} \, ,  \quad r \ne k \, ,
\\
& g''_{jj} = \varepsilon^2 \sum _{r=1}^n\  (e_{jr})^2  = \varepsilon^2 \,  , \quad j = 1, 2, \dots, n_2 \, ,
\end{aligned}
```
and for space ``H^{s_2}_\varepsilon (R^n)``:

```math
\tag{16}
\begin{aligned}
& h_i (\eta) =  \exp (-\varepsilon |\eta - p_i |) (3 + 3 \varepsilon |\eta - p_i | +  \varepsilon^2  |\eta - p_i |^2) )
 \, , \qquad \quad i = 1, 2, \dots, n_1 \, , \\
& h'_j (\eta) =\varepsilon^2 \exp (-\varepsilon |\eta - s_j | ) (1 + \varepsilon |\eta - s_j |) \sum _{k=1}^n (\eta_k - s_{jk}) e_{jk} \, , \quad  j = 1, 2, \dots, n_2 \, , \\
& g_{il}= \exp (-\varepsilon |p_i - p_l |) (3 + 3 \varepsilon |p_i - p_l | +  \varepsilon^2  |p_i - p_l |^2) ) \, ,
\\
& \qquad \qquad \qquad \qquad \qquad \qquad \qquad i = 1, 2, \dots, n_1 \, , \ \ l = 1, 2, \dots, n_1 \, , \\
& g'_{ij} = \varepsilon^2 \exp (-\varepsilon |p_i - s_j | ) (1 + \varepsilon |p_i - s_j |) \sum _{k=1}^n (p_{ik} - s_{jk}) e_{jk}  \, ,
\\
& \qquad \qquad \qquad \qquad \qquad \qquad \qquad  i = 1, 2, \dots, n_1 \, , \ \ j = 1, 2, \dots, n_2 \, ,
\\
& g''_{jm} =  \sum_{r=1}^n \sum_{k=1}^n  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_k}} e_{jk} e_{mr}
\\
& \quad \qquad j \ne m \, , \quad j = 1, 2, \dots, n_2 \, , \ \ m = 1, 2, \dots, n_2 \, ,
\\
& \text{where}
\\
&  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_r}} = \varepsilon^2 \exp (-\varepsilon | s_j - s_m |) (1 + \varepsilon | s_j - s_m | - \varepsilon^2 (s_{jr} - s_{mr})^2) \, ,
\\
&  \frac{ \partial^2 {V(s_j, s_m)} }{\partial{\eta_r} \partial{\xi_k}} = -\varepsilon^4 \exp (-\varepsilon | s_j - s_m |) (s_{jr} - s_{mr})(s_{jk} - s_{mk}) \, ,  \quad r \ne k \, ,
\\
& g''_{jj} = \varepsilon^2 \sum _{r=1}^n\  (e_{jr})^2  = \varepsilon^2 \,  , \quad j = 1, 2, \dots, n_2 \, ,
\end{aligned}
```

In a case when there is no information of function ``f`` derivatives the Problem ``(1)`` is reducing to the simplest interpolation problem:

*Problem:* Given points ``\{p_i, p_i \in R^n\}_{i=1}^{n_1}`` find a function ``f`` such that

```math
\tag{17}
\begin{aligned}
& f(p_i) =  u_i \, , \quad  i = 1, 2, \dots, n_1 \, ,
\\
& n_1 \gt 0 \, .
\end{aligned}
```
We assume that ``f`` is a continuous function. It can be treated as an element of Bessel Potential space ``H_\varepsilon^{s_0} (R^n) \, , s_0 =  n/2 + 1/2``, this space is continuously embedded in Hölder space ``C_b(R^n)`` of continuous and bounded functions.

Reproducing kernel of Bessel Potential space ``H_\varepsilon^{s_0}(R^n)``
can be written as:

```math
V(\eta , \xi) = \exp (-\varepsilon |\eta - \xi|) \, .
```
``(``with accuracy to a constant multipliler``)``, and expressions for $h_i, g_{il}$ are defined by:

```math
\tag{18}
\begin{aligned}
& h_i (\eta) = \exp (-\varepsilon |\eta - p_i |) \, , \quad  i = 1, 2, \dots, n_1 \, ,
 \\
& g_{il} = \exp (-\varepsilon | p_i - p_l | )) \ ,
 \ \quad
i = 1, 2, \dots, n_1 \, , \ \ l = 1, 2, \dots, n_1 \, .
\end{aligned}
```

When value of the parameter $\varepsilon$ is small this normal spline is similar to multivariate generalization of the one dimensional linear spline.

We now consider the choice of value for parameter $\varepsilon$. Approximating properties of the normal spline are getting better with smaller value of $\varepsilon$,
and if the value of this parameter is small enough the normal spline become similar to Duchon's $D^m -$spline [12]. However with decreasing value of $\varepsilon$ the condition number of the corresponding problem Gram matrix is increasing and the problem becomes numerically unstable. Therefore, when choosing the value of parameter $\varepsilon$, a compromise is needed. In practice, it is necessary to choose such value of the $\varepsilon$ that condition number of Gram matrix is small enough. Numerical procedures of the matrix condition number estimation are well known.

As well, it is useful to preprocess the source data of the problem by transforming the domain where interpolation nodes are located into the unit hypercube.

The normal splines method for one-dimensional function interpolation and linear ordinary differential and integral equations was proposed in [8] and [9] and developed in [10]. Multivariate generalization of the normal splines method was developed for two-dimensional problem of low-range computerized tomography in [15] and applied for solving a mathematical economics problem in [11]. Further results were reported on seminars and conferences [14,21,22].

**References**

[1]  R. Adams, J. Fournier, Sobolev Spaces. Pure and Applied Mathematics. (2nd ed.). Boston, MA: Academic Press, 2003.

[2] D. Adams, L. Hedberg, Function spaces and potential theory. Berlin, Springer, 1996.

[3] M. Agranovich, Sobolev Spaces, Their Generalizations and Elliptic Problems in Smooth and Lipschitz Domains, Springer, Switzerland, 2015.

[4] N. Aronszajn, Theory of reproducing kernels, Tranzactions of the AMS, Vol. 68, No. 3, 1950.

[5] N. Aronszajn, K.T. Smith, Theory of bessel potentials I, Ann.Inst.Fourier, Vol. 11, 1961.

[6] A. Balakrishnan, Applied Functional Analysis, New York, Springer-Verlag, 1976.

[7] A. Bezhaev, V. Vasilenko, Variational Theory of Splines, Springer US, 2001.

[8] V. Gorbunov, The method of normal spline collocation, USSR Computational Mathematics and Mathematical Physics, Vol. 29, No. 1, 1989.

[9] V. Gorbunov, Extremum Problems of Measurements Data Processing, Ilim, 1990 (in Russian).

[10] V. Gorbunov, V. Petrishchev, Improvement of the normal spline collocation method for linear differential equations, Comput. Math. Math. Phys., Vol. 43, No. 8, 2003.

[11] V. Gorbunov, I. Kohanovsky, K. Makedonsky, Normal splines in reconstruction of multi-dimensional dependencies. [Papers of WSEAS International Conference on Applied Mathematics, Numerical Analysis Symposium, Corfu, 2004](http://www.wseas.us/e-library/conferences/corfu2004/papers/488-312.pdf)

[12] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Springer, Berlin, Vol. 571, 1977.

[13] A. Ioffe, V. Tikhomirov, Theory of extremal problems, North-Holland, Amsterdam, 1979.

[14] I. Kohanovsky, Multidimensional Normal Splines and Problem of Physical Field Approximation, International Conference on Fourier Analysis and its Applications, Kuwait, 1998.

[15] I. Kohanovsky, Normal Splines in Computing Tomography, [Avtometriya, No.2, 1995](https://www.iae.nsk.su/images/stories/5_Autometria/5_Archives/1995/2/84-89.pdf)

[16] P.-J. Laurent, Approximation et optimization, Paris, 1972.

[17] H. Triebel, Interpolation. Function Spaces. Differential Operators. North-Holland, Amsterdam, 1978.

[18] R. Schaback, Kernel-based Meshless Methods, Lecture Notes, Goettingen, 2011.

[19] H. Wendland, Scattered Data Approximation. Cambridge University Press, 2005.

[20] [Reproducing Kernel of Bessel Potential space](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Reproducing-Kernel-of-Bessel-Potential-space).

[21] V. Gorbunov, I. Kohanovsky, Heterogeneous Parallel Method for the Construction of Multi-dimensional Smoothing Splines. [ESCO 2014 4th European Seminar on Computing, 2014](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)

[22] I. Kohanovsky, Inequality-Constrained Multivariate Normal Splines with Some Applications in Finance. [27th GAMM-Seminar Leipzig on Approximation of Multiparametric functions, 2011](https://www.ana.iusiani.ulpgc.es/proyecto2015-2017/pdfnew/ESCO2014_Book_of_Abstracts.pdf)

[23] G. Fasshauer, M. McCourt, Kernel-Based Approximation Methods Using Matlab, World Scientific, Singapore, 2015.

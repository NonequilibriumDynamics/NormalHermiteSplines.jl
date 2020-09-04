# Normal Splines Method

A problem of reconstruction of multivariate dependency under incomplete data set arises in many areas of research. Often there is a finite set of experimental measurements ``\{ \tilde u_i\}_{i=1}^m`` and we may treat an unknown function ``\varphi (x)`` as an element of a Hilbert space ``H(R^n)``. Consider the following data model

```math
\tag{1}
(f_i ,\varphi) =  u_i  \, , \qquad   \tilde  u_i\ = u_i + \epsilon_i \, \qquad        1 \le i \le  \ m,
```
where ``\{f_i\}`` are linear continuous functionals, ``\{u_i\}`` - "exact value of measurement", and ``\{\epsilon_i\}`` - random values uniformly distributed in ``\left[ -\delta , \delta \right]``. Hence we have a system of constraints

```math
\tag{2}
\left| (f_i ,\varphi) - \tilde u_i \right| \le \delta \, , \qquad 1 \le i \le \ m \, .
```
Problem of the function ``\varphi`` reconstruction is under-determined. We will approximate ``\varphi`` by constructing an element of minimal norm from the set of the system (2) solutions. Let's introduce a penalty functional ``J``

```math
(J, \varphi) = { \| \varphi\| }_H ^2 \ ,
```
and find a function ``\varphi`` in the Hilbert space ``H`` to minimize ``J`` subject to (2) (here ``\| \cdot \|_H`` denotes the ``H`` norm). Solution of this problem is a generalized spline of Atteia-Laurent [1]. We name it a normal spline following to the V. Gorbunov work [2] where the normal spline-collocation method for linear ordinary differential and integral equations was developed.

The normal spline method consists of a Hilbert norm minimization on the set of a collocation system solutions. In contrast to the classical collocation methods the basis system here is not given a priory, instead it is constructed in accordance with the chosen Hilbert space norm. The base functions are canonical images of the continuous linear functionals presented as inner product in the Hilbert space. Such functional presentation can be found if the Hilbert space ``H`` is a reproducing kernel Hilbert space and the corresponding reproducing kernel is known. In order to construct a reproducing kernel useful for a case of one-dimensional problems treated in [2] it was sufficient to suppose the Hilbert space ``H`` is a classical Sobolev space with integer order. Multivariate generalization of the normal spline method described in this blog was done ([6] — [10]) with usage of the Bessel potential spaces [3].

Essentially the normal spline method  is based on classical functional analysis results: the Sobolev and Bessel Potential space embedding theorems [4], the F.Riesz representation theorem [5] and Reproducing kernel properties.

Here and further we assume that ``\{f_i\}`` are linearly independent functionals therefore the set of the system (2) solutions is not empty, convex and closed one. It is known that every closed convex set in a Hilbert space has a unique element of minimal norm [1, 5] - here it is a uniform smoothing normal spline:
```math
\tag{3}
      \sigma = {\mathop{\rm arg\,min}\nolimits} \lbrace {\| \varphi \| }^2_H : (2) \rbrace \ .
```
In general a normal spline can be treated as a projection of an element of a Hilbert space to a closed convex set in that Hilbert space. Thereby the problem (3) always has the unique solution.

In accordance with Riesz representation theorem [5] every linear continuous functional ``f_i`` on a Hilbert space ``H`` can be represented as inner product of some element ``h_i \in H`` and ``\varphi \in H``, for any ``\varphi \in H`` :
```math
      (f_i ,\varphi) = {\langle h_i , \varphi \rangle}_H \, , \qquad  \forall \varphi \in H \  ,
```
where ``{\langle \cdot , \cdot \rangle}_H`` - inner product in ``H``. Then (2) can be written in form:
```math
\tag{4}
    | {\langle h_i , \varphi \rangle}_H - \tilde u_i | \le \delta \, ,     \qquad 1 \le i \le m \  .
```
and problem (2) is reduced to:
```math
\tag{5}
      \sigma = {\mathop{\rm arg\,min}\nolimits} \lbrace {\langle \varphi , \varphi \rangle}_H : (4) \rbrace \ .
```
In accordance with extremum conditions for the problem (5) its solution can be presented in form [1]:
```math
    \sigma =  \sum _{j=1} ^m (\mu _j - \mu _{j+m}) h_j \ , \quad  \quad \mu _j  \le 0 , \quad \mu _{j+m} \le 0 , \quad 1\le j \le m .
```
It allows to reduce the initial problem (5) of constructing a uniform smoothing normal spline to solving a finite-dimensional quadratic programming problem:
```math
\tag{6}
\begin{aligned}
   &   \sigma = {\mathop{\rm arg\,min}} \Big\lbrace {\sum _{i=1} ^m \sum _{j=1} ^m   (\mu _i - \mu _{i+m}) (\mu _j - \mu _{j+m}) g_{ij}} : (10), (11) \Big\rbrace \ ,
\\  &  \Big| \sum _{j=1} ^m g_{ij} (\mu _j - \mu _{j+m}) - \tilde u_i \Big| \, \le \, \delta \, ,    \qquad 1 \le i \le m \  ,
\\  & \mu _j  \le 0 , \quad \mu _{j+m} \le 0 , \quad 1\le j \le m \ ,
\end{aligned}
```
here ``g_{ij}={\langle h_i , h_j \rangle}_H `` - coefficients of the symmetric Gram matix of the set of the elements ``\{h_i\}, h_i \in H``. Notice that elements ``\{h_i\}`` (images of functionals ``\{f_i\}``) are linearly independent therefore the Gram matrix is a positive definite one.

It was shown [2] that it is not necessary to formulate the problem (6) in its explicit form. A special version of algorithm for solving this simple quadratic programming problem will be described in the next posts.

Results related to multidimensional normal spline method and its applications were published in works [6, 10] and presented at conferences [8, 9, 11, 12].

**References**

[1] P.-J. Laurent, Approximation et optimization, Paris, 1972.

[2] V. Gorbunov, A Method of  Normal Spline Collocation, Computational Mathematics and Mathematical Physics. 1989. Vol.29. No.2.

[3] N. Aronszajn, K. Smith, Theory of bessel potentials I, Ann.Inst.Fourier, 11, 1961. 

[4] S. Sobolev, Some applications of functional analysis in mathematical physics, AMS, 2008.

[5] A. Balakrishnan, Applied Functional Analysis, New York, Springer-Verlag, 1976.

[6] I. Kohanovsky, Normal Splines in Computing Tomography. Avtometriya, No. 2, 1995 (https://www.iae.nsk.su/images/stories/5_Autometria/5_Archives/1995/2/84-89.pdf).

[7] I. Kohanovsky, Data approximation using multidimensional normal splines, Unpublished manuscript, 1996.

[8] I. Kohanovsky, Multidimensional Normal Splines and Problem of Physical Field Approximation, International Conference on Fourier Analysis and its Applications, Kuwait, 1998.

[9] I. Kohanovsky, Normal splines in fractional order Sobolev spaces and some of its
applications, The Third Siberian Congress on Applied and Industrial mathematics 
(INPRIM-98), Novosibirsk, 1998.

[10] V. Gorbunov, I. Kohanovsky, K. Makedonsky, Normal splines in reconstruction of multi-dimensional dependencies. Papers of WSEAS International Conference on Applied Mathematics, Numerical Analysis Symposium, Corfu, 2004. (http://www.wseas.us/e-library/conferences/corfu2004/papers/488-312.pdf)

[11] I. Kohanovsky, Inequality-Constrained Multivariate Normal Splines with Some Applications in Finance. 27th GAMM-Seminar Leipzig on Approximation of Multiparametric functions, Max-Planck-Institute for Mathematics in the Sciences, Leipzig, Germany, 2011.

[12] V. Gorbunov, I. Kohanovsky, Heterogeneous Parallel Method for the Construction of Multi-dimensional Smoothing Splines. ESCO 2014 4th European Seminar on Computing. University of West Bohemia, Plzen, Czech Republic, 2014.


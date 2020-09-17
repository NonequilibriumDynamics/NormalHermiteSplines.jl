# Relation to Polyharmonic Splines

Interpolating normal spline ``\sigma`` is the solution of the variational problem:

```math
\tag{1}
  \| f \|^2_ {H^s_\varepsilon} = \int ( \varepsilon ^2  + | \xi |^2 )^s | {\mathcal F} [f(\xi)] |^2  \, d \xi \ \to min  \, , \qquad \forall f \in H^s_\varepsilon (R^n) \ , \quad s > \frac{n}{2} \ , 
```
```math
\tag{2}
 f(p_i) =  u_i \, , \quad  p_i \in R^n \, , \qquad i = 1, 2, \dots, n  \qquad \qquad\qquad\qquad\qquad\qquad\qquad 
```
here ``H^s_\varepsilon (R^n)`` is the Bessel Potential space, which is defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \quad  s > \frac{n}{2}  \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. Space
 ``H^s_\varepsilon`` is a Hilbert space (see [Interpolating Normal Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Splines/)).


It is known that for any positive integer ``m``, the space ``H^m_\varepsilon`` consists of all
square integrable functions whose derivatives in the sense of distributions up to
order ``m`` are square integrable and therefore in such case this space coincides with Sobolev space ([1], [3]).

Polyharmonic ``D^m`` spline ``\sigma_{D^m}`` is the the result of minimization of the quadratic functional (Sobolev semi-norm) ([2], [4])
```math
\tag{3}
 \int \sum_{|\alpha| = m} \frac{m!}{\alpha!} |D^\alpha f(x) |^2   \, d x \  \to min  \, , \qquad \forall f \in W^m_2 (R^n) \ , \quad m > \frac{n}{2}
```
under interpolation constraints
```math
\tag{4}
 f(p_i) =  u_i \, , \quad  p_i \in R^n \, , \qquad i = 1, 2, \dots, n  \qquad \qquad\qquad\qquad\qquad 
```
here ``W^m_2 (R^n)`` is the Sobolev space.

As it is pointed out in [4] we may replace the minimizing functional (3) with 

```math
 \int \sum_{|\alpha| = m} \frac{m!}{\alpha!} |{\mathcal F} [D^\alpha f] |^2   \, d \xi \  \to min  \, , \qquad \forall f \in W^m_2 (R^n) \ , \quad m > \frac{n}{2} 
```
``(`` in view of Parseval’s idenitity``)``. Further, as ``{\mathcal F} [D^\alpha f] = (iξ)^\alpha{\mathcal F} [f] `` we can write (5) as

```math
 \int |{\mathcal F} [f] |^2 \sum_{|\alpha| = m} \frac{m!}{\alpha!}  \xi^{2 \alpha}   \, d \xi \  \to min  \, , \qquad \forall f \in W^m_2 (R^n) \ , \quad m > \frac{n}{2} 
```
taking into account that ``\sum_{|\alpha| = m} \frac{m!}{\alpha!}  \xi^{2 \alpha} = |\xi|^{2m}`` we get that polyharmonic ``D^m`` spline ``\sigma_{D^m}`` can be found as the solution of the variational problem

```math
\tag{5}
  \int | \xi |^{2m} | {\mathcal F} [f(\xi)] |^2  \, d \xi \ \to min  \, , \qquad \forall f \in W^m_2 (R^n) \ , \quad m > \frac{n}{2}

```
```math
\tag{6}
 f(p_i) =  u_i \, , \quad  p_i \in R^n \, , \qquad i = 1, 2, \dots, n  \qquad \qquad\qquad\qquad 
```






**References**

[1] M. Agranovich, Sobolev Spaces, Their Generalizations and Elliptic Problems in Smooth and Lipschitz Domains, Springer, Switzerland, 2015.

[2] A. Bezhaev, V. Vasilenko, Variational Theory of Splines, Springer US, 2001.

[3] M. Buhmann. Radial Basis Functions. Cambridge University Press, 2003.

[4] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Vol. 571, Springer, Berlin, 1977


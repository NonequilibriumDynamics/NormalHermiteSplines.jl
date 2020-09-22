# Comparison with Polyharmonic Splines

Interpolating normal spline ``\sigma`` is solution of the variational problem:

```math
\tag{1}
  \| f \|^2 = \int ( \varepsilon ^2  + | \xi |^2 )^s | {\mathcal F} [f(\xi)] |^2  \, d \xi \ \to min  \, , \qquad \forall f \in H^s_\varepsilon (R^n) \ , \quad s > \frac{n}{2} \ , 
```
```math
\tag{2}
 f(p_i) =  u_i \, , \quad  p_i \in R^n \, , \qquad i = 1, 2, \dots, n  \qquad \qquad\qquad\qquad\qquad\qquad\qquad 
```
here ``H^s_\varepsilon (R^n)`` is Bessel Potential space, which is defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \quad  s > \frac{n}{2}  \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. Space
 ``H^s_\varepsilon`` is a Hilbert space (see [Interpolating Normal Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Splines/)).


It can be shown that for any positive integer ``m`` the space ``H^m_\varepsilon (R^n)`` consists of all square integrable functions whose derivatives in the sense of distributions up to order m are square integrable [1]. The norm on ``H^m_\varepsilon (R^n)`` can be defined by
```math
\| \varphi \|' = \left( \int \Big [ \varepsilon^2 | \varphi (x) |^2  + \sum_{|\alpha| = m} \frac{m!}{\alpha!} |D^\alpha \varphi (x) |^2 \Big ]
  \, d x \ 
 \right)^{1/2}  .
```
The corresponding inner product has the form
```math
\langle \varphi , \psi \rangle' =
\left( \int \Big [ \varepsilon^2 \varphi (x) \overline{\psi (x)}  + \sum_{|\alpha| = m} \frac{m!}{\alpha!} D^\alpha \varphi (x) \overline{D^\alpha \psi (x)}  \Big ]
  \, d x \ 
 \right)^{1/2} ,
```
here ``\alpha = (\alpha_1, \dots, \alpha_n )`` is multi-index with nonnegative integral
entries, ``| \alpha | = \alpha_1 + \dots + \alpha_n``, ``\, \alpha ! = \alpha_1 ! \dots \alpha_n !`` and ``D^\alpha \varphi (x) = \frac{ \partial^{|\alpha|}{\varphi} }{ \partial{x_1^{\alpha_1}} \dots x_n^{\alpha_n}}``.

The norms ``\| \varphi \|`` and ``\| \varphi \|'`` are equivalent and space ``H^m_\varepsilon (R^n)`` coincides with Sobolev space. Therefore the problem (1), (2) can be written as

```math
\tag{3}
  \| f \|'^2 = \int \Big [ \varepsilon^2 | f(x) |^2  + \sum_{|\alpha| = m} \frac{m!}{\alpha!} |D^\alpha f(x) |^2 \Big ]
  \, d x  \ \to min  \, , \ \  \forall f \in W^m_2 (R^n) \ , \ \varepsilon \gt 0 \, , \ m \gt \frac{n}{2} \ , 
```
```math
\tag{4}
 f(p_i) =  u_i \, , \quad  p_i \in R^n \, , \qquad i = 1, 2, \dots, n  \qquad \qquad \qquad \qquad\qquad \qquad\qquad\qquad  
```
here ``W^m_2 (R^n)`` is Sobolev space. 

Normal spline ``\sigma`` always exist if all points ``\{p_i\}`` are different and it is unique.

Polyharmonic ``D^m`` spline ``\sigma_{D^m}`` is the the result of minimization of the quadratic functional (Sobolev semi-norm) ([3], [4])
```math
\tag{5}
 \int\limits_\Omega \sum_{|\alpha| = m} \frac{m!}{\alpha!} |D^\alpha f(x) |^2   \, d x \  \to min  \, , \qquad \forall f \in W^m_2 (\Omega) \ , \quad m > \frac{n}{2}
```
under interpolation constraints
```math
\tag{6}
 f(p_i) =  u_i \, , \quad  p_i \in \Omega \, , \qquad i = 1, 2, \dots, n  \qquad \qquad\qquad\qquad\qquad 
```
here ``W^m_2 (\Omega)`` is Sobolev space (``\Omega \subset R^n`` is a bounded and
sufficiently regular ([3, 4]) domain in ``R^n``).

``D^m`` spline ``\sigma_{D^m}`` always exist if all points ``\{p_i\}`` are different and it is unique if set ``\{p_i\}`` contains ``P_{m-1}``- unisolvent set [3].

Comparing the variational problems (3),(4) and (5),(6) we can expect that their solutions – the normal spline ``\sigma`` and ``D^m`` spline ``\sigma_{D^m}`` will have the similar properties when ``\varepsilon \to 0`` in (3).

Also, in general case (when ``s`` is not an integer number), we could expect that once the normal spline scaling ("shape") parameter ``\varepsilon`` is small enough the normal spline ``\sigma`` will have similar properties as ``D^{m,s}`` spline ``\sigma_{D^{m,s}}`` of corresponding smoothness constructed in the appropriate intermediate Beppo-Levi space ([2, 4]).

**References**

[1] M. Agranovich, Sobolev Spaces, Their Generalizations and Elliptic Problems in Smooth and Lipschitz Domains, Springer, Switzerland, 2015.

[2] F. Utreras, Recent results on multivariate smoothing splines, in Multivariate
Approximation and Interpolation, ISNM 94, W. Haussmann and K. Jetter (eds.), Birkhauser, Basel, 299–312, 1990.

[3] A. Bezhaev, V. Vasilenko, Variational Theory of Splines, Springer US, 2001.

[4] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Vol. 571, Springer, Berlin, 1977


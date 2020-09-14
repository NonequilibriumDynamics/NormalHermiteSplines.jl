# Relation to Polyharmonic Splines

Interpolating normal spline ``\sigma`` is  the solution of the variational problem:

```math
\tag{1}
\begin{aligned}
 & \| f \|^2_ {H^s_\varepsilon} = \int ( \varepsilon ^2  + | \xi |^2 )^s \mathcal | F [f ] |^2  \, d \xi \ \to min  \, , \qquad \forall f \in H^s_\varepsilon (R^n) \ ,
\\ 
 &  f(p_i) =  u_i \, , \quad  i = 1, 2, \dots, n \, ,
\\
\end{aligned}
```
here ``H^s_\varepsilon`` is the Bessel Potential space, which is defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \ \ s = n/2 + 1/2 + r \, , \quad r = 1,2,\dots \, .
```
where ``| \cdot |`` is the Euclidean norm, ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. Space
 ``H^s_\varepsilon`` is a Hilbert space (see [Interpolating Normal Splines](https://igorkohan.github.io/NormalHermiteSplines.jl/stable/Interpolating-Normal-Splines/)).







**References**

[1] M. Agranovich, Sobolev Spaces, Their Generalizations and Elliptic Problems in Smooth and Lipschitz Domains, Springer, Switzerland, 2015.

[2] M. Buhmann. Radial Basis Functions. Cambridge University Press, 2003.

[3] J. Duchon, Interpolation des fonctions de deux variables suivante le principede la flexion des plaques minces, Rev. Franc¸aise Automat. Informat. Rech. Oper. Anal. Numer.,  Vol. 10, No.3, 1976.

[4] J. Duchon, Splines minimizing rotation-invariant semi-norms in Sobolev spaces, Lect. Notes in Math., Vol. 571, Springer, Berlin, 1977


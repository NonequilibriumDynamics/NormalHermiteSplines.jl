
# Reproducing Kernel of Bessel Potential space

The standard definition of Bessel Potential space ``H^s`` can be found in ([1], [2], [6], [11], [12]). Here the normal splines will be constructed in the Bessel Potential space ``H^s_\varepsilon`` defined as:

```math
   H^s_\varepsilon (R^n) = \left\{ \varphi | \varphi \in S' ,
  ( \varepsilon ^2 + | \xi |^2 )^{s/2}{\mathcal F} [\varphi ] \in L_2 (R^n) \right\} , \quad
  \varepsilon \gt 0 , \  s \gt \frac{n}{2} .
```
where ``S'  (R^n)`` is space of L. Schwartz tempered distributions, parameter ``s`` may be treated as a fractional differentiation order and ``\mathcal F [\varphi ]`` is a Fourier transform of the ``\varphi``. The parameter ``\varepsilon`` introduced here may be considered as a "scaling parameter". It allows to control approximation properties of the normal spline which usually are getting better with smaller values of ``\varepsilon``, also it may be used to reduce the ill-conditioness of the related computational problem (in traditional theory ``\varepsilon = 1``).

Theoretical properties of spaces ``H^s_\varepsilon`` at ``\varepsilon \gt 0`` are identical — they are Hilbert spaces with inner product

```math
\langle \varphi , \psi \rangle _{H^s_\varepsilon} =
\int ( \varepsilon ^2  + | \xi |^2 )^s
\mathcal F [\varphi ] \overline{\mathcal F [\psi ] } \, d \xi
```
and  norm

```math
\| \varphi \|_ {H^s_\varepsilon} = \left( \int ( \varepsilon ^2  + | \xi |^2 )^s
\mathcal | F [\varphi ] |^2  \, d \xi \ 
 \right)^{1/2} \ .
```
It is easy to see that all ``\| \varphi \|_{H^s_\varepsilon}`` norms are equivalent. It means that space ``H^s_\varepsilon (R^n)`` is equivalent to ``H^s (R^n) =  H^s_1 (R^n)``.

Let's describe the Hölder spaces ``C_b^t(R^n), t \gt 0`` ([9], [2]).

*Definition 1.* We denote the space

```math
   S(R^n) = \left\{ f | f \in C^\infty (R^n) ,  \sup_{x \in R^n} | x^\alpha D^\beta f(x) |  \lt \infty , \forall \alpha, \beta \in \mathbb{N}^n  \right\}
```
as Schwartz space (or space of complex-valued rapidly decreasing infinitely differentiable functions defined on ``R^n``) ([6], [7]).

Below is a definition of Hölder space ``C^t_b(R^n)`` [9]:

*Definition 2.* If ``0 \lt t = [t] + \{t\}, [t]`` is non-negative integer, ``0 \lt \{t\} \lt 1``, then ``C^t_b(R^n)`` denotes the completion of ``S(R^n)`` in the norm

```math
\begin{aligned}
   C^t_b (R^n) &= \left\{ f | f \in C^{[t]}_b (R^n) ,  \| f  \|_{C^t_b}  \lt \infty \right\} , \\

  \| f  \|_{C^t_b} &= \| f  \|_{C_b^{[t]}} \, + \,

\sum _{|\alpha | = [t]} \sup _{x \ne y} \frac {| D^\alpha  f(x)  - D^\alpha  f(y) | } { | x - y |^{\{t\}}} \ , \\

  \| f  \|_{C_b^{[t]}} &= \sup _{x \in R^n} | D^\alpha f(x) |, \, \forall \alpha : | \alpha | \le [t].
\end{aligned}
```
Space ``C^{[t]}_b (R^n)`` consists of all functions having bounded continuous derivatives up to order ``[t]``. It is easy to see that ``C_b^t(R^n)`` is Banach space [9].

Connection of Bessel Potential spaces ``H^s(R^n)`` with the spaces ``C_b^t(R^n)`` is expressed in Embedding theorem ([9], [2]).

*Embedding Theorem:* If ``s = n/2+t``, where ``t`` non-integer, ``t \gt 0``, then space ``H^s(R^n)`` is continuously embedded in ``C_b^t(R^n)``.

Particularly from this theorem follows that if ``f \in H^{n/2 + 1/2}_\varepsilon (R^n)``, corrected if necessary on a set of Lebesgue measure zero, then it is uniformly continuous and bounded. Further if ``f \in H^{n/2 + 1/2 + r}_\varepsilon (R^n)``, ``r`` — integer non-negative number, then it can be treated as ``f \in C^r (R^n)``, where ``C^r (R^n)`` is a class of functions with ``r`` continuous derivatives.

It can be shown ([3], [11], [8], [4], [5]) that function 

```math
\begin{aligned}
 & V_s ( \eta , x, \varepsilon ) = c_V (n,s,\varepsilon) (\varepsilon |\eta - x | )^{s - \frac{n}{2}}

          K_{s - \frac{n}{2}} (\varepsilon |\eta - x | )     \ ,
\\
 & c_V (n,s,\varepsilon) = \frac{\varepsilon ^{n-2s}} { 2^{s-1} (2 \pi )^{n/2} \Gamma (s) }, \ \eta \in R^n, \  x \in R^n, \ \varepsilon \gt 0 ,  s \gt \frac{n}{2}
\end{aligned}
```
is a reproducing kernel of ``H^s_\varepsilon (R^n)`` space. Here ``K_{\gamma}`` is modified Bessel function of the second kind [10]. The exact value of ``c_V (n,s,\varepsilon)`` is not important here and will be set to ``\sqrt{\frac{2}{\pi}}`` for ease of further calculations. 

This reproducing kernel is known as Matérn kernel [4,13].

The kernel ``K_{\gamma}`` becomes especially simple when ``\gamma``  is half-integer. 

```math

 \gamma =  r  + \frac{1}{2} \ , (r = 0, 1, \dots ).
```
In this case it is expressed via elementary functions (see [10]):

```math
\begin{aligned}
K_{r+1/2}(t) &=
\sqrt{\frac{\pi} {2t}} t^{r+1} \left (
- \frac{1}{t} \frac{d}{dt} \right )^{r+1} \exp (-t) \ ,
\\
K_{r+1/2}(t) &=
\sqrt{\frac{\pi} {2t}} \exp (-t) \sum_{k=0}^r
\frac{(r+k)!}{k! (r-k)! (2t)^k} \ ,
\  (r = 0, 1, \dots ) \  .
\end{aligned}
```
Let ``s_r =  r + \frac{n}{2} + \frac{1}{2}, \  r = 0, 1, \dots``, then ``H^{s_r}_\varepsilon(R^n)`` is continuously embedded in ``C_b^r(R^n)`` and its reproducing kernel with accuracy to constant multiplier can be presented as follows

```math
\begin{aligned}
V_{r + \frac{n}{2} + \frac{1}{2}}(\eta , x, \varepsilon) &= \exp (-\varepsilon |\eta - x |) 
\sum_{k=0}^{r} \frac{(r+k)!}{2^k k! (r-k)!} (\varepsilon |\eta - x |)^{r-k} \ ,
\\
&  (r = 0, 1, \dots ) \  .
\end{aligned}
```

In particular we have:

```math
\begin{aligned}
& V_{\frac{n}{2} + \frac{1}{2}}(\eta , x, \varepsilon) = \exp (-\varepsilon |\eta - x |) \ ,
\\
& V_{1 + \frac{n}{2} + \frac{1}{2}}(\eta , x, \varepsilon) = \exp (-\varepsilon |\eta - x |)
(1 + \varepsilon |\eta - x |) \ ,
\\
& V_{2 + \frac{n}{2} + \frac{1}{2}}(\eta , x, \varepsilon) = \exp (-\varepsilon |\eta - x |)
(3 + 3\varepsilon |\eta - x | + \varepsilon ^2 |\eta - x | ^2 ) \ .
\end{aligned}
```

**References**

[1] D. Adams, L. Hedberg, Function spaces and potential theory. Berlin, Springer, 1996.

[2] M. Agranovich, Sobolev Spaces, Their Generalizations and Elliptic Problems in Smooth and Lipschitz Domains, Springer, Switzerland, 2015.

[3] N. Aronszajn, K. Smith, Theory of bessel potentials I, Ann.Inst.Fourier, 11, 1961.

[4] G. Fasshauer, Green’s Functions: Taking Another Look at Kernel Approximation, Radial Basis Functions, and Splines. In: Neamtu M., Schumaker L. (eds) Approximation Theory XIII: San Antonio 2010. Springer Proceedings in Mathematics, vol 13. Springer, New York, 2012.

[5] I. Kohanovsky, Multidimensional Normal Splines and Problem of Physical Field Approximation, International Conference on Fourier Analysis and its Applications, Kuwait, 1998.

[6] S. Nikol'skiĭ, Approximation of functions of several variables and imbedding theorems, Grundl. Math. Wissensch., 205, Springer-Verlag, New York, 1975.

[7] M. Reed, B. Simon, Methods of Modern Mathematical Physics, I: Functional Analysis, Academic Press, 1972.

[8] R. Schaback, Kernel-based Meshless Methods, Lecture Notes, Goettingen, 2011.

[9] H. Triebel, Interpolation. Function Spaces. Differential Operators, North-Holland, Amsterdam, 1978.

[10] G. Watson, A Treatise on the Theory of Bessel Functions ( 2nd.ed.), Cambridge University Press, 1966.

[11] H. Wendland, Scattered Data Approximation. Cambridge University Press, 2005.

[12] J. Lions, E. Magenes, Problemes Aux Limites Non-Homogenes et Applications Vol. 1, Dunod, Paris, 1968.

[13] G. Fasshauer, M. McCourt, Kernel-Based Approximation Methods Using Matlab, World Scientific, Singapore, 2015.


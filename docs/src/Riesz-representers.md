# The Riesz representation of functionals and a reproducing kernel Hilbert space

In this post we discuss a way of constructing the Riesz representers of the continuous linear functionals $f_i$ assuming the Hilbert space $H$ is a reproducing kernel Hilbert space. Let's recall the Riesz representation theorem and the reproducing kernel Hilbert space definition.

Riesz representation theorem ([1]): If $f$ is a linear continuous functional on a Hilbert space $H$ then there exists some $h \in H$ such that for every $\varphi \in H$ we have
```math
     (f, \varphi) = {\langle \varphi , h \rangle}_H
```
Reproducing Kernel Hilbert space ([2]): A Hilbert space $H(R^n)$ is called a reproducing kernel Hilbert space (RKHS) if there is a reproducing kernel function $K(\eta, x)$ of $\eta$ and $x$ in $R^n$ such that:
1) For any $x \in R^n$ the function $K(\eta, x)$ belongs to $H(R^n)$ as a function of the $\eta$.
2) The reproducing property. For any $x \in R^n$ and any $\varphi \in H(R^n)$, the following equality is valid:
$\qquad \qquad  \qquad  \varphi (x) = {\langle \varphi(\eta) , K(\eta , x) \rangle}_H$

Inner product here applies to functions of $\eta$. It is known, if a reproducing kernel exists it is unique and it is symmetric with respect to the $\eta$ and $x$ ([2]): $K(\eta , x) = K(x, \eta)$.

Let $K(\eta, x)$ is a reproducing kernel of the Hilbert space $H(R^n)$, then we can find the Riesz representers (images) $h_i$ of the functionals $f_i$. Namely:
```math
      h_i (x) = (f_i, K (\cdot, x)) \, .
```
Indeed, by reproducing property:
```math
      h_i (x) = {\langle h_i , K(\cdot , x) \rangle}_H
```
but since the $h_i$ is a Riesz representer of the $f_i$ (see (5) in the previous post):
```math
      {\langle h_i , K(\cdot , x) \rangle}_H = (f_i , K(\cdot , x) ) \,
```
(here $K(\eta, x) \in H$ as a function of the $\eta$).
       Then the coefficients $g_{ij}$ of the Gram matrix of the set of elements $\{h_i\}$ are defined as follows:
```math
  g_{ij} = {\langle h_i , h_j \rangle}_H = (f_i , h_j) =  \bigl(f_i , ( f_j , K ) \bigr) \ .
```
In the next post a reproducing kernel of the Bessel Potential space will be presented.

**References**

[1] A. Balakrishnan. Applied Functional Analysis. // New York: Springer-Verlag, 1976.

[2] N. Aronszajn, Theory of reproducing kernels, Tranzactions of the AMS.– 950 – Vol.68.



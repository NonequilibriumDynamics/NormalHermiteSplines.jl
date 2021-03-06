# The Riesz representation of functionals and a reproducing kernel Hilbert space

In this post we discuss a way of constructing a Riesz representer of the continuous linear functional ``f`` on a Hilbert space ``H``, assuming the space ``H`` is a reproducing kernel Hilbert space. Let's recall the Riesz representation theorem and the reproducing kernel Hilbert space definition.

Riesz representation theorem ([1]): If ``f`` is a linear continuous functional on a Hilbert space ``H`` then there exists some ``h \in H`` such that for every ``\varphi \in H`` we have
```math
     (f, \varphi) = {\langle \varphi , h \rangle}_H
```
Reproducing Kernel Hilbert space definition ([2]): A Hilbert space ``H(R^n)`` is called a reproducing kernel Hilbert space (RKHS) if there is a reproducing kernel function ``K(\eta, x)`` of ``\eta`` and ``x`` in ``R^n`` such that:
1) For any ``x \in R^n`` the function ``K(\eta, x)`` belongs to ``H(R^n)`` as a function of the ``\eta``.
2) The reproducing property: for any ``x \in R^n`` and any ``\varphi \in H(R^n)``, the following equality is valid:
``\qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad  \qquad  \varphi (x) = {\langle \varphi(\eta) , K(\eta , x) \rangle}_H``

Inner product here applies to functions of ``\eta``. It is known, if a reproducing kernel exists it is unique and it is symmetric with respect to the ``\eta`` and ``x`` ([2]): ``K(\eta , x) = K(x, \eta)``.

Let ``K(\eta, x)`` is a reproducing kernel of the Hilbert space ``H(R^n)``, then we can find the Riesz representer ``h`` of the functional ``f``. Namely:
```math
      h (x) = (f, K (\cdot, x)) \, .
```
Indeed, by reproducing property:
```math
      h (x) = {\langle h , K(\cdot , x) \rangle}_H
```
but since the ``h`` is a Riesz representer of the ``f``:
```math
      {\langle h , K(\cdot , x) \rangle}_H = (f , K(\cdot , x) ) \,
```
here ``K(\eta, x) \in H`` as a function of the ``\eta``.

 Let's we have a set of the continuous linear functionals ``\{f_i\}`` on Hilbert space ``H`` and the corresponding set of their Riesz representers ``\{h_i\}``. Then coefficients ``g_{ij}`` of the Gram matrix of elements ``\{h_i\}`` can be found as follows:
```math
 \qquad \qquad  g_{ij} = {\langle h_i , h_j \rangle}_H = (f_i , h_j) =  \bigl(f_i , ( f_j , K ) \bigr) \ .
```

**References**

[1] A. Balakrishnan. Applied Functional Analysis, New York, Springer-Verlag, 1976.

[2] N. Aronszajn, Theory of reproducing kernels, Transactions of the American Mathematical Society, Vol. 68, No. 3, 1950.


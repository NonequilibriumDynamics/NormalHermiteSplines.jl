# Numerical Tests

## 2D interpolation case

### Models

To test the accuracy of the normal spline method, we take three test functions defined in the region ``{Ω=[0,1]\times[0,1]}``:

1. Franke Exponential function ``f_1`` [1]

```math
\begin{aligned}
f_1(x,y) =& \ 0.75*\exp(-((9*x - 2)^2 + (9*y - 2)^2)/4) + 0.75*\exp(-((9*x + 1)^2)/49 - (9*y + 1)/10) \, +
\\
& \ 0.5*\exp(-((9*x - 7)^2 + (9*y - 3)^2)/4) - 0.2*\exp(-(9*x - 4)^2 - (9*y - 7)^2)
\end{aligned}
```

2. Franke Cliff function ``f_2`` [1]

```math
    f_2 (x,y) = \frac{\tanh(9*(y - x)) + 1}{\tanh(9)+ 1} \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad \qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad
```

3. Franke nad Nielson "faults and creases" model ``f_3`` [2] 

```math
f_3 (x,y) = \begin{cases}
            0.5 \, ,  \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad \ \ y \le 0.4  \cr
            0.5(1 - ((y - 0.4)/0.6)^2)\, ,  \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad \  x \le 0.1 \, , \ y \gt 0.4  \cr
            0.5((y - 1)/0.6)^2 (1 - x) / 0.8 \, ,  \qquad\qquad\qquad\qquad\qquad\qquad\qquad \ \ \,   x \ge 0.2 \, , \ y \gt 0.4  \cr
            0.5((y - 1)/0.6)^2 (x - 0.1) + (1 - ((y - 0.4)/0.6)^2) (0.2 - x)) \, ,  \quad  0.1 \lt x \le 0.2 \, , \ y \gt 0.4  \cr
            \end{cases} 
            \quad
```


**References**

[1] R. Franke, A critical comparison of some methods for interpolation of scattered
data, NPS-53-79-003, Dept. of Mathematics, Naval Postgraduate School, Monterey, CA, 1979.

[2] R. Franke and G. Nielson, Surface approximation with imposed conditions, Computer Aided Geometric Design, North Holland Pubi. Co., 1983. 

### Grids

The test functions are sampled on two kinds of the data sets: Halton [1] and non-uniform random points.


**References**

[1] J.Halton, On the efficiency of certain quasi-random sequences of points in evaluating multi-dimensional integrals, Numer. Math., Vol.2, 1960.



### Tests



## 3D interpolation case

### Models

### Grids

### Tests


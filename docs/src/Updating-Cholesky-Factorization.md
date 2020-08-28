# Algorithms for updating Cholesky factorization

A problem of updating Cholesky factorization was treated in [1], [2] and [3]. Here algorithms taken from these works are described. These algorithms are for computing a factorization ``\tilde A = \tilde L \tilde L^T`` where ``\tilde A`` is a matrix received of the matrix ``A = L L^T`` after a row and the symmetric column were added or deleted from ``A``.
Let ``A \in R^{n \times n}`` is a symmetric positive definite matrix. Applying Cholesky method to ``A`` yields the factorization ``A = L L^T`` where ``L`` is a lower triangular matrix with positive diagonal elements.
Suppose ``\tilde A`` is a positive definite matrix created by adding a row and the symmetric 
column to ``A``:

```math
    \tilde A = \left( \begin{array}{cc}
       A & d \\
       d^T & \gamma
      \end{array}
      \right) \ , \qquad d^T  \in R^n
```
Then its Cholesky factorization is [2]:

```math
    \tilde L = \left( \begin{array}{cc}
       L  \\
       e^T & \alpha
      \end{array}
      \right) \ ,
```
where

```math
     e = L^{-1} d \ ,
```
and

```math
    \alpha = \sqrt{ \tau } \ , \quad \tau = \gamma - e^T e, \quad \tau > 0.
```
If ``\tau \le 0`` then ``\tilde A`` is not positive definite.

Now assume ``\tilde A`` is obtained by deleting the ``r^{th}`` row and column from ``A``. Matrix ``\tilde A`` is the positive definite matrix (as any principal square submatrix of the positive definite matrix). It is shown in [1],[3] how to get ``\tilde L`` from ``L`` in such case. Let us partition ``A`` and ``L`` along the ``r^{th}`` row and column:

```math
    A = \left( \begin{array}{ccc}
          A_{11} & a_{1r} & A_{12} \\
          a_{1r}^T & a_{rr} & a_{2r}^T \\
          A_{21} & a_{2r} & A_{22}
      \end{array}
     \right) \ , \qquad
    L = \left( \begin{array}{ccc}
   L_{11}   \\
   l_{1r}^T & l_{rr}   \\
   L_{21} & l_{2r} & L_{22}
    \end{array}
     \right) \ .
```
Then ``\tilde A`` can be written as

```math
    \tilde A = \left( \begin{array}{cc}
    A_{11}  & A_{12} \\
    A_{21}  & A_{22}
    \end{array}
     \right) \ .
```
By deleting the ``r^{th}`` row from ``L`` we get the matrix

```math
    H = \left( \begin{array}{ccc}
  L_{11}     \\
  L_{21} & l_{2r} & L_{22}
    \end{array}
     \right) 
```
with the property that ``H H^T = \tilde A``. Factor ``\tilde L`` can be received from ``H`` by applying Givens rotations to the matrix elements ``h_{r, r+1}, h_{r+1, r+2},`` ``\dots , h_{n-1, n}`` and deleting the last column of the obtained matrix ``\tilde H``:

```math
    \tilde H = H R = \left( \begin{array}{ccc}
  L_{11} \\
  L_{21} &  \tilde L_{22} & 0
    \end{array}
     \right) =  \Bigl( \tilde L \ 0 \Bigr) \ ,
```
where ``R=R_0 R_1 \dots R_{n - 1 - r}`` is the matrix of Givens rotations, ``L_{11}`` and ``\tilde L_{22}\, - `` lower triangle matrices, and

```math
\begin{aligned}
&  \tilde L = \left( \begin{array}{cc}
  L_{11}    \\
  L_{21} &  \tilde L_{22}
    \end{array}  \right) \, ,  
\\
& \tilde H  \tilde H^T = H R R^T H^T = H H^T = \tilde A \, , 
\\
 & \tilde H  \tilde H^T = \tilde L \tilde L^T \, , \quad \tilde L \tilde L^T = \tilde A \, .
\end{aligned}
```
Orthogonal matrix ``R_t \ (0 \le t \le n - 1 -r), R_t \in R^{n \times n}`` which defines the Givens rotation annulating ``(r+t , r+t+1)``th element of the ``H_k R_0 \dots R_{t-1}`` matrix is defined as

```math
    R_t = \left( \begin{array}{cccccc}
   1      & \ldots & 0      & 0      & \ldots & 0 \\
   \vdots & \ldots & \ldots & \ldots & \ldots & \vdots \\
   0      & \ldots & c      & s      & \ldots & 0 \\
   0      & \ldots & -s     & c      & \ldots & 0 \\
   \vdots & \ldots & \ldots & \ldots & \ldots & \vdots \\
   0      & \ldots & 0      & 0      & \ldots & 1
    \end{array}
     \right) \ .
```
where entries ``( r+t , r+t )`` and ``( r+t+1 , r+t+1 )`` equal ``c``,  ``( r+t , r+t+1 )`` entry equals ``s``, and ``( r+t+1 , r+t )`` entry equals ``-s``, here ``c^2 + s^2 = 1``. Let's ``\tilde l_{ij}^{t-1} - `` coefficients of matrix ``H_k R_0 \dots R_{t-1}`` (``\tilde l_{ij}^{-1} - `` coefficients of ``H_k``) and

```math
  c = \frac{ \tilde l_{r+t,r+t}^{t-1} }
{ \sqrt{ (\tilde l_{r+t,r+t}^{t-1})^2 + (\tilde l_{r+t,r+t+1}^{t-1})^2 } } \ ,
\quad
  s = \frac{ \tilde l_{r+t,r+t+1}^{t-1} }
{ \sqrt{ (\tilde l_{r+t,r+t}^{t-1})^2 + (\tilde l_{r+t,r+t+1}^{t-1})^2 } } \ ,
```

Then matrix ``H_k R_0 \dots R_t`` will differ from ``H_k R_0 \dots R_{t-1}`` with entries of ``( r+t )`` и ``( r+t+1 )`` columns only, thereby
```math
\begin{aligned}
&  \tilde l_{i,r+t}^t = \tilde l_{i,r+t}^{t-1} = 0 \, , \ \ 
     \tilde l_{i,r+t+1}^t = \tilde l_{i,r+t+1}^{t-1} = 0 \, ,
     \qquad \qquad \qquad \ \   1 \le i \le r + t  - 1  \, ,
\\
&  \tilde l_{i,r+t}^t = c \tilde l_{i,r+t}^{t-1} +  s \tilde l_{i,r+t+1}^{t-1} \, ,
 \  \ 
\tilde l_{i,r+t+1}^t = -s \tilde l_{i,r+t}^{t-1} +  c \tilde l_{i,r+t+1}^{t-1} \,
     \quad  r + t \le i \le n - 1 \, ,
\end{aligned}
```
where

```math
  \tilde l_{r+t,r+t}^t =
\sqrt{ (\tilde l_{r+t,r+t}^{t-1})^2 + (\tilde l_{r+t,r+t+1}^{t-1})^2 } \ ,
\qquad
  \tilde l_{r+t,r+t+1}^t = 0   \,  .
```
Also, coefficient ``\tilde l_{r+t,r+t}^t`` is a nonnegative one.
In order to avoid unnecessary overflow or underflow during computation of ``c`` and ``s``, it was recommended [3] to calculate value of the square root ``w = \sqrt{x^2 + y^2}`` as folows:

```math
\begin{aligned}
&  v = \max \{ |x| , |y| \} \ , \qquad u = \min \{ |x| , |y| \} \ ,
\\
& w = \begin{cases} v \sqrt{ 1 + \left( \frac{u}{v} \right)^2 } \ ,
                    \quad v \ne 0 \cr \cr
              0 \ , \quad v = 0 \cr
      \end{cases} \ .

\end{aligned}
```
Applying the this technique to Cholesky factorization allows significally reduce the complexity of calculations. So, in case of adding a row and the symmetric column to the original matrix it will be necessary to carry out about ``n^2`` flops instead of about ``\frac {(n + 1) ^ 3} {3}`` flops for the direct calculation of the new Cholesky factor. In the case of deleting a row and the symmetric column from the original matrix, the new Cholesky factor can be obtained with about ``3(n - r)^2`` flops (the worst case requires about ``3 (n - 1) ^ 2`` operations) instead of about ``\frac {(n - 1) ^ 3} {3} `` flops required for its direct calculation.

## References

[1] T. F. Coleman, L. A. Hulbert, A direct active set  algorithm for large sparse quadratic programs with simple bounds, Mathematical Programming, 45, 1989.

[2] J. A. George, J. W-H. Liu, Computer Solution of Large Sparse Positive Definite Systems, Prentice-Hall, 1981.

[3] C. L. Lawson, R. J. Hanson, Solving least squares problems, SIAM, 1995.



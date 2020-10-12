# Quadratic Programming in Hilbert space. Part I.

Here an algorithm for finding a normal solution of linear inequalities system in Hilbert space is described. An original version of this algorithm was developed by V. Gorbunov and published in [3]. A modified version of that algorithm was presented in [6].

Let ``H`` be a Hilbert space with inner product ``{\langle \cdot \, , \cdot \rangle}_H`` and the induced norm ``\| \cdot \|_H``, there are elements ``h_i \in H``, numbers ``u_i, \, 1 \le i \le M+L`` and positive numbers  ``\delta _i, \, 1 \le i \le M``.
It is required to minimize functional ``J``:

```math
\tag{1}
   (J, \varphi) =  {\| \varphi \|}_H ^2 \to \min \ ,
```
subject to constraints

```math
\tag{2}
    {\langle h_i , \varphi \rangle}_H = u_i \, ,
    \qquad\qquad\qquad 1 \le i \le L \ ,
```

```math
\tag{3}
    | {\langle h_{i+L} , \varphi \rangle}_H - u_{i+L} | \le \delta _i \, ,
    \quad \ \ 1 \le i \le M \, , \ \varphi \in H.
```
The elements ``h_i, \, 1 \le i \le M+L`` are assumed to be linearly independent (thereby the system (2), (3) is compatible). Solution of the problem (1)—(3) obviously exists and it is unique as a solution of the problem of finding a projection of zero element onto a nonempty convex closed set in Hilbert space [1, 5]. Expanding the modules in (3) we rewrite the system (2), (3) in form:

```math
\tag{4}
    {\langle h_i , \varphi \rangle}_H = b_i \,\quad 1 \le i \le L \ ,
```
```math
\tag{5}
\qquad  {\langle h_i , \varphi \rangle}_H \le b_i \,\quad L+1 \le i \le N \ ,
```
where:

```math
\begin{aligned}
&  N = L + 2M \, ,
\\ 
& S = L + M \,  ,
\\ 
 &  b_i = u_i \,  , \quad  1 \le i \le L \, ,
\\ 
 &  b_{i+L} = u_{i+L} + \delta _i \  , 
 \\ 
 &  b_{i+S} = -u_{i+L} + \delta _i \  , 
 \\
 &  h_{i+S} = -h_{i+L} \ , \quad  1 \leq i \leq M \ .
\end{aligned}
```
Let ``\Phi`` be a convex set defined by constraints (4) and (5). We denote

```math
\begin{aligned}
 & I_1 = \lbrace 1, \dots , L \rbrace \,  , \quad I_2 = \lbrace L+1, \dots , N \rbrace \,  ,
\\
 &  I = I_1 \cup I_2 = \lbrace 1, \dots , N \rbrace \, ,
\\
 &  A(\varphi) = \lbrace i \in I : \langle h_i , \varphi \rangle_H = b_i \rbrace \, ,
 \  P(\varphi) = I \setminus A(\varphi) \, ,
\\
 &  \Psi(\varphi) = \{ \psi \in H : \langle h_i, \psi \rangle_H = b_i \, , \ i \in A(\varphi) \} \ ,
\\
 &  g_{ij} = \langle h_i , h_j \rangle_H \, , \quad 1 \le i,j \le S \ .
\end{aligned}
```
Feasible region of the ``\Psi(\varphi), \, \varphi \in \Phi`` is a face of set ``\Phi`` containing ``\varphi``. If ``A(\varphi)`` is empty (it is possible when ``L =0``) then ``\Psi(\varphi) = H``. 

In accordance with optimality conditions [4, 5] the solution of the problem (1), (4), (5) can be represented as: 

```math
\begin{aligned}
\tag{6}
& \qquad \sigma = \sum _{i=1} ^L \mu _i h_i + \sum _{i=L+1} ^{M+L} (\mu _i - \mu _{i+M}) h_i \ ,
\\ 
  & \qquad  \mu_i \le 0 \ , \quad \mu_{i+M} \le 0 \ , \quad  L+1 \le i \le L+M \, ,
\end{aligned}
```
```math
\tag{7}
 \mu_i (\langle h_i , \sigma \rangle_H - b_i ) = 0 \, ,   \quad  L+1 \le i \le N \, ,
```
```math
\tag{8}
 \mu_i \, \mu_{i+M} = 0 \ , \qquad \qquad \ \ \  L+1 \le i \le S \, .
```
Here complementary slackness conditions (7) means that ``\mu_k = 0`` for ``k \in P(\sigma)`` and the relations (8) reflect the fact that any pair of constraints (5) with indices ``i`` and ``i + M`` cannot be simultaneously active on the solution. Thereby there are at most ``S`` non-trivial coefficients ``\mu_i`` in formula (6) and actual dimension of the problem (1), (4), (5) is ``S``.

Let ``\hat\sigma`` be the normal solution of the system (4), then it can be written as

```math
    \hat\sigma =  \sum _{i=1} ^L \mu _i h_i \ ,
```
where coefficients ``\mu_i`` are determined from the system of linear equations with symmetric positive definite Gram matrix ``\{g_{ij}\}`` of the linear independent elements ``h_i``:

```math
    \sum _{j=1} ^L g_{ij} \mu _j = b_i  \ , \quad 1 \le i \le L \ .
```
If ``L = 0`` then there are no constraints (4) and ``\hat\sigma = 0``. If ``\hat\sigma`` satisfies constraints (5), that is, the inequalities

```math
    {\langle h_i , \hat\sigma \rangle}_H \le b_i \  ,  \qquad L+1 \le i \le N \ ,
```
which can be written like this:
```math
    \sum _{j=1} ^L g_{ij} \mu _j \le b_i  \ ,    \qquad L+1 \le i \le N \ ,
```
then ``\hat\sigma`` is a solution to the original problem (1), (4), (5), i.e. ``\sigma = \hat\sigma``.

Consider nontrivial case when ``\hat\sigma`` does not satisfy constraints (5). In this case ``\sigma`` belongs to the boundary of the set ``\Phi`` and it is a projection of zero element of space ``H`` onto the set ``\Psi (\sigma)``.
Projection ``\vartheta`` of zero element of space ``H`` onto the set ``\Psi (\varphi)`` can be presented as

```math
  \vartheta  =  \sum _{i \in A(\varphi)} \mu _i h_i \ ,
```
where factors ``\mu_i`` are defined from the system of linear equations with symmetric positive definite matrix

```math
    \sum _{j \in A(\varphi)} g_{ij} \mu _j = b_i  \ ,    \qquad i \in A(\varphi) \ ,
```
If factors ``\mu_i``, ``i \in I_2 \cap A(\varphi)`` corresponding to the inequality constraints are nonpositive, then we can set ``\mu_i = 0`` for ``i \in P(\varphi)`` and get all conditions (8)—(10) satisfied, thereby ``\vartheta  = \sigma`` is a solution of the problem under consideration. The following algorithm is based on this remark.  Let's describe its iteration.

Let's ``\sigma^ k`` be a feasible point of the system (4), (5):

```math
\tag{9}
  \sigma^k = \sum _{i = 1}^N \mu_i^k h_i \ ,
```
where there are at most ``S`` non-zero multipliers ``\mu_i^k``.  A projection of zero element of space ``H`` onto the set ``\Psi (\vartheta^k)`` we denote as ``\vartheta^k``,  and ``A_k = A(\sigma^k)``, ``P_k = P(\sigma^k) = I \setminus A_k``. Then ``\vartheta^k`` can be presented as

```math
 \qquad\qquad\qquad\qquad\qquad \vartheta^k = \sum _{i \in A_k} \lambda_i^k h_i = \sum _{i = 1}^N \lambda_i^k h_i \ , \quad \lambda_i^k = 0  \, , \  \forall i \in P_k \ ,
```
```math
\tag{10}
  \sum _{j \in A_k} g_{ij} \lambda_j^k = b_i  \, ,  \quad i \in A_k \ .
```
There are two possible cases: the element ``\vartheta^k`` is feasible or it is not feasible. In the first case we check the optimality conditions, namely: if ``\lambda_i^k \le 0, \  \forall i \in I_2`` then  ``\vartheta^k`` is the solution of the problem. If the optimality conditions are not satisfied then we set ``\sigma^{k+1} = \vartheta^k``, find an index ``i \in A_k`` such that ``\lambda_i^k > 0`` and remove it from ``A_k``.
In the second case ``\sigma^{k+1}`` will be defined as a feasible point of the ray ``\vartheta (t)``

```math
\tag{11}
  \vartheta (t) = \sigma^k + t (\vartheta^k - \sigma^k)  \ ,
```math
such that it is closest to the ``\vartheta^k``. Denote ``t^k_{min}`` the corresponding value of ``t`` and ``i_k \in P_k`` — related number of the violated constraint. This index ``i_k`` will be added to ``A_k`` forming ``A_{k+1}``.
Since all ``\sigma^k`` are feasible points, the minimization process proceeds within the feasible region and value of ``\| \sigma^k \|_H`` is not increasing. The minimization process is finite because of linear independence of the constraints, it eliminates the possibility of the algorithm cycling. The feasibility of the ``\vartheta^k`` can be checked as follows. Introduce the values ``{e_i^k}, \,  k \in P_k``:

```math
 e_i^k = \langle h_i , \vartheta^k - \sigma^k \rangle_H =  \sum _{j=1}^N g_{ij}(\lambda_j^k - \mu _j^k) \ ,  \quad i \in P_k \ ,
```
if ``{e_i^k} > 0`` then constraint with number ``i`` can be violated at the transition from ``\sigma^k`` to point ``\vartheta^k``, in a case when ``{e_i^k} < 0`` the constraint with number ``i + M`` can be violated. For all ``i \in P_k`` compute values

```math
  t_i^k = \begin{cases}
           \frac {b_i - \langle h_i , \sigma^k \rangle_H}{e_i^k} \ ,
                          \quad e_i^k > 0  \cr \cr
           \frac {-b_{i+M} - \langle h_{i+M} , \sigma^k \rangle_H}{e_i^k} \ ,
                          \quad e_i^k < 0  \cr \cr
           1 \ , \quad e_i^k = 0 \cr
               \end{cases} \ ,
```
where

```math
\langle h_i , \sigma^k \rangle_H =   \sum _{j=1}^N g_{ij} \mu _j^k \ ,
\qquad 
 \langle h_{i+M} , \sigma^k \rangle_H =
                 \sum _{j=1}^N g_{i+M,j} \mu _j^k \ , \quad i \in P_k \ ,
```
here all `` t_i^k \ge 0``. Now the maximum feasible step ``t^k_{min}`` is computed as

```math
    t^k_{min} = \min_{i \in P_k} \{ t_i^k \} 
```
and ``\vartheta^k`` is feasible if and only if ``t^k_{min} \ge 1`` (see (11)).

Thereby the algorithm's iteration consist of seven steps:

1. Initialization. Let ``\sigma^0`` be a feasible point of the system (4), (5) and ``\mu_i^0``  are corresponding multipliers (9).
2. Compute multipliers ``\lambda_j^k, \ i \in A_k`` as solution of the system (10).
3. If ``| A_k | = S`` then go to Step 6.
4. Compute ``t_i^k, \ \forall i \in P_k``. Find ``t^k_{min}`` and the corresponding index ``i_k``. 
5. If ``t^k_{min} < 1`` (projection ``\vartheta^k`` is not feasible) then set

```math
\begin{aligned}
& \mu_i^{k+1} = \mu_i^k + t^k_{min} (\lambda_i^k - \mu_i^k) \ , \quad i \in A_k \ ,
\qquad \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad \quad 
\\ 
& A_{k+1} = A_k \cup \{ i_k \} \ , 
\\
 \text{and return to Step 1.} & \, 
\end{aligned}
``` 
6. Projection ``\vartheta^k`` is feasible. If exists index ``i_p, \ i_p \in A_k`` such that ``\lambda_{i_p}^k > 0`` then set

```math
\begin{aligned}
 & A_{k+1} = A_k \setminus \{ i_p \} \ , 
 \qquad \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad \quad
\\
 & \mu_i^{k+1} = \lambda_i^k  \ , \quad  i \in A_{k+1} \ , 
 \\
 \text{and return to Step 1.} & \, 
\end{aligned}
``` 
7. Set ``\sigma = \vartheta^k``. Stop.

The algorithm starts from an initial feasible point of the system (4), (5). Such point ``\sigma^0`` can be defined as the normal solution of the system

```math
 {\langle h_i , \varphi \rangle}_H = u_i \, \quad 1 \le i \le S \ ,
```
and can be presented as

```math
  \sigma^0 = \sum _{i = 1}^S \mu_i^0 h_i \ ,
```
where multipliers ``\mu_i^0`` are defined from the system

```math
  \sum _{j=1} ^S g_{ij} \mu _j^0  = u_i  \ ,    \qquad 1 \le i \le  S \ .
```
The algorithm described here implements a variant of the active set method [2] with the simplest form of the minimized functional (1). It can also be treated as a variant of the conjugate gradient method [7] and allows to solve the considered problem by a finite number of operations. The conjugate gradient method in this case is reduced into the problem of finding the projection ``\vartheta^ k``, that is, into solving of system (10). Infinite dimensionality of the space ``H`` for this algorithm is irrelevant.
The most computationally costly part of the algorithm is calculation of projection ``\vartheta^k`` by solving system (10). Matrices of all systems (10) are positive definite because of linear independence of elements ``h_i``, and it is convenient to use Cholesky decomposition to factorize them . At every step of the algorithm a constraint is added to or removed from the current active set and the system (10) matrix gets modified by the corresponding row and symmetric column addition or deletion. It is possible to get Cholesky factor of the modified matrix without computing the full matrix factorization, it allows greatly reduce the total amount of computations.

This algorithm can be also applied in a case when instead of modular constraints (3) there are one-sided inequality constraints

```math
    {\langle h_{i+L} , \varphi \rangle}_H \le u_{i+L} \, , \quad 1 \le i \le M \, ,
```
for that it is enough to set 

```math
  b_{i+L} = u_{i+L} \ , \quad  b_{i+S} = -\infty \, , \quad 1 \le i \le M \, .
```

**References**

[1] A. Balakrishnan, Applied Functional Analysis, Springer-Verlag, New York, 1976.

[2]  P. Gill, W. Murray, M. Wright, Practical Optimization, Academic Press, London, 1981.

[3] V. Gorbunov, The method of reduction of unstable computational problems (Metody reduktsii neustoichivykh vychislitel'nkh zadach), Ilim, 1984.

[4] A. Ioffe, V. Tikhomirov, Theory of extremal problems, North-Holland, Amsterdam, 1979.

[5] P.-J. Laurent, Approximation et optimization, Paris, 1972.

[6] I. Kohanovsky, Data approximation using multidimensional normal splines. Unpublished manuscript. 1996.

[7] B. Pshenichnyi, Yu. Danilin, Numerical methods in extremal problems,  Mir Publishers, Moscow, 1978.

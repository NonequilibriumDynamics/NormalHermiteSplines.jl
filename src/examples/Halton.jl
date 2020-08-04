### References:
# 1. M. Berblinger, C. Schlier, Monte Carlo integration with quasi-random numbers: Some experience, Computer Physics Communications 66, 157-166, 1991.
# 2. M. Kolář, S. O'Shea, Fast, portable, and reliable algorithm for the calculation of Halton numbers Computers & Mathematics with Applications 25(7), 1993.
# 3. Halton sequence (https://en.wikipedia.org/wiki/Halton_sequence)
# 4. R. Renka, R. Brown, Algorithm 792: accuracy test of ACM algorithms for interpolation of scattered data in the plane, ACM Transactions on Mathematical Software (TOMS) 25 (1),1999
# 5. S. De Marchi, W. Erb, F. Marchetti, E. Perracchione, M. Rossini, Shape-Driven Interpolation with Discontinuous Kernels: Error Analysis, Edge Extraction and Applications in MPI, arXiv: Numerical Analysis, 2019.
# 6. S. De Marchi, A. Martinez, E. Perracchione, Fast and stable rational RBF-based partition of unity interpolation. Journal of Computational and Applied Mathematics, 2018.
# 7. [W. Erb, C. Kaethner, M. Ahlborg, T.M. Buzug, Bivariate Lagrange interpolation at the node points of non-degenerate Lissajous nodes, Numer. Math. 133, 1, 2016.](https://www.researchgate.net/publication/268988454_Bivariate_Lagrange_interpolation_at_the_node_points_of_non-degenerate_Lissajous_curves)

###

_b = 0
_mn = 0.0
_dn = 1.0

"""
Initiate the Halton sequence calculation.
This function must be called prior to `get_halton()` which will
then start with n0-th Halton number.
"""
function init_halton(base::Int, n0::Int = 66)
 global _b = base
 while n0 > 1
     n0 -= 1
     get_halton_node()
 end
end

function reset_halton()
 global _b = 0
 global _mn = 0.0
 global _dn = 1.0
end

"""
Return next Halton number.
"""
function get_halton_node()
     x = _dn - _mn
     if x == 1.0
         global _mn = 1.0
         global _dn *= _b
     else
         y = _dn / _b
         while x < y
             y /= _b
         end
         global _mn = (1.0 + _b) * y - x
     end
     return _mn / _dn
end

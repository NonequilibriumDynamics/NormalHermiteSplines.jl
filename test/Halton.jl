### References:
# 1. M. Berblinger, C. Schlier, "Monte Carlo integration with quasi-random numbers: Some experience", Computer Physics Communications 66, 157-166, 1991.
# 2. M. Kolář, S. O'Shea, "Fast, portable, and reliable algorithm for the calculation of Halton numbers" Computers & Mathematics with Applications 25(7), 1993.
# 3. Halton sequence (https://en.wikipedia.org/wiki/Halton_sequence)
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
     get_halton()
 end
end

"""
Return next Halton number.
"""
function get_halton()
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

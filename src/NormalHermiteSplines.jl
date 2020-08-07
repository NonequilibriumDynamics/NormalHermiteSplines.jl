module NormalHermiteSplines

#### Inteface deinition
export prepare, construct, interpolate
export evaluate, evaluate_one, evaluate_gradient
export NormalSpline, RK_H0, RK_H1, RK_H2
export get_epsilon, estimate_epsilon, get_cond, estimate_interpolation_quality
# -- 1D case --
export evaluate_derivative
# --
####

#include("./Examples/Main.jl")

using LinearAlgebra

abstract type ReproducingKernel end
abstract type ReproducingKernel_0 <: ReproducingKernel end
abstract type ReproducingKernel_1 <: ReproducingKernel_0 end
abstract type ReproducingKernel_2 <: ReproducingKernel_1 end

abstract type AbstractSpline end

@doc raw"
`struct NormalSpline{T, RK} <: AbstractSpline where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Define a structure containing full information of a normal spline
# Fields
- `_kernel`: a reproducing kernel spline was built with.
- `_compression`: factor of transforming the original node locations into unit hypercube
- `_nodes`: transformed function value nodes
- `_values`: function values at interpolation nodes
- `_d_nodes`: transformed function derivative nodes
- `_es`: normalized derivative directions
- `_d_values`: function directional derivative values
- `_min_bound`: minimal bounds of the original node locations area
- `_gram`: Gram matrix of the problem
- `_chol`: Cholesky factorization of the Gram matrix
- `_mu`: spline coefficients
- `_cond`: estimation of the Gram matrix condition number
"
struct NormalSpline{T, RK} <: AbstractSpline where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    _kernel::RK
    _compression::T
    _nodes::Union{Matrix{T}, Nothing}
    _values::Union{Vector{T}, Nothing}
    _d_nodes::Union{Matrix{T}, Nothing}
    _es::Union{Matrix{T}, Nothing}
    _d_values::Union{Vector{T}, Nothing}
    _min_bound::Union{Vector{T}, Nothing}
    _gram::Union{Matrix{T}, Nothing}
    _chol::Union{Cholesky{T, Matrix{T}}, Nothing}
    _mu::Union{Vector{T}, Nothing}
    _cond::T
end

include("./ReproducingKernels.jl")
include("./GramMatrix.jl")
include("./Interpolate.jl")

"""
`prepare(nodes::Matrix{T}, kernel::RK = RK_H0()) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Prepare the normal spline by constructing and factorizing Gram matrix of the interpolation problem.
Initialize the `NormalSpline` object.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space and
           `n_1` is the number of function value nodes. It means that each column in the matrix defines one node.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H0` if the spline is constructing as a continuous function,
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the partly initialized `NormalSpline` object that must be passed to `construct` function
        in order to complete the spline initialization.
"""
function prepare(nodes::Matrix{T},
                 kernel::RK = RK_H0()
                ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
     spline = _prepare(nodes, kernel)
     return spline
end

"""
`construct(spline::NormalSpline{T, RK}, values::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Construct an interpolating spline by calculating coefficients of the spline and completely initializing
the `NormalSpline` object.
# Arguments
- `spline`: the partly initialized `NormalSpline` object returned by `prepare` function.
- `values`: function values at `nodes` nodes.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function construct(spline::NormalSpline{T, RK},
                   values::Vector{T}
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    spline = _construct(spline, values)
    return spline
end

"""
`interpolate(nodes::Matrix{T}, values::Vector{T}, kernel::RK = RK_H0()) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Create the normal spline by `values` of function defined at `nodes`.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space
           and `n_1` is the number of function value nodes.
           It means that each column in the matrix defines one node.
- `values`: function values at `nodes` nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H0` if the spline is constructing as a continuous function,
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function interpolate(nodes::Matrix{T},
                     values::Vector{T},
                     kernel::RK = RK_H0()
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
     spline = _prepare(nodes, kernel)
     spline = _construct(spline, values)
     return spline
end

"""
`evaluate(spline::NormalSpline{T, RK}, points::Matrix{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Evaluate normal spline at the locations defined in `points`.

# Arguments
- `spline: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `points`: locations at which spline values are evaluating.
            This should be an `n×m` matrix, where `n` is dimension of the sampled space
            and `m` is the number of locations where spline values are evaluating.
            It means that each column in the matrix defines one location.

Return: `Vector{T}` of spline values at the locations defined in `points`.
"""
function evaluate(spline::NormalSpline{T, RK},
                  points::Matrix{T}
                 ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    return _evaluate(spline, points)
end

"""
`evaluate_one(spline::NormalSpline{T, RK}, point::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Evaluate normal spline at the `point` location.

# Arguments
- `spline`: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `point`: location at which spline value is evaluating.
           This should be a vector of size `n`, where `n` is dimension of the sampled space.

Return: spline value at the location defined in `point`.
"""
function evaluate_one(spline::NormalSpline{T, RK},
                      point::Vector{T}
                     ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    return _evaluate(spline, reshape(point, length(point), 1))[1]
end

"""
`evaluate_gradient(spline::NormalSpline{T, RK}, point::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Evaluate gradient of the normal spline at the location defined in `point`.

# Arguments
- `spline`: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `point`: location at which gradient value is evaluating.
           This should be a vector of size `n`, where `n` is dimension of the sampled space.

Return: `Vector{T}` - gradient of the normal spline at the location defined in `point`.
"""
function evaluate_gradient(spline::NormalSpline{T, RK},
                           point::Vector{T}
                          ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    return _evaluate_gradient(spline, point)
end

########

"""
`prepare(nodes::Matrix{T}, d_nodes::Matrix{T}, es::Matrix{T}, kernel::RK = RK_H1()) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Prepare the normal spline by constructing and factorizing Gram matrix of the interpolation problem.
Initialize the `NormalSpline` object.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space and
           `n_1` is the number of function value nodes.
            It means that each column in the matrix defines one node.
- `d_nodes`: The function directional derivatives nodes.
             This should be an `n×n_2` matrix, where `n` is dimension of the sampled space and
             `n_2` is the number of function directional derivative nodes.
- `es`: Directions of the function directional derivatives.
        This should be an `n×n_2` matrix, where `n` is dimension of the sampled space and
        `n_2` is the number of function directional derivative nodes.
        It means that each column in the matrix defines one direction of the function directional derivative.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the partly initialized `NormalSpline` object that must be passed to `construct` function
        in order to complete the spline initialization.
"""
function prepare(nodes::Matrix{T},
                 d_nodes::Matrix{T},
                 es::Matrix{T},
                 kernel::RK = RK_H1()
                ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
     spline = _prepare(nodes, d_nodes, es, kernel)
     return spline
end

"""
`construct(spline::NormalSpline{T, RK}, values::Vector{T}, d_values::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Construct an interpolating spline by calculating coefficients of the spline and completely initializing
the `NormalSpline` object.
# Arguments
- `spline`: the partly initialized `NormalSpline` object returned by `prepare` function.
- `values`: function values at `nodes` nodes.
- `d_values`: function directional derivative values at `d_nodes` nodes.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function construct(spline::NormalSpline{T, RK},
                   values::Vector{T},
                   d_values::Vector{T}
                  ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
     spline = _construct(spline, values, d_values)
     return spline
end

"""
`interpolate(nodes::Matrix{T}, values::Vector{T}, d_nodes::Matrix{T}, es::Matrix{T}, d_values::Vector{T}, kernel::RK = RK_H1()) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Create the normal spline by `values` of function defined at `nodes` and
`d_values` of function directional derivatives defined at `d_nodes`.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space
           and `n_1` is the number of function value nodes.
           It means that each column in the matrix defines one node.
- `values`: function values at `nodes` nodes.
- `d_nodes`: The function directional derivative nodes.
            This should be an `n×n_2` matrix, where `n` is dimension of the sampled space and
            `n_2` is the number of function directional derivative nodes.
- `es`: Directions of the function directional derivatives.
       This should be an `n×n_2` matrix, where `n` is dimension of the sampled space and
       `n_2` is the number of function directional derivative nodes.
       It means that each column in the matrix defines one direction of the function directional derivative.
- `d_values`: function directional derivative values at `d_nodes` nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function interpolate(nodes::Matrix{T},
                     values::Vector{T},
                     d_nodes::Matrix{T},
                     es::Matrix{T},
                     d_values::Vector{T},
                     kernel::RK = RK_H1()
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
     spline = _prepare(nodes, d_nodes, es, kernel)
     spline = _construct(spline, values, d_values)
     return spline
end

"""
`get_epsilon(spline::NormalSpline{T, RK}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Get a `ε` - 'scaling parameter' of Bessel Potential space the normal spline was built in.
# Arguments
- `spline`: the `NormalSpline` object returned by `prepare`, `construct` or `interpolate` function.

Return: `ε`.
"""
function get_epsilon(spline::NormalSpline{T, RK}
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    return spline._kernel.ε
end

"""
`estimate_epsilon(nodes::Matrix{T}) where T <: AbstractFloat`

Get an the estimation of the 'scaling parameter' of Bessel Potential space the normal spline was built in.
This coinsides with result returned by `get_epsilon` function if all `nodes` are located in a unit hypercube,
othewise the estimated value of `ε` can be significantly higher then necessary.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space
           and `n_1` is the number of function value nodes.
           It means that each column in the matrix defines one node.

Return: estimation of `ε`.
"""
function estimate_epsilon(nodes::Matrix{T}
                         ) where T <: AbstractFloat
    ε = _estimate_ε(nodes)
    return ε
end

"""
`estimate_epsilon(nodes::Matrix{T}, d_nodes::Matrix{T}) where T <: AbstractFloat`

Get an the estimation of the 'scaling parameter' of Bessel Potential space the normal spline was built in.
This coinsides with result returned by `get_epsilon` function if all `nodes` and `d_nodes` are located in a unit hypercube,
othewise the estimated value of `ε` can be significantly higher then necessary.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space
           and `n_1` is the number of function value nodes.
           It means that each column in the matrix defines one node.
- `d_nodes`: The function directional derivative nodes.
           This should be an `n×n_2` matrix, where `n` is dimension of the sampled space and
           `n_2` is the number of function directional derivative nodes.

Return: estimation of `ε`.
"""
function estimate_epsilon(nodes::Matrix{T},
                          d_nodes::Matrix{T}
                         ) where T <: AbstractFloat
    ε = _estimate_ε(nodes, d_nodes)
    return ε
end

"""
`get_cond(spline::NormalSpline{T, RK}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Get an estimation of the Gram matrix condition number.
(C. Brás, W. Hager, J. Júdice, An investigation of feasible descent algorithms for estimating the condition number of a matrix. TOP Vol.20, No.3, 2012.)
# Arguments
- `spline`: the `NormalSpline` object returned by `prepare`, `construct` or `interpolate` function.

Return: an estimation of the Gram matrix condition number.
"""
function get_cond(spline::NormalSpline{T, RK}
                 ) where {T <: AbstractFloat, RK <: ReproducingKernel}
    return spline._cond
end

# Return the Root Mean Square Error (RMSE) of interpolation
@inline function get_RMSE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return norm(f .- σ) / sqrt(length(f))
end

# Return the Maximum Absolute Error (MAE) of interpolation
@inline function get_MAE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return maximum(abs.(f .- σ))
end

"""
`estimate_interpolation_quality(spline::NormalSpline{T, RK}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Estimate the interpolation quality
# Arguments
- `spline`: the `NormalSpline` object returned by `construct` or `interpolate` function.

Return: RMSE of interpolation at function value nodes.
"""
function estimate_interpolation_quality(spline::NormalSpline{T, RK}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    σ = evaluate(spline, spline._nodes)
    return get_RMSE(σ, spline._values)
end

############################## One-dimensional case

"""
`prepare(nodes::Vector{T}, kernel::RK = RK_H0()) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Prepare 1D normal spline by constructing and factorizing Gram matrix of the interpolation problem.
Initialize the `NormalSpline` object.
# Arguments
- `nodes`: The function value nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H0` if the spline is constructing as a continuous function,
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the partly initialized `NormalSpline` object that must be passed to `construct` function
        in order to complete the spline initialization.
"""
function prepare(nodes::Vector{T},
                 kernel::RK = RK_H0()
                ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
     spline = _prepare(Matrix(nodes'), kernel)
     return spline
end

"""
`interpolate(nodes::Vector{T}, values::Vector{T}, kernel::RK = RK_H0()) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Create 1D normal spline by `values` of function defined at `nodes`.
# Arguments
- `nodes`: The function value nodes.
           This should be an `n×n_1` matrix, where `n` is dimension of the sampled space
           and `n_1` is the number of function value nodes.
           It means that each column in the matrix defines one node.
- `values`: function values at `n_1` interpolation nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H0` if the spline is constructing as a continuous function,
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function interpolate(nodes::Vector{T},
                     values::Vector{T},
                     kernel::RK = RK_H0()
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
     spline = _prepare(Matrix(nodes'), kernel)
     spline = _construct(spline, values)
     return spline
end

"""
`evaluate(spline::NormalSpline{T, RK}, points::Vector{T}) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Evaluate 1D normal spline at the `points` locations.

# Arguments
- `spline`: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `points`: locations at which spline values are evaluating.
                       This should be a vector of size `m` where `m` is the number of evaluating points.

Return: spline value at the `point` location.
"""
function evaluate(spline::NormalSpline{T, RK},
                  points::Vector{T}
                 ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    return _evaluate(spline, Matrix(points'))
end

"""
`evaluate_one(spline::NormalSpline{T, RK}, point::T) where {T <: AbstractFloat, RK <: ReproducingKernel_0}`

Evaluate 1D normal spline at the `point` location.

# Arguments
- `spline`: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `point`: location at which spline value is evaluating.

Return: spline value at the `point` location.
"""
function evaluate_one(spline::NormalSpline{T, RK},
                      point::T
                     ) where {T <: AbstractFloat, RK <: ReproducingKernel_0}
    v_points = Vector{T}(undef, 1)
    v_points[1] = point
    return _evaluate(spline, Matrix(v_points'))[1]
end

"""
`evaluate_derivative(spline::NormalSpline{T, RK}, point::T) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Evaluate 1D normal spline derivative at the `point` location.

# Arguments
- `spline`: the `NormalSpline` object returned by `interpolate` or `construct` function.
- `point`: location at which spline value is evaluating.

Return: spline derivative value at the `point` location.
"""
function evaluate_derivative(spline::NormalSpline{T, RK},
                             point::T
                            ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
    v_points = Vector{T}(undef, 1)
    v_points[1] = point
    return _evaluate_gradient(spline, v_points)[1]
end

"""
`prepare(nodes::Vector{T}, d_nodes::Vector{T}, kernel::RK = RK_H1()) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Prepare the 1D normal spline by constructing and factorizing Gram matrix of the interpolation problem.
Initialize the `NormalSpline` object.
# Arguments
- `nodes`: The function value nodes.
- `d_nodes`: The function derivatives nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the partly initialized `NormalSpline` object that must be passed to `construct` function
        in order to complete the spline initialization.
"""
function prepare(nodes::Vector{T},
                 d_nodes::Vector{T},
                 kernel::RK = RK_H1()
                ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}

     es = ones(T, size(d_nodes, 2))
     spline = _prepare(Matrix(nodes'), Matrix(d_nodes'), Matrix(es'), kernel)
     return spline
end

"""
`interpolate(nodes::Vector{T}, values::Vector{T}, d_nodes::Vector{T}, d_values::Vector{T}, kernel::RK = RK_H1()) where {T <: AbstractFloat, RK <: ReproducingKernel_1}`

Create the 1D normal spline by `values` of function defined at `nodes` and
`d_values` of function derivatives defined at `d_nodes`.
# Arguments
- `nodes`: The function value nodes.
- `values`: function values at `nodes` nodes.
- `d_nodes`: The function derivative nodes.
- `d_values`: function derivative values at `d_nodes` nodes.
- `kernel`: reproducing kernel of Bessel potential space the normal spline is constructed in.
            It must be a struct object of the following type:
              `RK_H1` if the spline is constructing as a differentiable function,
              `RK_H2` if the spline is constructing as a twice differentiable function.

Return: the completely initialized `NormalSpline` object that can be passed to `evaluate` function
        to interpolate the data to required points.
"""
function interpolate(nodes::Vector{T},
                     values::Vector{T},
                     d_nodes::Vector{T},
                     d_values::Vector{T},
                     kernel::RK = RK_H1()
                    ) where {T <: AbstractFloat, RK <: ReproducingKernel_1}
     es = ones(T, size(d_nodes, 2))
     spline = _prepare(Matrix(nodes'), Matrix(d_nodes'), Matrix(es'), kernel)
     spline = _construct(spline, values, d_values)
     return spline
end

end # module
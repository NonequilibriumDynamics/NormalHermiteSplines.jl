export test_1D, test_2D, test_3D
export init_halton, get_halton_node, get_Lissajous_nodes

export big_sur

export get_2D_grid, get_2D_grid2, get_2D_test1_nodes, get_2D_Lissajous_nodes, get_2D_halton_nodes, get_2D_uniformrandom_grid
export get_2D_rect_grid, get_2D_eps_grid, get_2D_border_nodes
export get_2D_model1, get_2D_model2, get_2D_model3, get_2D_model4, get_2D_model5, get_2D_model6, get_2D_model6_grad
export get_2D_model10, get_2D_model11, get_2D_model12, get_2D_model12_Grad
export get_2D_model13, get_2D_model13_Grad
export readme_1, readme_2
export _gram

export get_1D_model1, get_1D_model1_grad
export get_1D_grid, get_1D_eps_grid, get_1D_halton_nodes

export get_3D_model1, get_3D_model2, get_3D_model3, get_3D_model4, get_3D_model5, get_3D_model6
export get_3D_model1_grad, get_3D_model2_grad, get_3D_model3_grad, get_3D_model4_grad, get_3D_model5_grad, get_3D_model6_grad
export get_3D_grid, get_3D_eps_grid, get_3D_halton_nodes, get_3D_plot_grid, get_3D_random_grid

export get_separation_distance, get_fill_distance

export readme_3
export demo
export usage1, usage2

export param1, param10


# Return the Root Mean Square Error (RMSE) of interpolation
@inline function get_RMSE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return norm(f .- σ) / sqrt(length(f))
end

# Return the Maximum Absolute Error (MAE) of interpolation
@inline function get_MAE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return maximum(abs.(f .- σ))
end

# Return the Relative Root Mean Square Error (RRMSE) of interpolation
@inline function get_RRMSE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return norm(f .- σ) / (sqrt(length(f)) * maximum(abs.(f)))
end

# Return the Relative Maximum Absolute Error (RMAE) of interpolation
@inline function get_RMAE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
    return maximum(abs.(f .- σ)) / maximum(abs.(f))
end

# # Return the Relative Root Mean Square Error (RRMSE) of interpolation
# @inline function get_RRMSE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
#     n = length(f)
#     del = similar(f)
#     @inbounds for i = 1:n
#         if f[i] <= T(1.0)
#             del[i] = f[i] - σ[i]
#         else
#             del[i] = (f[i] - σ[i]) / f[i]
#         end
#     end
#     return norm(del) / sqrt(n)
# end
#
# # Return the Relative Maximum Absolute Error (RMAE) of interpolation
# @inline function get_RMAE(f::Vector{T}, σ::Vector{T}) where T <: AbstractFloat
#     n = length(f)
#     del = similar(f)
#     @inbounds for i =1:n
#         if f[i] <= T(1.0)
#             del[i] = abs(f[i] - σ[i])
#         else
#             del[i] = abs(f[i] - σ[i]) / abs(f[i])
#         end
#     end
#     return maximum(del)
# end

include("Halton.jl")
include("Grids.jl")
include("Models.jl")
include("Examples_1D.jl")
include("Examples_2D.jl")
include("Examples_3D.jl")
include("Demo.jl")

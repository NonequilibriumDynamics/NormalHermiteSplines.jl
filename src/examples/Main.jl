export test_2D
export init_halton, get_halton_node, get_Lissajous_nodes
export get_2D_grid, get_2D_grid2, get_2D_test1_nodes, get_2D_Lissajous_nodes, get_2D_halton_nodes
export get_2D_rect_grid, get_2D_eps_grid, get_2D_border_nodes
export get_2D_model1, get_2D_model2, get_2D_model3, get_2D_model4, get_2D_model5, get_2D_model6
export get_2D_test1_nodes, get_2D_Lissajous_nodes, get_2D_halton_nodes
export get_2D_model10, get_2D_model11, get_2D_model12, get_2D_model12_Grad
export get_2D_model13, get_2D_model13_Grad
export readme_1, readme_2
export _gram
export get_1D_model1, get_1D_model1_grad
export big_sur

#include("Demo.jl")
include("Halton.jl")
include("Grids.jl")
include("Models.jl")
include("Examples_1D.jl")
include("Examples_2D.jl")

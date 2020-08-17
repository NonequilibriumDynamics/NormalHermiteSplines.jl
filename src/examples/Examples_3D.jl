using Printf
using PyPlot

function test_3D(model_id::Int,
                 use_derivatives::Bool = false,
                 type_of_samples::Int = 1,
                 n_of_samples::Int = 4,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 plot_grid_size::Int = 50
#                 plot_grid_size::Int = 25
#                 ,do_parallel::Bool = false
                )
    if use_derivatives && type_of_kernel == 0
        error("Cannot use derivative data when type_of_kernel is `0` (`RK_H0` kernel)")
    end

    if type_of_samples == 1
        samples_size = [1, 100, 500, 1000, 4000, 8000, 16000]
        nodes = get_3D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 2
        samples_size = [1, 4, 7, 9, 15, 19, 25]
        nodes = get_3D_grid(samples_size[n_of_samples]) #1(8), 4(125), 7(512), 9(1000), 15(4096)(d), 19(8000), 24(15625)
    elseif type_of_samples == 3
        samples_size = [1, 4, 7, 9, 15, 19, 25]
        nodes = get_3D_eps_grid(samples_size[n_of_samples])
    else
        error("Incorrect value of 'type_of_samples'")
    end

    if type_of_kernel == 0
        if eps == 0.0
            rk = RK_H0()
        else
            rk = RK_H0(eps)
        end
    elseif type_of_kernel == 1
        if eps == 0.0
            rk = RK_H1()
        else
            rk = RK_H1(eps)
        end
    elseif type_of_kernel == 2
        if eps == 0.0
            rk = RK_H2()
        else
            rk = RK_H2(eps)
        end
    else
        error("Incorrect value of 'type_of_kernel'")
        return
    end

    n_1 = size(nodes, 2)
    #grid = get_3D_grid(plot_grid_size)
    grid = get_3D_plot_grid(plot_grid_size)
    m = size(grid, 2)
    if model_id == 1
        wnodes = similar(nodes)
        k = 0
        for i = 1:n_1
            if get_3D_model1(nodes[:, i]) < -0.1
               continue
            end
            k += 1
            wnodes[1,k] = nodes[1,i]
            wnodes[2,k] = nodes[2,i]
            wnodes[3,k] = nodes[3,i]
        end
        nodes = wnodes[:,1:k]

        wgrid =similar(grid)
        k = 0
        for i = 1:m
            if get_3D_model1(grid[:, i]) < -0.1
               continue
            end
            k += 1
            wgrid[1,k] = grid[1,i]
            wgrid[2,k] = grid[2,i]
            wgrid[3,k] = grid[3,i]
        end
        grid = wgrid[:,1:k]

        n_1 = size(nodes, 2)
        m = size(grid, 2)
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model1(nodes[:, i])
            k += 1
            grad = get_3D_model1_grad()
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[2]
            es[1,k] = 0.0
            es[2,k] = 1.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[3]
            es[1,k] = 0.0
            es[2,k] = 0.0
            es[3,k] = 1.0
        end

        for i = 1:m
            f[i] = get_3D_model1(grid[:, i])
        end
    elseif model_id == 2
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model2(nodes[:, i])
            k += 1
            grad = get_3D_model2_grad()
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[2]
            es[1,k] = 0.0
            es[2,k] = 1.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[3]
            es[1,k] = 0.0
            es[2,k] = 0.0
            es[3,k] = 1.0
        end

        for i = 1:m
            f[i] = get_3D_model2(grid[:, i])
        end
    elseif model_id == 3
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model3(nodes[:, i])
            k += 1
            grad = get_3D_model3_grad()
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[2]
            es[1,k] = 0.0
            es[2,k] = 1.0
            es[3,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            d_nodes[3,k] = nodes[3,i]
            du[k] = grad[3]
            es[1,k] = 0.0
            es[2,k] = 0.0
            es[3,k] = 1.0
        end

        for i = 1:m
            f[i] = get_3D_model3(grid[:, i])
        end
    else
        error("Incorrect value of 'model_id'")
        return
    end

    if use_derivatives
        @printf "nodes#: %d  d_nodes#: %d (total nodes: %d)\n" n_1 n_1 (n_1+n_1)
    else
        @printf "nodes#: %d\n" n_1
    end

    @printf "model_id type_of_samples  n_of_samples  type_of_kernel   plot_grid_size  use_derivatives\n"
    @printf "%2d      %2d             %4d             %1d               %3d               %s\n" model_id type_of_samples n_of_samples type_of_kernel  plot_grid_size use_derivatives
#
    @printf "Creating spline..\n"
    ts = time_ns()
    if use_derivatives
        if rk.ε == 0
             epsilon = estimate_epsilon(nodes, d_nodes)
             @printf "Estimated EPSILON:%0.1e\n" epsilon
        end
        spline = interpolate(nodes, u, d_nodes, es, du, rk)
    else
        if rk.ε == 0
             epsilon = estimate_epsilon(nodes)
             @printf "Estimated EPSILON:%0.1e\n" epsilon
        end
        spline = interpolate(nodes, u, rk)
    end
    te = time_ns()
    c_time = (te - ts) / 10^9
    @printf "Spline created. time: %0.1e sec\n" c_time
    cond = get_cond(spline)
    ε = get_epsilon(spline)
    @printf "EPSILON:%0.1e   COND: %0.1e \n" ε cond

    iq = assess_quality(spline)
    @printf "interpolation quality: %0.1e\n" iq

    @printf "Evaluating spline..\n"
    ts = time_ns()
    σ = evaluate(spline, grid)
#    σ = evaluate(spline, grid, do_parallel)
    te = time_ns()
    e_time = (te - ts) / 10^9
    @printf "Spline evaluated. time: %0.1e sec\n" e_time
#
    rmse = get_RMSE(f, σ)
    delta = σ .- f
    mae = maximum(abs.(delta))
    spline_min = minimum(σ)
    spline_max = maximum(σ)
    delta_min = minimum(delta)
    delta_max = maximum(delta)
    @printf "RMSE: %0.1e  MAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e delta_min:%0.1e delta_max:%0.1e\n" rmse mae spline_min spline_max delta_min delta_max
    open("c:/0/$model_id.txt","a") do io
        @printf io "model_id type_of_samples  n_of_samples  type_of_kernel   plot_grid_size  use_derivatives\n"
        @printf io "%2d      %2d             %4d             %1d               %3d               %s\n" model_id type_of_samples n_of_samples type_of_kernel  plot_grid_size use_derivatives
        @printf io "RMSE: %0.1e  MAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e   EPS:%0.1e   COND: %0.1e\n" rmse mae spline_min spline_max ε cond
        @printf io "c_time: %0.1e  e_time: %0.1e\n\n" c_time e_time
    end

    @printf "Creating plots..\n"

    PyPlot.clf()
    pygui(false)
    #PyPlot.view_init(30,-60)
    if model_id == 1
        PyPlot.view_init(30,30)
    end
    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=f, s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o)
    savefig("c:/0/m_t_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")
    PyPlot.clf()


    @printf "Plots created.\n"
#    return spline
    return Nothing
end
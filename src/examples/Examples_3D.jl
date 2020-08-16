using Printf
using PyPlot

function test_3D(model_id::Int,
                 use_derivatives::Bool = false,
                 type_of_samples::Int = 1,
                 n_of_samples::Int = 4,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 plot_grid_size::Int = 50
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
    u = Vector{Float64}(undef, n_1)
    grid = get_3D_plot_grid(plot_grid_size)
    m = size(grid, 2)
    f = Vector{Float64}(undef, m)
    if model_id == 1
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
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

        for i = 1:n_1
            u[i] = get_3D_model1(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_3D_model1(grid[:, i])
        end
    elseif model_id == 2
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            k += 1
            grad = get_2D_model2_grad(nodes[:, i])
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            du[k] = grad[2]
            es[1,k] = 0.0
            es[2,k] = 1.0
        end

        for i = 1:n_1
            u[i] = get_2D_model2(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model2(grid[:, i])
        end
    elseif model_id == 3
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            k += 1
            grad = get_2D_model3_grad(nodes[:, i])
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            k += 1
            d_nodes[1,k] = nodes[1,i]
            d_nodes[2,k] = nodes[2,i]
            du[k] = grad[2]
            es[1,k] = 0.0
            es[2,k] = 1.0
        end

        for i = 1:n_1
            u[i] = get_2D_model3(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model3(grid[:, i])
        end
    else
        error("Incorrect value of 'model_id'")
        return
    end

    return

    if use_derivatives
        @printf "nodes#: %d  d_nodes#: %d (total nodes: %d)\n" n_1 n_1 (n_1+n_1)
    else
        @printf "nodes#: %d\n" n_1
    end

    @printf "model_id type_of_samples  n_of_samples  type_of_kernel  regular_grid_size  use_derivatives\n"
    @printf "%2d      %2d             %4d             %1d               %3d               %s\n" model_id type_of_samples n_of_samples type_of_kernel regular_grid_size use_derivatives
#
    @printf "Creating spline..\n"
    ts = time_ns()
    if use_derivatives
        if rk.ε == 0
             epsilon = estimate_epsilon(nodes, d_nodes)
             @printf "Estimated EPSILON:%0.1e\n" epsilon
        end
        spline = interpolate(nodes, u, d_nodes, du, rk)
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
        @printf io "model_id type_of_samples  n_of_samples  type_of_kernel  regular_grid_size  use_derivatives\n"
        @printf io "%2d      %2d             %4d             %1d               %3d               %s\n" model_id type_of_samples n_of_samples type_of_kernel regular_grid_size use_derivatives
        @printf io "RMSE: %0.1e  MAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e   EPS:%0.1e   COND: %0.1e\n" rmse mae spline_min spline_max ε cond
        @printf io "c_time: %0.1e  e_time: %0.1e\n\n" c_time e_time
    end

    @printf "Creating plots..\n"
    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    ss = 40
    if n_of_samples == 3
        ss = 20
    end
    if n_of_samples == 4
        ss = 10
    end
    if n_of_samples == 5
        ss = 7
    end
    scatter(nodes, u, alpha=1.0, s=ss, c ="Red")
    if type_of_samples == 1
        PyPlot.title("Samples: $(samples_size[n_of_samples]) Halton nodes")
    elseif type_of_samples == 3
        PyPlot.title("Samples: $(samples_size[n_of_samples]) regular nodes")
    end
    savefig("c:/0/sca_$model_id,$type_of_samples,$n_of_samples,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    plot(grid, f, color="red")
    PyPlot.xlabel("x")
    PyPlot.ylabel("f")
    PyPlot.title("Function f(x)=1.0+x^2+log(abs(3.0*(1.0-x)+1.0))/3.3")
    savefig("c:/0/fun_$model_id,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    plot(grid, σ, color="blue", label="Spline")
    PyPlot.xlabel("x")
    PyPlot.ylabel("σ")
    if type_of_samples == 1
        PyPlot.title("Spline σ: $(samples_size[n_of_samples]) Halton nodes")
    elseif type_of_samples == 3
        PyPlot.title("Spline σ: $(samples_size[n_of_samples]) regular nodes")
    end
    scatter(nodes, u, alpha=1.0, s=ss, c ="Red", label="Nodes")
    legend()
    savefig("c:/0/spl_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,$use_derivatives,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    plot(grid, delta, color="green", label="Error")
    scatter(nodes, u, alpha=1.0, s=ss, c ="Red", label="Nodes")
    PyPlot.xlabel("x")
    PyPlot.ylabel("Error")
    if type_of_samples == 1
        PyPlot.title("Error (σ-f): $(samples_size[n_of_samples]) Halton nodes")
    elseif type_of_samples == 3
        PyPlot.title("Error (σ-f): $(samples_size[n_of_samples]) regular nodes")
    end
    legend()
    savefig("c:/0/err_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,$use_derivatives,_.png")

    PyPlot.clf()
    @printf "Plots created.\n"
#    return spline
    return Nothing
end

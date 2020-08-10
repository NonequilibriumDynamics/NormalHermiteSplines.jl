using Printf
using PyPlot

function test_1D(model_id::Int,
                 use_grad::Bool = true,
                 type_of_samples::Int = 3,
                 n_of_samples::Int = 1,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 regular_grid_size::Int = 1000
#                 ,do_parallel::Bool = false
                )
    if use_grad && type_of_kernel == 0
        error("Cannot use derivative data when type_of_kernel is `0` (`RK_H0` kernel)")
    end

    samples_size = [50, 100, 200, 400]
    if type_of_samples == 1
        nodes = get_1D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 2
        nodes = get_1D_eps_grid(samples_size[n_of_samples])
    elseif type_of_samples == 3
        nodes = get_1D_grid(samples_size[n_of_samples])
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

    n_1 = length(nodes)
    u = Vector{Float64}(undef, n_1)
    grid = get_1D_grid(regular_grid_size)
    m = length(grid)
    f = Vector{Float64}(undef, m)
    d_nodes = Vector{Float64}(undef, n_1)
    du = Vector{Float64}(undef, n_1)
    if model_id == 1
        k = 0
        for i = 1:n_1
            d_nodes[i] = nodes[i]
            du[i] = get_1D_model1_grad(nodes[i])
        end
        for i = 1:n_1
            u[i] = get_1D_model1(nodes[i])
        end
        for i = 1:m
            f[i] = get_1D_model1(grid[i])
        end
    else
        error("Incorrect value of 'model_id'")
        return
    end

    if use_grad
        @printf "nodes#: %d  d_nodes#: %d (total nodes: %d)\n" n_1 n_1 (n_1+n_1)
    else
        @printf "nodes#: %d\n" n_1
    end

    @printf "model_id type_of_samples  n_of_samples  type_of_kernel  regular_grid_size\n"
    @printf "%2d      %2d             %4d             %1d               %3d\n" model_id type_of_samples n_of_samples type_of_kernel regular_grid_size
#
    @printf "Creating spline..\n"
    ts = time_ns()
    if use_grad
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
        @printf io "model_id type_of_samples  n_of_samples  type_of_kernel  regular_grid_size\n"
        @printf io "%2d      %2d             %4d             %1d               %3d\n" model_id type_of_samples n_of_samples type_of_kernel regular_grid_size
        @printf io "RMSE: %0.1e  MAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e   EPS:%0.1e   COND: %0.1e\n" rmse mae spline_min spline_max ε cond
        @printf io "c_time: %0.1e  e_time: %0.1e\n\n" c_time e_time
    end

    @printf "Creating pictures..\n"

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    ss = 4
    if n_of_samples >= 3
        ss = 1
    end
    scatter(nodes, u, alpha=1.0, s=ss, c ="Red")
    savefig("c:/0/sca_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    plot(grid, f, color="red")
    savefig("c:/0/fun_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 2.5)
    plot(grid, σ, color="blue")
    savefig("c:/0/spl_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    pygui(false)
    PyPlot.xlim(1.3, 1.4)
    PyPlot.ylim(0.0, 1.1*mae)
    plot(grid, delta, color="green")
    savefig("c:/0/dlt_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    @printf "Pictures created.\n"
    return spline
end

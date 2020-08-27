using Printf
using PyPlot

function test_2D(model_id::Int,
                 use_grad::Bool,
                 type_of_samples::Int = 1,
                 n_of_samples::Int = 1,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 regular_grid_size::Int = 100
#                 ,do_parallel::Bool = false
                )
    if use_grad && type_of_kernel == 0
        error("Cannot use derivative data when type_of_kernel is `0` (`RK_H0` kernel)")
    end
    if type_of_samples == 1
        samples_size = [50, 100, 1000, 2500, 5000, 10000]
        nodes = get_2D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 2
        samples_size = [9, 15, 30, 50, 75, 100]
        nodes = get_2D_random_grid(samples_size[n_of_samples])
    elseif type_of_samples == 3
        samples_size = [9, 15, 30, 50, 75, 100]
        nodes = get_2D_grid(samples_size[n_of_samples])
    elseif type_of_samples == 4
        samples_size = [9, 15, 32, 49, 70, 99]
        nodes = get_2D_eps_grid(samples_size[n_of_samples])
    elseif type_of_samples == 5
        sets = [4; 6; 21; 34; 49; 70]
        nodes = get_2D_Lissajous_nodes(sets[n_of_samples])
    elseif type_of_samples == 6
        samples_size = [1, 11, 24, 250, 1250]
        nodes = get_2D_test1_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 32
        samples_size = [6, 9, 13, 19, 31, 49, 70, 99] # 49, 100, 196, 400, 1024, 2500, 5041, 10000
        # samples_size = [50, 100, 200, 400, 1000, 2500, 5000, 10000]
        nodes = get_2D_uniformrandom_grid(samples_size[n_of_samples])
        # nodes = get_2D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 33
        samples_size = [6, 9, 13, 19, 31, 49, 70, 99] # 49, 100, 196, 400, 1024, 2500, 5041, 10000
        # samples_size = [50, 100, 200, 400, 1000, 2500, 5000, 10000]
        nodes = get_2D_uniformrandom_grid(samples_size[n_of_samples])
        # nodes = get_2D_halton_nodes(samples_size[n_of_samples])
        bnodes = get_2D_border_nodes(19)
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
    grid = get_2D_grid(regular_grid_size)
    m = size(grid, 2)
    f = Vector{Float64}(undef, m)
    if model_id == 1
        # d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        # es = Matrix{Float64}(undef, 2, 2 * n_1)
        # du = Vector{Float64}(undef, 2 * n_1)
        # k = 0
        # for i = 1:n_1
        #     k += 1
        #     d_nodes[1,k] = nodes[1,i]
        #     d_nodes[2,k] = nodes[2,i]
        #     du[k] = 0.0
        #     es[1,k] = 1.0
        #     es[2,k] = 0.0
        #     k += 1
        #     d_nodes[1,k] = nodes[1,i]
        #     d_nodes[2,k] = nodes[2,i]
        #     du[k] = 0.0
        #     es[1,k] = 0.0
        #     es[2,k] = 1.0
        # end
        # d_nodes = d_nodes[:,1:k]
        # es = es[:,1:k]
        # du = du[1:k]

        d_nodes = get_2D_grid(100)
        d_n_1 = size(d_nodes, 2)
        d_nodes = [d_nodes d_nodes]
        es = Matrix{Float64}(undef, 2, 2 * d_n_1)
        du = Vector{Float64}(undef, 2 * d_n_1)
        k = 0
        for i = 1:d_n_1
            k += 1
            du[k] = 0.0
            es[1,k] = 1.0
            es[2,k] = 0.0
            k += 1
            du[k] = 0.0
            es[1,k] = 0.0
            es[2,k] = 1.0
        end

        for i = 1:n_1
            u[i] = get_2D_model1(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model1(grid[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:n_1
            u[i] = get_2D_model3(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model3(grid[:, i])
        end
    elseif model_id == 32
        if type_of_samples != 32
            error("Incorrect value of 'type_of_samples' for model #32.")
        end
        use_grad = false # never use the grad values for this test
        for i = 1:n_1
            u[i] = get_2D_model3(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model3(grid[:, i])
        end
    elseif model_id == 33
        if type_of_samples != 33
            error("Incorrect value of 'type_of_samples' for model #33.")
        end
        use_grad = true # always use the grad values for this test
        bn_1 = size(bnodes, 2)
        d_nodes = Matrix{Float64}(undef, 2, 2 * bn_1)
        es = Matrix{Float64}(undef, 2, 2 * bn_1)
        du = Vector{Float64}(undef, 2 * bn_1)
        k = 0
        for i = 1:bn_1
            k += 1
            grad = get_2D_model3_grad(bnodes[:, i])
            d_nodes[1,k] = bnodes[1,i]
            d_nodes[2,k] = bnodes[2,i]
            du[k] = grad[1]
            es[1,k] = 1.0
            es[2,k] = 0.0
            k += 1
            d_nodes[1,k] = bnodes[1,i]
            d_nodes[2,k] = bnodes[2,i]
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
    elseif model_id == 4
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            k += 1
            grad = get_2D_model4_grad(nodes[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:n_1
            u[i] = get_2D_model4(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model4(grid[:, i])
        end
    elseif model_id == 5
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            k += 1
            grad = get_2D_model5_grad(nodes[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:n_1
            u[i] = get_2D_model5(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model5(grid[:, i])
        end
    elseif model_id == 6
        grid = get_2D_grid2(regular_grid_size)
        m = size(grid, 2)
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            nodes[1, i] = nodes[1, i] * 2.0 - 1.0
            nodes[2, i] = nodes[2, i] * 2.0 - 1.0
            u[i] = get_2D_model6(nodes[:, i])
            k += 1
            grad = get_2D_model6_grad(nodes[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:m
            f[i] = get_2D_model6(grid[:, i])
        end
    elseif model_id == 10
        for i = 1:n_1
            u[i] = get_2D_model10(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model10(grid[:, i])
        end
    elseif model_id == 11
        for i = 1:n_1
            u[i] = get_2D_model11(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model11(grid[:, i])
        end
    elseif model_id == 12
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        grid = get_2D_rect_grid(regular_grid_size)
        m = size(grid, 2)
        k = 0
        for i = 1:n_1
            nodes[1, i] = nodes[1, i] * 8.0 - 4.0
            nodes[2, i] = nodes[2, i] * 4.0 - 2.0
            u[i] = get_2D_model12(nodes[:, i])
            k += 1
            grad = get_2D_model12_grad(nodes[:, i])
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

        f = Vector{Float64}(undef, m)
        for i = 1:m
            f[i] = get_2D_model12(grid[:, i])
        end
    elseif model_id == 13
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            k += 1
            grad = get_2D_model13_grad(nodes[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:n_1
            u[i] = get_2D_model13(nodes[:, i])
        end
        for i = 1:m
            f[i] = get_2D_model13(grid[:, i])
        end
    elseif model_id == 14
        #  (x,y) \in [-2;2]x[-1;3]
        grid = get_2D_grid3(regular_grid_size)
        m = size(grid, 2)
        d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)
        es = Matrix{Float64}(undef, 2, 2 * n_1)
        du = Vector{Float64}(undef, 2 * n_1)
        k = 0
        for i = 1:n_1
            nodes[1, i] = nodes[1, i] * 4.0 - 2.0
            nodes[2, i] = nodes[2, i] * 4.0 - 1.0
            u[i] = get_2D_model14(nodes[:, i])
            k += 1
            grad = get_2D_model14_grad(nodes[:, i])
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
        d_nodes = d_nodes[:,1:k]
        es = es[:,1:k]
        du = du[1:k]

        for i = 1:m
            f[i] = get_2D_model14(grid[:, i])
        end
    else
        error("Incorrect value of 'model_id'")
        return
    end

    if use_grad
        @printf "nodes#: %d  d_nodes#: %d (total nodes: %d)\n" n_1 k (n_1+k)
    else
        @printf "nodes#: %d\n" n_1
    end

    @printf "model_id type_of_samples  n_of_samples  type_of_kernel  regular_grid_size\n"
    @printf "%2d      %2d             %4d             %1d               %3d\n" model_id type_of_samples n_of_samples type_of_kernel regular_grid_size

    @printf "Creating spline..\n"
    ts = time_ns()
    if use_grad
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

    rmse = get_RMSE(f, σ)
    delta = f .- σ
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

    @printf "Creating plots..\n"
    gx = grid[1,:]
    gy = grid[2,:]
    x = unique(grid[1,:])
    y = unique(grid[2,:])
    gσ = reshape(σ, length(y), length(x))
    gf = reshape(f, length(y), length(x))
    ss = (6 - n_of_samples) > 1 ? (6 - n_of_samples)/n_of_samples : 2.0/n_of_samples

    if model_id == 1
        lvls = [0.0;0.05;0.1;0.3;0.5;0.7;0.8;0.9;0.95;1.0]
        lvls2 = [-0.001;0.0;0.05;0.1;0.3;0.5;0.7;0.8;0.9;0.95;1.0;1.001;]
    end
    if model_id == 2
        lvls = [-0.1;0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.1]
        lvls2 = lvls
    end
    if model_id == 3 || model_id == 32 || model_id == 33
        lvls=[-0.1;-0.05;0.0;0.05;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1.0;1.05;1.1]
        lvls2 = lvls
    end
    if model_id == 4
        lvls=[-0.6;-0.55;-0.5;-0.4;-0.3;-0.2;-0.1;0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.05;1.1]
        lvls2 = lvls
    end
    if model_id == 5
        lvls=[-1.2;-1.15;-1.05;-1.0;-0.95;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;0.95;1.0;1.05;1.1;1.15;1.2]
        lvls2 = lvls
    end
    if model_id == 6
        lvls=[-1.0;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;1.0]
        lvls2 = lvls
#        lvls=[-1.1;-1.0;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;1.0;1.1]
    end

    if model_id == 10
        lvls=[0.8; 0.85;0.90;0.95;1.0;1.05;1.1;1.15;1.2]
        lvls2 = lvls
    end
    if model_id == 11
        lvls=[0.0;0.1;0.2;0.3;0.39]
        lvls2 = lvls
    end
    if model_id == 12
        lvls=[-0.1;-0.05;0.0;0.05;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1.0;1.05;1.1]
        lvls2 = lvls
    end
    if model_id == 13
        lvls=[-0.1;0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.1;1.2;1.3]
        lvls2 = lvls
    end
    if model_id == 14
        lvls=[-0.1;0.0;0.01;0.02;0.03;0.04;0.05;0.06;0.07;0.08;0.09;0.1;0.12;0.14;0.16;0.18;0.2;0.22;0.24;0.26;0.28;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.05;1.1]
        lvls2 = lvls
    end

    PyPlot.clf()
    pygui(false)
    o = contourf(x, y, gσ, levels=lvls2, cmap=ColorMap("gnuplot"))
    axis("equal")
    # if n_of_samples <= 2
    #     scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
    # end
    if model_id == 6
        PyPlot.xlim(-1.0, 1.0)
        PyPlot.ylim(-1.0, 1.0)
    end
    if model_id == 12
        PyPlot.xlim(-4.0, 4.0)
        PyPlot.ylim(-2.0, 2.0)
    end
    if model_id == 14
        PyPlot.xlim(-2.0, 2.0)
        PyPlot.ylim(-1.0, 3.0)
    end

    colorbar(o)
    savefig("c:/0/s_cf_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")

    PyPlot.clf()
    pygui(false)
    scatter(nodes[1,:], nodes[2,:], s=ss)
    if model_id == 33
        scatter(bnodes[1,:], bnodes[2,:], s=(2*ss), c="red")
    end
    gca().set_aspect("equal")
    savefig("c:/0/m_grid_$type_of_samples,$n_of_samples.png", dpi=150, bbox_inches="tight")

    PyPlot.clf()
    pygui(false)
    o = contourf(x, y, gf, levels=lvls, cmap=ColorMap("gnuplot"))
#    o = contourf(x, y, gf, levels=lvls, cmap=ColorMap("gnuplot"))
    axis("equal")
    # if n_of_samples <= 2
    #     scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
    # end
    if model_id == 6
        PyPlot.xlim(-1.0, 1.0)
        PyPlot.ylim(-1.0, 1.0)
    end
    if model_id == 12
        PyPlot.xlim(-4.0, 4.0)
        PyPlot.ylim(-2.0, 2.0)
    end
    if model_id == 14
        PyPlot.xlim(-2.0, 2.0)
        PyPlot.ylim(-1.0, 3.0)
    end
    colorbar(o)
    savefig("c:/0/m_cf_$model_id.png", dpi=150, bbox_inches="tight")

    PyPlot.clf()
    pygui(false)
    o = contour(x, y, gf, levels=lvls, cmap=ColorMap("gnuplot"))
    #    o = contourf(x, y, gf, levels=lvls, cmap=ColorMap("gnuplot"))
    axis("equal")
    # if n_of_samples <= 2
    #     scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
    # end
    if model_id == 6
        PyPlot.xlim(-1.0, 1.0)
        PyPlot.ylim(-1.0, 1.0)
    end
    if model_id == 12
        PyPlot.xlim(-4.0, 4.0)
        PyPlot.ylim(-2.0, 2.0)
    end
    if model_id == 14
        PyPlot.xlim(-2.0, 2.0)
        PyPlot.ylim(-1.0, 3.0)
    end
    colorbar(o)
    savefig("c:/0/m_c_$model_id.png", dpi=150, bbox_inches="tight")

    PyPlot.clf()
    pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    #PyPlot.view_init(30,-60) #default
    if model_id == 13
        PyPlot.view_init(30,30)
    end
    if model_id == 14
        PyPlot.view_init(30,-120)
    end
    o = scatter3D(grid[1,:],grid[2,:], f, c=f,  s=1, cmap=ColorMap("gnuplot"))
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/0/m_t_$model_id.png", dpi=150, bbox_inches="tight")

#     PyPlot.clf()
#     pygui(false)
#     o = contour(x, y, gσ, levels=lvls)
#     axis("equal")
#     if n_of_samples <= 2
#         scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
#     end
#     if model_id == 6
#         PyPlot.xlim(-1.0, 1.0)
#         PyPlot.ylim(-1.0, 1.0)
#     end
#     if model_id == 12
# #        axis("equal")
#         PyPlot.xlim(-4.0, 4.0)
#         PyPlot.ylim(-2.0, 2.0)
#     end
#     colorbar(o)
#     savefig("c:/0/s_c_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")
#     PyPlot.clf()

    # pygui(false)
    # # if model_id == 3
    # #     PyPlot.view_init(20,30)
    # # end
    # #o = surf(gx, gy, σ, cmap=ColorMap("viridis"), alpha=0.75)
    # o = surf(gx, gy, σ, cmap=ColorMap("viridis"), linewidth=0, antialiased=false, alpha=1.0)
    # tick_params(axis="both", which="major", labelsize=6)
    # tick_params(axis="both", which="minor", labelsize=6)
    # colorbar(o)
    #
    # #scatter(nodes[1,:], nodes[2,:], c="red", s= ss, zdir="z")
    # #colorbar(o)
    # savefig("c:/0/s_s_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")
    # PyPlot.clf()
    #
    # pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    # if model_id == 13
    #     PyPlot.view_init(30,30)
    # end
    # #o = surf(gx, gy, σ, cmap=ColorMap("viridis"), alpha=0.75)
    # o = surf(gx, gy, f, cmap=ColorMap("viridis"), linewidth=0, antialiased=false, alpha=1.0)
    # tick_params(axis="both", which="major", labelsize=6)
    # tick_params(axis="both", which="minor", labelsize=6)
    # # colorbar(o)
    # savefig("c:/0/m_s_$model_id.png")
    #

    PyPlot.clf()
    pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    #PyPlot.view_init(30,-60)
    if model_id == 13
        PyPlot.view_init(30,30)
    end
    if model_id == 14
        PyPlot.view_init(30,-120)
    end
    o = scatter3D(grid[1,:],grid[2,:], σ, c=σ, s=1, cmap=ColorMap("gnuplot"))
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/0/s_t_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")
    PyPlot.clf()
######
    gd = reshape(delta, length(y), length(x))
    PyPlot.clf()
    pygui(false)
    o = contourf(x, y, gd, cmap=ColorMap("jet"))
    axis("equal")
    # if n_of_samples <= 2
    #     scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
    # end
    if model_id == 6
        PyPlot.xlim(-1.0, 1.0)
        PyPlot.ylim(-1.0, 1.0)
    end
    if model_id == 12
    #        axis("equal")
        PyPlot.xlim(-4.0, 4.0)
        PyPlot.ylim(-2.0, 2.0)
    end
    if model_id == 14
        PyPlot.xlim(-2.0, 2.0)
        PyPlot.ylim(-1.0, 3.0)
    end
    colorbar(o)
    savefig("c:/0/delta_cf_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")

    #PyPlot.clf()
    # pygui(false)
    # o = contour(x, y, gd)
    # axis("equal")
    # if n_of_samples <= 2
    #     scatter(nodes[1,:], nodes[2,:], c="red", s= ss)
    # end
    # if model_id == 6
    #     PyPlot.xlim(-1.0, 1.0)
    #     PyPlot.ylim(-1.0, 1.0)
    # end
    # if model_id == 12
    # #        axis("equal")
    #     PyPlot.xlim(-4.0, 4.0)
    #     PyPlot.ylim(-2.0, 2.0)
    # end
    # colorbar(o)
    # savefig("c:/0/delta_c_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    #o = surf(gx, gy, σ, cmap=ColorMap("viridis"), alpha=0.75)
    if model_id == 13
        PyPlot.view_init(30,30)
    end
    if model_id == 14
        PyPlot.view_init(30,-120)
    end
    o = surf(gx, gy, delta, cmap=ColorMap("jet"), linewidth=0, antialiased=false, alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/0/delta_s_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")

    #PyPlot.clf()
    # pygui(false)
    # # if model_id == 3
    # #     PyPlot.view_init(20,30)
    # # end
    # o = scatter3D(grid[1,:],grid[2,:], delta, c=delta,  s=1)
    # tick_params(axis="both", which="major", labelsize=6)
    # tick_params(axis="both", which="minor", labelsize=6)
    # colorbar(o)
    # savefig("c:/0/delta_t_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    @printf "Plots created.\n"
    return nothing
#    return spline
end

function readme_1()
    nodes = get_2D_halton_nodes(1000)             # generates Halton data set in [0, 1] x [0, 1]
    n_1 = size(nodes, 2)
    u = Vector{Float64}(undef, n_1)               # function values
    grid = get_2D_grid2(100)                       # regular grid covering [-1, 1] x [-1, 1]
    d_nodes = Matrix{Float64}(undef, 2, 2 * n_1)  # directional derivative nodes
    es = Matrix{Float64}(undef, 2, 2 * n_1)       # derivative directions
    du = Vector{Float64}(undef, 2 * n_1)          # directional derivative values

    k = 0
    for i = 1:n_1
        nodes[1, i] = nodes[1, i] * 2.0 - 1.0     # transforming Halton nodes to [-1, 1] x [-1, 1]
        nodes[2, i] = nodes[2, i] * 2.0 - 1.0
        u[i] = sin(4.0 * sqrt(nodes[1, i]^2 + nodes[2, i]^2))
        k += 1
        grad = [0.0; 0.0]
        d = sqrt(nodes[1, i]^2 + nodes[2, i]^2)
        if d > 0.0
            grad[1] = 4.0 * nodes[1, i] * cos(4.0 * d) / d
            grad[2] = 4.0 * nodes[2, i] * cos(4.0 * d) / d
        end
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

    # Hermite spline must be constructed with ```RK_H1``` or ```RK_H2``` kernel,
    # the 'scaling parameter' ```ε``` here is estimated in the '''interpolate''' procedure
    rk = RK_H1()
    #
    spline = interpolate(nodes, u, d_nodes, es, du, rk)
    σ = evaluate(spline, grid)
    return σ
end

function readme_2()
    nodes = get_2D_Lissajous_nodes(49)            # generates 4999 LS_2^{(49,50)} Lissajous nodes
    n_1 = size(nodes, 2)
    u = Vector{Float64}(undef, n_1)               # function values
    grid = get_2D_grid(100)                       # uniform Cartesian grid of size 101x101 in [0, 1] x [0, 1]
    for i = 1:n_1
        val = 0.0
        x = nodes[1,i]
        y = nodes[2,i]
        d = (x - 0.75)^2 + (y - 0.75)^2
        r2 = 0.04
        if d <= r2
            val = 1.0
        end
        if x <= 0.25 && x >= 0.05 && y <= 0.25 && y >= 0.05
            val = 1.0
        end
        u[i] = val
    end

    # Here spline is being constructed with ```RK_H0``` kernel,
    # the 'scaling parameter' ```ε``` is defined explicitly.
    rk = RK_H0(1.0)
    #
    spline = interpolate(nodes, u, rk)
    σ = evaluate(spline, grid)
    return σ
end

function big_sur()
    data = [
            1 92.2 2.0 10.39
            2 88.8 2.3 10.23
            3 87.8 3.3 10.34
            4 71.0 6.0 11.65
            5 67.3 6.6 12.25
            6 62.2 5.3 12.84
            7 48.8 2.3 13.58
            8 39.6 5.4 14.11
            9 31.0 4.9 14.21
            10 21.0 6.1 14.19
            11 12.8 8.0 14.16
            12 87.2 34.8 10.20
            13 86.7 34.1 10.45
            14 85.6 33.8 10.97
            15 84.9 33.9 10.81
            16 84.0 33.9 10.67
            17 77.3 35.7 11.53
            18 73.2 34.7 12.38
            19 68.0 35.1 12.50
            20 63.3 33.8 12.85
            21 52.3 35.5 13.09
            22 44.5 34.0 13.68
            23 32.0 33.9 14.01
            24 20.9 35.9 13.65
            25 12.4 35.5 14.08
            26 2.3 32.9 14.25
            27 86.7 64.9 11.41
            28 85.9 64.8 11.02
            29 84.8 64.9 11.09
            30 83.5 64.4 10.62
            31 77.9 64.9 11.09
            32 74.1 63.7 12.13
            33 67.4 65.0 12.78
            34 62.9 63.9 12.73
            35 51.8 64.9 11.66
            36 41.3 65.6 12.63
            37 29.7 62.0 13.18
            38 22.3 65.2 13.42
            39 13.0 65.9 13.22
            40 2.7 64.8 13.01
            41 3.0 94.0 13.09
            42 11.5 96.8 12.84
            43 22.3 95.9 12.99
            44 32.0 94.0 13.30
            45 42.1 95.8 12.94
            46 53.1 94.2 12.25
            47 63.1 94.2 10.95
            48 66.0 94.2 10.66
            49 66.0 93.0 10.58
            50 67.3 95.2 10.73
            51 68.8 94.4 10.48
            52 69.9 94.0 10.38
            53 71.8 94.5 10.45
            54 73.8 94.9 10.27
            55 75.2 94.5 10.22
            56 83.0 124.5 13.41
            57 72.0 126.0 13.04
            58 62.5 123.3 13.38
            59 51.9 124.8 12.89
            60 32.5 126.5 13.22
            61 27.1 126.0 13.03
            62 20.2 125.3 11.86
            63 7.9 127.6 11.99
            64 2.0 123.0 12.06
           ]
    return data
end

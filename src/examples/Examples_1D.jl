using Printf
using PyPlot

function test_1D(model_id::Int,
                 type_of_samples::Int = 2,
                 n_of_samples::Int = 1,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 regular_grid_size::Int = 40,
                 use_grad::Bool = true
#                 ,do_parallel::Bool = false
                )
    if use_grad && type_of_kernel == 0
        error("Cannot use derivative data when type_of_kernel is `0` (`RK_H0` kernel)")
    end
    if type_of_samples == 1
        samples_size = [1, 11, 24, 250, 1250]
        nodes = get_2D_test1_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 2
        samples_size = [50, 100, 1000, 2500, 5000, 10000]
        nodes = get_2D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 3
        sets = [4; 6; 21; 34; 49; 70]
        nodes = get_2D_Lissajous_nodes(sets[n_of_samples])
    elseif type_of_samples == 4
        samples_size = [9, 15, 32, 49, 70, 99]
        nodes = get_2D_eps_grid(samples_size[n_of_samples])
    elseif type_of_samples == 5
        samples_size = [9, 15, 30, 50, 75, 100]
        nodes = get_2D_grid(samples_size[n_of_samples])
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

    @printf "Creating pictures..\n"
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
    if model_id == 3
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
    colorbar(o)
    savefig("c:/0/s_cf_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

    PyPlot.clf()
    pygui(false)
    scatter(nodes[1,:], nodes[2,:], s= ss)
    gca().set_aspect("equal")
    savefig("c:/0/m_grid_$model_id,$type_of_samples,$n_of_samples.png")

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
    colorbar(o)
    savefig("c:/0/m_cf_$model_id.png")

    PyPlot.clf()
    pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    #PyPlot.view_init(30,-60)
    o = scatter3D(grid[1,:],grid[2,:], f, c=σ,  s=1, cmap=ColorMap("gnuplot"))
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o)
    savefig("c:/0/m_t_$model_id.png")

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
    # # if model_id == 3
    # #     PyPlot.view_init(20,30)
    # # end
    # #o = surf(gx, gy, σ, cmap=ColorMap("viridis"), alpha=0.75)
    # o = surf(gx, gy, f, cmap=ColorMap("viridis"), linewidth=0, antialiased=false, alpha=1.0)
    # tick_params(axis="both", which="major", labelsize=6)
    # tick_params(axis="both", which="minor", labelsize=6)
    # colorbar(o)
    # savefig("c:/0/m_s_$model_id.png")
    #

    PyPlot.clf()
    pygui(false)
    # if model_id == 3
    #     PyPlot.view_init(20,30)
    # end
    #PyPlot.view_init(30,-60)
    o = scatter3D(grid[1,:],grid[2,:], σ, c=σ, s=1, cmap=ColorMap("gnuplot"))
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o)
    savefig("c:/0/s_t_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")
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
    colorbar(o)
    savefig("c:/0/delta_cf_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

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
    o = surf(gx, gy, delta, cmap=ColorMap("jet"), linewidth=0, antialiased=false, alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o)
    savefig("c:/0/delta_s_$model_id,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

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
    @printf "Pictures created.\n"
    return spline
end

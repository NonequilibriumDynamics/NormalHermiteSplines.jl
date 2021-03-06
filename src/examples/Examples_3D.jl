using Printf
using PyPlot
using PyCall

function test_3D(model_id::Int,
                 use_grad::Bool = false,
                 type_of_samples::Int = 1,
                 n_of_samples::Int = 4,
                 type_of_kernel::Int = 0,
                 eps::Float64 = 0.0,
                 plot_grid_type::Int = 1,
                 plot_grid_size = 25, # plot_grid_type = 1: 50, 25, 20, 16, 30
                                      # plot_grid_type = 2: 50, 40, 75
                 #,do_parallel::Bool = false
                )
    if use_grad && type_of_kernel == 0
        error("Cannot use derivative data when type_of_kernel is `0` (`RK_H0` kernel)")
    end

    if type_of_samples == 1
        samples_size = [1, 100, 500, 1000, 2000, 4000, 8000, 16000]
        nodes = get_3D_halton_nodes(samples_size[n_of_samples])
    elseif type_of_samples == 2
        samples_size = [1, 4, 7, 9, 12, 15, 19, 24]
        nodes = get_3D_random_grid(samples_size[n_of_samples])
    elseif type_of_samples == 3
        samples_size = [1, 4, 7, 9, 12, 15, 19, 24]
        nodes = get_3D_grid(samples_size[n_of_samples]) #1(8), 4(125), 7(512), 9(1000), 12(2197), 15(4096)(d), 19(8000), 24(15625)
    elseif type_of_samples == 4
        samples_size = [1, 4, 7, 9, 12, 15, 19, 24]
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

    if plot_grid_type == 1
        grid = get_3D_grid(plot_grid_size)
    else
        grid = get_3D_plot_grid(plot_grid_size)
    end
    m = size(grid, 2)
    n_1 = size(nodes, 2)
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
            grad = get_3D_model3_grad(nodes[:, i])
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
    elseif model_id == 4
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model4(nodes[:, i])
            k += 1
            grad = get_3D_model4_grad(nodes[:, i])
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
            f[i] = get_3D_model4(grid[:, i])
        end
    elseif model_id == 5
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model5(nodes[:, i])
            k += 1
            grad = get_3D_model5_grad(nodes[:, i])
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
            f[i] = get_3D_model5(grid[:, i])
        end
    elseif model_id == 6
        u = Vector{Float64}(undef, n_1)
        f = Vector{Float64}(undef, m)
        d_nodes = Matrix{Float64}(undef, 3, 3*n_1)
        es = Matrix{Float64}(undef, 3, 3*n_1)
        du = Vector{Float64}(undef, 3*n_1)
        k = 0
        for i = 1:n_1
            u[i] = get_3D_model6(nodes[:, i])
            k += 1
            grad = get_3D_model6_grad(nodes[:, i])
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
            f[i] = get_3D_model6(grid[:, i])
        end
    else
        error("Incorrect value of 'model_id'")
        return
    end

    if use_grad
        @printf "nodes#: %d  d_nodes#: %d (total nodes: %d)  grid %d\n" n_1 k (n_1+k) m
    else
        @printf "nodes#: %d  grid %d\n" n_1 m
    end

    @printf "model_id type_of_samples  n_of_samples  type_of_kernel   plot_grid_size  use_grad\n"
    @printf "%2d      %2d             %4d             %1d                 %5d            %s\n" model_id type_of_samples n_of_samples type_of_kernel m use_grad
#
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

    iq = estimate_accuracy(spline)
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
    rrmse = get_RRMSE(f, σ)
    rmae = get_RMAE(f, σ)
    spline_min = minimum(σ)
    spline_max = maximum(σ)
    delta_min = minimum(delta)
    delta_max = maximum(delta)
    @printf "RMSE: %0.1e  MAE:%0.1e RRMSE: %0.1e RMAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e delta_min:%0.1e delta_max:%0.1e\n" rmse mae rrmse rmae spline_min spline_max delta_min delta_max
    open("c:/000/$model_id.txt","a") do io
        @printf io "model_id type_of_samples  n_of_samples  type_of_kernel   plot_grid_size  use_grad\n"
        @printf io "%2d      %2d             %4d             %1d                 %5d            %s\n" model_id type_of_samples n_of_samples type_of_kernel m use_grad
        @printf io "RMSE: %0.1e  MAE:%0.1e RRMSE: %0.1e RMAE:%0.1e  SPLINE_MIN:%0.1e  SPLINE_MAX:%0.1e   EPS:%0.1e   COND: %0.1e\n" rmse mae rmse rmae spline_min spline_max ε cond
        @printf io "c_time: %0.1e  e_time: %0.1e\n\n" c_time e_time
    end

    @printf "Creating plots..\n"

    PyPlot.clf()
    pygui(false)
    if model_id == 1
        PyPlot.view_init(30,30)
    else
        PyPlot.view_init(20,30)
    end
    o = scatter3D(nodes[1,:], nodes[2,:], nodes[3,:], c=u, s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
#    o = scatter3D(nodes[1,:], nodes[2,:], nodes[3,:], c=(u./maximum(u)), s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/000/m_nodes_$model_id,$type_of_samples,$n_of_samples.png", dpi=150, bbox_inches="tight")

    PyPlot.clf()
    pygui(false)
    #PyPlot.view_init(30,-60)
    if model_id == 1
        PyPlot.view_init(30,30)
    else
        PyPlot.view_init(20,30)
    end
    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=f, s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
#    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=(f./maximum(f)), s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/000/m_grid_$model_id,$plot_grid_type,$plot_grid_size,_.png", dpi=150, bbox_inches="tight") # dpi=300 for bettr quality

    PyPlot.clf()
    PyPlot.clf()
    pygui(false)
    #PyPlot.view_init(30,-60)
    if model_id == 1
        PyPlot.view_init(30,30)
    else
        PyPlot.view_init(20,30)
    end
    if model_id == 4
        hmax = 2.0*maximum(σ)
        o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=(σ./hmax) , s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    else
        o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=σ, s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    end
#    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=(σ./maximum(σ)), s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
    savefig("c:/000/s_grid_$model_id,$use_grad,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,$plot_grid_type,$plot_grid_size,_.png", dpi=150, bbox_inches="tight")
    PyPlot.clf()

    PyPlot.clf()
    pygui(false)
    #PyPlot.view_init(30,-60)
    if model_id == 1
        PyPlot.view_init(30,30)
    else
        PyPlot.view_init(20,30)
    end
    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=delta, s=1, cmap=ColorMap("jet"), alpha=1.0)
#    o = scatter3D(grid[1,:],grid[2,:], grid[3,:], c=(delta./maximum(delta)) , s=1, cmap=ColorMap("gnuplot"), alpha=1.0)
    tick_params(axis="both", which="major", labelsize=6)
    tick_params(axis="both", which="minor", labelsize=6)
    colorbar(o, shrink=0.75)
#    gca().set_axis_off() # hide grid

    savefig("c:/000/delta_t_$model_id,$use_grad,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,$plot_grid_type,$plot_grid_size,_.png", dpi=150, bbox_inches="tight")
    PyPlot.clf()

    if plot_grid_type == 2 && model_id != 1
       zk = [0.0; 0.25; 0.5; 0.75; 1.0]
       for k = 1:5
           ik = grid[3,:] .== zk[k]
           gridk = grid[:,ik]
           fk = f[ik]
           σk = σ[ik]

           gx = gridk[1,:]
           gy = gridk[2,:]
           x = unique(gridk[1,:])
           y = unique(gridk[2,:])
           gf = reshape(fk, length(y), length(x))
           gσ = reshape(σk, length(y), length(x))

       #     if model_id == 4
       #         lvls=[-0.6;-0.55;-0.5;-0.4;-0.3;-0.2;-0.1;0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.05;1.1]
       #         lvls2 = lvls
       #     end
       #     if model_id == 5
       #         lvls=[-1.2;-1.15;-1.05;-1.0;-0.95;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;0.95;1.0;1.05;1.1;1.15;1.2]
       #         lvls2 = lvls
       #     end
       #     if model_id == 6
       #         lvls=[-1.0;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;1.0]
       #         lvls2 = lvls
       # #        lvls=[-1.1;-1.0;-0.9;-0.7;-0.5;-0.3;-0.1;0.0;0.1;0.3;0.5;0.7;0.9;1.0;1.1]
       #     end
           PyPlot.clf()
           pygui(false)
           o = contourf(x, y, gσ, cmap=ColorMap("gnuplot"))
#           o = contourf(x, y, gσ, levels=lvls2, cmap=ColorMap("gnuplot"))
           axis("equal")

           # if model_id == 6
           #     PyPlot.xlim(-1.0, 1.0)
           #     PyPlot.ylim(-1.0, 1.0)
           # end
           # if model_id == 12
           #     PyPlot.xlim(-4.0, 4.0)
           #     PyPlot.ylim(-2.0, 2.0)
           # end
           colorbar(o)
           savefig("c:/000/sc_$model_id,_$k,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")

           PyPlot.clf()
           pygui(false)
           o = contourf(x, y, gf, cmap=ColorMap("gnuplot"))
#           o = contourf(x, y, gf, levels=lvls2, cmap=ColorMap("gnuplot"))
           axis("equal")
           # if model_id == 6
           #     PyPlot.xlim(-1.0, 1.0)
           #     PyPlot.ylim(-1.0, 1.0)
           # end
           # if model_id == 12
           #     PyPlot.xlim(-4.0, 4.0)
           #     PyPlot.ylim(-2.0, 2.0)
           # end
           colorbar(o)
           savefig("c:/000/fc_$model_id,_$k,$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png", dpi=150, bbox_inches="tight")

           # PyPlot.clf()
           # pygui(false)
           # # if model_id == 3
           # #     PyPlot.view_init(20,30)
           # # end
           # #PyPlot.view_init(30,-60) #default
           # # if model_id == 13
           # #     PyPlot.view_init(30,30)
           # # end
           # o = scatter3D(gridk[1,:],gridk[2,:], fk, c=fk,  s=1, cmap=ColorMap("gnuplot"))
           # tick_params(axis="both", which="major", labelsize=6)
           # tick_params(axis="both", which="minor", labelsize=6)
           # colorbar(o)
           # savefig("c:/000/m_t_$model_id,_$k.png", dpi=150, bbox_inches="tight")

           # PyPlot.clf()
           # pygui(false)
           # # if model_id == 3
           # #     PyPlot.view_init(20,30)
           # # end
           # #PyPlot.view_init(30,-60) #default
           # # if model_id == 13
           # #     PyPlot.view_init(30,30)
           # # end
           # o = scatter3D(gridk[1,:],gridk[2,:], σk, c=σk,  s=1, cmap=ColorMap("gnuplot"))
           # tick_params(axis="both", which="major", labelsize=6)
           # tick_params(axis="both", which="minor", labelsize=6)
           # colorbar(o)
           # savefig("c:/000/s_t_$model_id,_$k,_$type_of_samples,$n_of_samples,$type_of_kernel,_$eps,_.png")

########## Animation

           # PyPlot.clf()
           # pygui(false)
           # fig = figure(figsize=(5,5))

           if k == 5
               global gf5 = gf
               global gσ5 = gσ
#               imshow(gf5, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
           elseif k == 4
               global gf4 = gf
               global gσ4 = gσ
#               imshow(gf4, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
           elseif k == 3
               global gf3 = gf
               global gσ3 = gσ
#               imshow(gf3, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
           elseif k == 2
               global gf2 = gf
               global gσ2 = gσ
#               imshow(gf2, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
           elseif k == 1
               global gf1 = gf
               global gσ1 = gσ
#               imshow(gf1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
           end
           # savefig("c:/000/imshow_$model_id,_$k.png")
           # PyPlot.clf()
       end # for k = 1:5

       PyPlot.clf()
       pygui(false)
       PyPlot.rc("figure", max_open_warning = 0)
       anim = pyimport("matplotlib.animation")
       fig = figure(figsize=(2.75,2.75))
#       fig = figure(figsize=(3,3))
       withfig(fig) do
           anim = anim.FuncAnimation(fig, make_frame_gf, init_func=init_gf, frames=22, interval=500)
           anim.save("c:/000/model.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
       end

       PyPlot.clf()
       pygui(false)
       PyPlot.rc("figure", max_open_warning = 0)
       anim = pyimport("matplotlib.animation")
       fig = figure(figsize=(2.75,2.75))
#       fig = figure(figsize=(3,3))
       withfig(fig) do
           anim = anim.FuncAnimation(fig, make_frame_gs, init_func=init_gs, frames=22, interval=500, repeat=true)
           anim.save("c:/000/spline.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
       end
       PyPlot.clf()
   end # if plot_grid_type == 2

    @printf "Plots created.\n"
    # return gr
    return Nothing
#    return spline
end

function make_frame_gf(k)
    if k >= 17 && k <= 20
        imshow(gf5, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 13 && k <= 16
        imshow(gf4, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 9 && k <= 12
        imshow(gf3, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 5 && k <= 8
        imshow(gf2, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 1 && k <=4
        imshow(gf1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 21 && k <=22
        imshow(gf1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    end
end

function init_gf()
    imshow(gf1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
end

function init_gs()
    imshow(gσ1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
end

function make_frame_gs(k)
    if k >= 17 && k <= 20
        imshow(gσ5, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 13 && k <= 16
        imshow(gσ4, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 9 && k <= 12
        imshow(gσ3, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 5 && k <= 8
        imshow(gσ2, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 1 && k <=4
        imshow(gσ1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    elseif k >= 21 && k <=22
        imshow(gσ1, interpolation="none", extent=[0, 1, 0, 1], origin="lower", cmap=ColorMap("gnuplot"))
    end
end

global gf1
global gf2
global gf3
global gf4
global gf5
global gσ5
global gσ4
global gσ3
global gσ2
global gσ1

function readme_3()
    nodes = get_3D_random_grid(9)       # generates 1000 non-uniform random grid nodes
    n_1 = size(nodes, 2)
    u = Vector{Float64}(undef, n_1)     # function values
    grid = get_3D_grid(50)              # uniform Cartesian grid of size 51x51x51 in [0, 1] x [0, 1] x [0, 1]
    for i = 1:n_1
        x = nodes[1,i]
        y = nodes[2,i]
        z = nodes[3,i]
        u[i] = cos(π*x)*cos(y - 0.5)*sin(π*(z - 0.5))
    end

    # Here spline is being constructed with ```RK_H2``` kernel,
    # the 'scaling parameter' ```ε``` is defined explicitly.
    rk = RK_H2(5.0)
    #
    spline = interpolate(nodes, u, rk)
    σ = evaluate(spline, grid)
    return σ
end

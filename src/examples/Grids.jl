using Random

########### 1D grids

function get_1D_grid(m::Int)
    return collect(range(1.3, 1.4; step = 0.1/m))
end

function get_1D_eps_grid(m::Int)
    x = collect(range(1.3, 1.4; step = 0.1/m))
    m1 = m + 1
    Random.seed!(123)
    eps = 0.1/m
    for i = 1:m1
        x[i] = x[i] + eps * (rand() - 0.5)
        x[i] = x[i] < 1.0 ? 1.0 : x[i]
        x[i] = x[i] > 2.0 ? 2.0 : x[i]
    end
    return x
end

function get_1D_halton_nodes(m::Int)
    x0 = [1.3; 1.4]
    if m <= 2
        return x0
    end
    m -= 2
    x = Vector{Float64}(undef, m)
    reset_halton()
    init_halton(2)
    for i = 1:m
        x[i] = x0[1] + (x0[2] - x0[1]) * get_halton_node()
    end
    return [x0; x]
end

########### 2D grids

function get_2D_grid(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    m1 = m + 1
    ms = m1^2
    mat = Matrix{Float64}(undef, 2, ms)
    for i = 1:m1
        for j = 1:m1
            r = (i - 1) * m1 + j
            mat[1, r] = x[i]
            mat[2, r] = y[j]
        end
    end
    return mat
end

function get_2D_grid2(m::Int)
    x = collect(range(-1.0, 1.0; step = 2.0/m))
    y = collect(range(-1.0, 1.0; step = 2.0/m))
    m1 = m + 1
    ms = m1^2
    mat = Matrix{Float64}(undef, 2, ms)
    for i = 1:m1
        for j = 1:m1
            r = (i - 1) * m1 + j
            mat[1, r] = x[i]
            mat[2, r] = y[j]
        end
    end
    return mat
end

function get_2D_eps_grid(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    m1 = m + 1
    ms = m1^2
    mat = Matrix{Float64}(undef, 2, ms)
    Random.seed!(123)
    eps = 2.0 / (5.0 * m)
    for i = 1:m1
        for j = 1:m1
            r = (i - 1) * m1 + j
            tx = x[i] + eps * (rand() - 0.5)
            tx = tx < 0.0 ? 0.0 : tx
            tx = tx > 1.0 ? 1.0 : tx
            ty = y[j] + eps * (rand() - 0.5)
            ty = ty < 0.0 ? 0.0 : ty
            ty = ty > 1.0 ? 1.0 : ty
            mat[1, r] = tx
            mat[2, r] = ty
        end
    end
    return mat
end

function get_2D_rect_grid(m::Int)
    x = collect(range(-4.0, 4.0; step = 8.0/m))
    y = collect(range(-2.0, 2.0; step = 4.0/m))
    n1 = length(y)
    n2 = length(x)
    mat = Matrix{Float64}(undef, 2, n1*n2)
    for i = 1:n2
        for j = 1:n1
            r = (i - 1) * n1 + j
            mat[1, r] = x[i]
            mat[2, r] = y[j]
        end
    end
    return mat
end

function get_2D_test1_nodes(m::Int)
    mat0 = [0.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 1.0]'
    if m < 1
        return mat0
    end
    m1 = m + 1
    p = collect(range(1.0/m1, (1.0 - 1.0/m1); step = 1.0/m1))
    ms = m * 4
    mat = Matrix{Float64}(undef, 2, ms)
    for i = 1:m
        mat[1,i] = 0.0
        mat[2,i] = p[i]
    end
    for i = (m+1):(2*m)
        mat[1,i] = 1.0
        mat[2,i] = p[i-m]
    end
    for i = (2*m+1):(3*m)
        mat[1,i] = p[i-2*m]
        mat[2,i] = 0.0
    end
    for i = (3*m+1):(4*m)
        mat[1,i] = p[i-3*m]
        mat[2,i] = 1.0
    end
    w = hcat(mat0, mat)
    return hcat([0.5; 0.5], w)
end

function get_2D_border_nodes(m::Int)
    mat0 = [0.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 1.0]'
    if m < 1
        return mat0
    end
    m1 = m + 1
    p = collect(range(1.0/m1, (1.0 - 1.0/m1); step = 1.0/m1))
    ms = m * 4
    mat = Matrix{Float64}(undef, 2, ms)
    for i = 1:m
        mat[1,i] = 0.0
        mat[2,i] = p[i]
    end
    for i = (m+1):(2*m)
        mat[1,i] = 1.0
        mat[2,i] = p[i-m]
    end
    for i = (2*m+1):(3*m)
        mat[1,i] = p[i-2*m]
        mat[2,i] = 0.0
    end
    for i = (3*m+1):(4*m)
        mat[1,i] = p[i-3*m]
        mat[2,i] = 1.0
    end
    w = hcat(mat0, mat)
    return w
end

"""
Return Lissajous nodes.
"""
function get_2D_Lissajous_nodes(n::Int)
    n1 = n + 1
    n_all = 4 * n * n1
    n_nodes = 2 * n^2 + 4 * n + 1
    nodes = Matrix{Float64}(undef, 2, n_nodes)
    tol = sqrt(eps())
    m = 1
    for k = 1:n_all
        tk = π  * (k - 1) / (2.0 * n * n1)
        x = sin(n * tk)
        y = sin(n1 * tk)
        dup = false
        for i = 1:(m-1)
            if abs(x - nodes[1,i]) <= tol && abs(y - nodes[2,i]) <= tol
                dup = true
                continue
            end
        end
        if !dup
            nodes[1,m] = x
            nodes[2,m] = y
            m += 1
        end
    end
    for k = 1:n_nodes
        nodes[1,k] = (nodes[1,k] + 1.0) / 2.0
        nodes[2,k] = (nodes[2,k] + 1.0) / 2.0
    end
    return nodes
end

function get_2D_halton_nodes(m::Int)
    mat0 = [0.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 1.0]'
    if m <= 4
        return mat0
    end
    m -= 4
    x = Vector{Float64}(undef, m)
    reset_halton()
    init_halton(2)
    for i = 1:m
        x[i] = get_halton_node()
    end
    y = Vector{Float64}(undef, m)
    init_halton(3)
    for i = 1:m
        y[i] = get_halton_node()
    end
    mat = Matrix{Float64}(undef, 2, m)
    for i = 1:m
        mat[1, i] = x[i]
        mat[2, i] = y[i]
    end
    return hcat(mat0, mat)
end

########### 3D grids

function get_3D_grid(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    z = collect(range(0.0, 1.0; step = 1.0/m))
    m1 = m + 1
    ms = m1^3
    mat = Matrix{Float64}(undef, 3, ms)
    for i = 1:m1
        for j = 1:m1
            for k = 1:m1
                r = ((i - 1) * m1 + (j - 1)) * m1 + k
                mat[1, r] = x[i]
                mat[2, r] = y[j]
                mat[3, r] = z[k]
            end
        end
    end
    return mat
end

function get_3D_plot_grid(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    z = [0.0; 0.25; 0.5; 0.75; 1.0]
    mz = length(z)
    m1 = m + 1
    ms = mz * m1^2
    mat = Matrix{Float64}(undef, 3, ms)
    for i = 1:m1
        for j = 1:m1
            for k = 1:mz
                r = ((i - 1) * m1 + (j - 1)) * mz + k
                mat[1, r] = x[i]
                mat[2, r] = y[j]
                mat[3, r] = z[k]
            end
        end
    end
    return mat
end

function get_3D_eps_grid(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    z = collect(range(0.0, 1.0; step = 1.0/m))
    m1 = m + 1
    ms = m1^3
    mat = Matrix{Float64}(undef, 2, ms)
    Random.seed!(123)
    eps = 0.4 / m
    for i = 1:m1
        for j = 1:m1
            for k = 1:m1
                r = ((i - 1) * m1 + (j - 1)) * m1 + k
                tx = x[i] + eps * (rand() - 0.5)
                tx = tx < 0.0 ? 0.0 : tx
                tx = tx > 1.0 ? 1.0 : tx
                ty = y[j] + eps * (rand() - 0.5)
                ty = ty < 0.0 ? 0.0 : ty
                ty = ty > 1.0 ? 1.0 : ty
                tz = z[k] + eps * (rand() - 0.5)
                tz = tz < 0.0 ? 0.0 : tz
                tz = tz > 1.0 ? 1.0 : tz
                mat[1, r] = tx
                mat[2, r] = ty
                mat[3, r] = tz
            end
        end
    end
    return mat
end

function get_3D_halton_nodes(m::Int)
    mat0 = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0; 1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]'
    if m <= 8
        return mat0
    end
    m -= 8
    x = Vector{Float64}(undef, m)
    reset_halton()
    init_halton(2)
    for i = 1:m
        x[i] = get_halton_node()
    end
    y = Vector{Float64}(undef, m)
    init_halton(3)
    for i = 1:m
        y[i] = get_halton_node()
    end
    z = Vector{Float64}(undef, m)
    init_halton(5)
    for i = 1:m
        z[i] = get_halton_node()
    end
    mat = Matrix{Float64}(undef, 3, m)
    for i = 1:m
        mat[1, i] = x[i]
        mat[2, i] = y[i]
        mat[3, i] = z[i]
    end
    return hcat(mat0, mat)
end
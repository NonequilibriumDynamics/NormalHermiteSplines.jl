function get_2D_reg_knots(m::Int)
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

function get_3D_reg_knots(m::Int)
    x = collect(range(0.0, 1.0; step = 1.0/m))
    y = collect(range(0.0, 1.0; step = 1.0/m))
    z = collect(range(0.0, 1.0; step = 1.0/m))
    m1 = m + 1
    ms = m1^3
    mat = Matrix{Float64}(undef, 3, ms)
    for i = 1:m1
        xv = x[i]
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

function get_2D_test_knots(m::Int)
    mat0 = [0.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 1.0]'
    if m <= 4
        return mat0
    end
    m -= 4
    x = Vector{Float64}(undef, m)
    init_halton(2)
    for i = 1:m
        x[i] = get_halton()
    end
    y = Vector{Float64}(undef, m)
    init_halton(3)
    for i = 1:m
        y[i] = get_halton()
    end
    mat = Matrix{Float64}(undef, 2, m)
    for i = 1:m
        mat[1, i] = x[i]
        mat[2, i] = y[i]
    end
    return hcat(mat0, mat)
end

function get_3D_test_knots(m::Int)
    mat0 = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0; 1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]'
    if m <= 8
        return mat0
    end
    m -= 8
    x = Vector{Float64}(undef, m)
    init_halton(2)
    for i = 1:m
        x[i] = get_halton()
    end
    y = Vector{Float64}(undef, m)
    init_halton(3)
    for i = 1:m
        y[i] = get_halton()
    end
    z = Vector{Float64}(undef, m)
    init_halton(5)
    for i = 1:m
        z[i] = get_halton()
    end
    mat = Matrix{Float64}(undef, 3, m)
    for i = 1:m
        mat[1, i] = x[i]
        mat[2, i] = y[i]
        mat[3, i] = z[i]
    end
    return hcat(mat0, mat)
end

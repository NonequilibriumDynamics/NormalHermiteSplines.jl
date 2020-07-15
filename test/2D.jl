@testset "Test 2D" begin

    p = [0. 0. 1. 1. 0.5;  0. 1. 0. 1. 0.5] #  nodes
    u = [0.0; 0.0; 0.0; 0.0; 1.0]   # function values in nodes
    t = [0.5 0.5 0.499999; 0.5 0.499999 0.500001]  # evaluation points

    dp = [0.5 0.5]
    es = [1.0; 0.0; 0.0; 1.0]
    du = [0.]

    @testset "Test 2D-RK_H0 kernel" begin
        spl = prepare(p, RK_H0(0.001))                   # prepare spline
        c = get_cond(spl)                                # get estimation of the problem's Gram matrix condition number
        @test c ≈ 100000.0
        spl = construct(spl, u)                          # construct spline
        vt = [1.0, 3.0]
        σ = evaluate(spl, vt)                            # evaluate spline in the knot
        @test σ ≈ 1.0

        wt = [0.0, 3.0]
        σ1 = evaluate(spl, wt)

        u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
        spl = construct(spl, u2)
        σ2 = evaluate(spl, wt)
        @test σ2 ≈ 2.0 * σ1

        spl = interpolate(p, u, RK_H0(0.001))            # prepare and construct spline
        σ = evaluate(spl, vt)
        @test σ ≈ 1.0
    end

    @testset "Test 2D-1-RK_H1 kernel" begin
        spl = prepare(p, RK_H1(0.001))
        c = get_cond(spl)
        @test c ≈ 1.0e11
        spl = construct(spl, u)
        vt = [1.0, 3.0]
        σ = evaluate(spl, vt)
        @test σ ≈ 1.0

        wt = [0.0, 3.0]
        σ1 = evaluate(spl, wt)

        u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
        spl = construct(spl, u2)
        σ2 = evaluate(spl, wt)
        @test σ2 ≈ 2.0 * σ1

        spl = interpolate(p, u, RK_H1(0.001))
        σ = evaluate(spl, vt)
        @test σ ≈ 1.0
    end

    @testset "Test 2D-1-RK_H2 kernel" begin
        spl = prepare(p, RK_H2(0.001))
        c = get_cond(spl)
        @test c ≈ 1.0e15

        spl = construct(spl, u)
        vt = [1.0, 3.0]
        σ = evaluate(spl, vt)
        @test σ ≈ 1.0

        wt = [0.0, 3.0]
        σ1 = evaluate(spl, wt)

        u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
        spl = construct(spl, u2)
        σ2 = evaluate(spl, wt)
        @test σ2 ≈ 2.0 * σ1

        spl = interpolate(p, u, RK_H2(0.001))
        σ = evaluate(spl, vt)
        @test σ ≈ 1.0
    end


end

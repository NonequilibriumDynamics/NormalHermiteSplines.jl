@testset "Test 2D" begin

    p = [0. 0. 1. 1. 0.5;  0. 1. 0. 1. 0.5] # nodes
    u = [0.0; 0.0; 0.0; 0.0; 1.0]   # function values in nodes
    u2 = [0.0; 0.0; 0.0; 0.0; 2.0]  # second function values in nodes
    t = [0.5 0.5 0.499999; 0.5 0.499999 0.5]  # evaluation points

    dp = [0.5 0.5; 0.5 0.5]
    es = [1.0 0.0; 0.0 1.0]
    du = [0.; 100000.0]

    @testset "Test 2D-RK_H0 kernel" begin
        rk = RK_H0(0.001)
        s = interpolate(p, u, rk)
        σ = evaluate(s, t)
        @test σ[1] ≈ u[5]
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        σ1 = evaluate(s, t[:,1])
        @test σ1[1] ≈ u[5]

        rk = RK_H0()
        s = prepare(p, rk) # prepare spline
        s = construct(s, u) # construct spline
        σ1 = evaluate(s, t) # evaluate spline in points
        @test σ[1] ≈ u[5]
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test σ2[1] ≈ u2[5]
        @test isapprox(σ2[2], u2[5], atol = 1e-5)
        @test isapprox(σ2[3], u2[5], atol = 1e-5)
        @test isapprox(σ2[2], σ2[3], atol = 1e-5)

        cond = get_cond(s)
        @test cond ≈ 10.0

        eps = get_epsilon(s)
        @test isapprox(eps, 1.32, atol = 1e-2)

        est_eps = estimate_epsilon(p) # get estimation of the problem's Gram matrix condition number
        @test isapprox(est_eps, 1.32, atol = 1e-2)
    end

    @testset "Test 2D-RK_H1 kernel" begin
        rk = RK_H1(0.001)
        s = interpolate(p, u, rk)
        cond = get_cond(s)
        @test cond ≈ 1.e11

        σ = evaluate(s, t)
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        grad = evaluate_grad(s, p[:,5])
        @test abs(grad[1]) < 1.0e-5 && abs(grad[2]) < 1.0e-5

        rk = RK_H1()
        s = prepare(p, rk) # prepare spline
        cond = get_cond(s)
        @test cond ≈ 10.0

        s = construct(s, u) # construct spline

        σ1 = evaluate(s, t) # evaluate spline in points
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ2[2], u2[5], atol = 1e-5)
        @test isapprox(σ2[3], u2[5], atol = 1e-5)
        @test isapprox(σ2[2], σ2[3], atol = 1e-5)

        est_eps = estimate_epsilon(p, dp) # get estimation of the problem's Gram matrix condition number
        @test est_eps ≈ 1.37 atol = 1e-2
###
        eps = 0.0001
        rk = RK_H1(eps)
        s = interpolate(p, u, dp, es, du, rk)
        cond = get_cond(s)
        @test cond == 1.0e14

        σ = evaluate(s, p[:,5])
        @test !isapprox(σ[1], u[5], atol = 0.1)

        q = estimate_interpolation_quality(s)
        @test q > 1.0

# Same test with extended precision
        rk = RK_H1(Double64(eps))
        p = Double64.(p)
        dp = Double64.(dp)
        es = Double64.(es)
        u = Double64.(u)
        du = Double64.(du)
        t = Double64.(t)
        s = prepare(p, dp, es, rk)
        s = construct(s, u, du)
        σ = evaluate(s, p)
        @test all(isapprox.(σ, u, atol = 1e-5))

        q = estimate_interpolation_quality(s)
        @test q < 1e-5
    end

    @testset "Test 2D-RK_H2 kernel" begin
        p = [0. 0. 1. 1. 0.5;  0. 1. 0. 1. 0.5] # nodes
        u = [0.0; 0.0; 0.0; 0.0; 1.0]   # function values in nodes
        u2 = [0.0; 0.0; 0.0; 0.0; 2.0]  # second function values in nodes
        t = [0.5 0.5 0.499999; 0.5 0.499999 0.5]  # evaluation points

        dp = [0.5 0.5; 0.5 0.5]
        es = [1.0 0.0; 0.0 1.0]
        du = [0.; 100000.0]

        rk = RK_H2(0.01)
        s = interpolate(p, u, rk)
        cond = get_cond(s)
        @test cond ≈ 1e11

        σ = evaluate(s, t)
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        grad = evaluate_grad(s, p[:,5])
        @test abs(grad[1]) < 1.0e-5 && abs(grad[2]) < 1.0e-5

        rk = RK_H2()
        s = prepare(p, rk) # prepare spline
        cond = get_cond(s)
        @test cond ≈ 100.0

        s = construct(s, u) # construct spline

        σ1 = evaluate(s, t) # evaluate spline in points
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ[2], u[5], atol = 1e-5)
        @test isapprox(σ[3], u[5], atol = 1e-5)
        @test isapprox(σ[2], σ[3], atol = 1e-5)

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test isapprox(σ[1], u[5], atol = 1e-5)
        @test isapprox(σ2[2], u2[5], atol = 1e-5)
        @test isapprox(σ2[3], u2[5], atol = 1e-5)
        @test isapprox(σ2[2], σ2[3], atol = 1e-5)
###
        eps = 0.01
        rk = RK_H2(eps)
        s = interpolate(p, u, dp, es, du, rk)
        cond = get_cond(s)
        @test cond == 1.0e13

        σ = evaluate(s, p[:,5])
        @test !isapprox(σ[1], u[5], atol = 0.1)

        q = estimate_interpolation_quality(s)
        @test q > 1.0

# Same tests with extended precision
        rk = RK_H2(Double64(eps))
        p = Double64.(p)
        dp = Double64.(dp)
        es = Double64.(es)
        u = Double64.(u)
        du = Double64.(du)
        t = Double64.(t)
        s = prepare(p, dp, es, rk)
        s = construct(s, u, du)

        σ = evaluate(s, p)
        @test all(isapprox.(σ, u, atol = 1e-5))

        q = estimate_interpolation_quality(s)
        @test q < 1e-5
    end
end

@testset "Test 2D-Bis" begin

    p = collect([-1.0 2.0; -1.0 4.0; 3.0 2.0; 3.0 4.0; 1.0 3.0]') # function knots
    u = [0.0; 0.0; 0.0; 0.0; 1.0]                        # function values in knots

    t = collect([-1.0 3.0; 0.0 3.0; 1.0 3.0; 2.0 3.0; 3.0 3.0]')  # evaluation points

    @testset "Test 2D-Bis-RK_H0 kernel" begin
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

    @testset "Test 2D-Bis-RK_H1 kernel" begin
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

    @testset "Test 2D-Bis-RK_H2 kernel" begin
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

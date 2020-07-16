@testset "Test 2D" begin

    p = [0. 0. 1. 1. 0.5;  0. 1. 0. 1. 0.5] # nodes
    u = [0.0; 0.0; 0.0; 0.0; 1.0]   # function values in nodes
    u2 = [0.0; 0.0; 0.0; 0.0; 2.0]  # function values in nodes (2)
    t = [0.5 0.5 0.499999; 0.5 0.499999 0.5]  # evaluation points

    dp = [0.5 0.5; 0.5 0.5]
    es = [1.0 0.0; 0.0 1.0]
#    du = [0.; 0.0]
    du = [0.; 10000.0]

    @testset "Test 2D-RK_H0 kernel" begin
        rk = RK_H0(0.001)
        s = interpolate(p, u, rk)
        σ = evaluate(s, t)
        @test σ[1] ≈ u[5]
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))

        σ1 = evaluate(s, t[:,1])
        @test σ1[1] ≈ u[5]

        rk = RK_H0()
        s = prepare(p, rk) # prepare spline
        s = construct(s, u) # construct spline
        σ1 = evaluate(s, t) # evaluate spline in points
        @test σ[1] ≈ u[5]
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test σ2[1] ≈ u2[5]
        @test all(isapprox.(σ2[2], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[3], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[2], σ2[3], atol = 1e-5))

        cond = get_cond(s)
        @test cond ≈ 10.0

        eps = get_epsilon(s)
        @test all(isapprox.(eps, 1.32, atol = 1e-2))

        est_eps = estimate_epsilon(p) # get estimation of the problem's Gram matrix condition number
        @test all(isapprox.(est_eps, 1.32, atol = 1e-2))
    end

    @testset "Test 2D-RK_H1 kernel" begin
        rk = RK_H1(0.001)
        s = interpolate(p, u, rk)
        σ = evaluate(s, t)
#        println(σ)
        @test all(isapprox.(σ[1], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))

        rk = RK_H1()
        s = prepare(p, rk) # prepare spline
        s = construct(s, u) # construct spline
        σ1 = evaluate(s, t) # evaluate spline in points
#        println(σ1)

        @test all(isapprox.(σ[1], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test all(isapprox.(σ[1], u[5], atol = 1e-5))
        @test all(isapprox.(σ2[2], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[3], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[2], σ2[3], atol = 1e-5))

        σ3 = evaluate(s, p[:,5])
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))

        cond = get_cond(s)
        @test cond ≈ 100.0

        est_eps = estimate_epsilon(p, dp) # get estimation of the problem's Gram matrix condition number
        @test all(isapprox.(est_eps, 1.37, atol = 1e-2))
###
        rk = RK_H1(0.1)
        s = interpolate(p, u, dp, es, du, rk)
        # s = prepare(p, dp, es, rk)
        # s = construct(s, u, du)
        cond = get_cond(s)
        σ = evaluate(s, p[:,5])
        println("point #1")
        println(cond)
        println(σ)

        σ = evaluate(s, t)
        println("point #2")
        println(σ)

        q = get_interpolation_quality(s, p, u)
        println("point #3")
        println(q)
        println(p)
        println(u)
        println(evaluate(s, p))

# Same test with extended precision
        rk = RK_H1(Double64(0.1))
        p = Double64.(p)
        dp = Double64.(dp)
        es = Double64.(es)
        u = Double64.(u)
        du = Double64.(du)
        t = Double64.(t)
        s = interpolate(p, u, dp, es, du, rk)
        # s = prepare(p, dp, es, rk)
        # s = construct(s, u, du)
        cond = get_cond(s)
        σ = evaluate(s, p[:,5])
        println("point #4")
        println(cond)
        println(σ)

        σ = evaluate(s, t)
        println("point #5")
        println(σ)

        q = get_interpolation_quality(s, p, u)
        println("point #6")
        println(q)
        println(p)
        println(u)
        println(evaluate(s, p))


        # @test σ[1] ≈ u[5]
        # @test all(isapprox.(σ[1], u[5], atol = 1e-5))

        # σ = evaluate(s, t)
        # println(σ)
        # @test all(isapprox.(σ[1], u[5], atol = 1e-5))
        # @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        # @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        # @test all(isapprox.(σ[2], σ[3], atol = 1e-5))


    end



end

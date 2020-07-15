@testset "Test 2D" begin

    p = [0. 0. 1. 1. 0.5;  0. 1. 0. 1. 0.5] # nodes
    u = [0.0; 0.0; 0.0; 0.0; 1.0]   # function values in nodes
    u2 = [0.0; 0.0; 0.0; 0.0; 2.0]  # function values in nodes (2)
    t = [0.5 0.5 0.499999; 0.5 0.499999 0.5]  # evaluation points

    dp = [0.5 0.5]
    es = [1.0; 0.0; 0.0; 1.0]
    du = [0.; -1000.0]

    @testset "Test 2D-RK_H0 kernel" begin
        rk = RK_H0(0.001)
        s = interpolate(p, u, rk)
        σ = evaluate(s, t)
        @test σ[1] ≈ u[5]
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))
        println(σ)

        rk = RK_H0()
        s = prepare(p, rk)
        s = construct(s, u)
        σ1 = evaluate(s, t)
        @test σ[1] ≈ u[5]
        @test all(isapprox.(σ[2], u[5], atol = 1e-5))
        @test all(isapprox.(σ[3], u[5], atol = 1e-5))
        @test all(isapprox.(σ[2], σ[3], atol = 1e-5))
        println(σ)

        s = construct(s, u2)
        σ2 = evaluate(s, t)
        @test σ2[1] ≈ u2[5]
        @test all(isapprox.(σ2[2], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[3], u2[5], atol = 1e-5))
        @test all(isapprox.(σ2[2], σ2[3], atol = 1e-5))
        println(σ2)

    #     spl = prepare(p, RK_H0(0.001))                   # prepare spline
    #     c = get_cond(spl)                                # get estimation of the problem's Gram matrix condition number
    #     @test c ≈ 100000.0
    #     spl = construct(spl, u)                          # construct spline
    #     vt = [1.0, 3.0]
    #     σ = evaluate(spl, vt)                            # evaluate spline in the knot
    #     @test σ ≈ 1.0
    #
    #     wt = [0.0, 3.0]
    #     σ1 = evaluate(spl, wt)
    #
    #     u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
    #     spl = construct(spl, u2)
    #     σ2 = evaluate(spl, wt)
    #     @test σ2 ≈ 2.0 * σ1
    #
    #     spl = interpolate(p, u, RK_H0(0.001))            # prepare and construct spline
    #     σ = evaluate(spl, vt)
    #     @test σ ≈ 1.0
     end
    #
    # @testset "Test 2D-1-RK_H1 kernel" begin
    #     spl = prepare(p, RK_H1(0.001))
    #     c = get_cond(spl)
    #     @test c ≈ 1.0e11
    #     spl = construct(spl, u)
    #     vt = [1.0, 3.0]
    #     σ = evaluate(spl, vt)
    #     @test σ ≈ 1.0
    #
    #     wt = [0.0, 3.0]
    #     σ1 = evaluate(spl, wt)
    #
    #     u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
    #     spl = construct(spl, u2)
    #     σ2 = evaluate(spl, wt)
    #     @test σ2 ≈ 2.0 * σ1
    #
    #     spl = interpolate(p, u, RK_H1(0.001))
    #     σ = evaluate(spl, vt)
    #     @test σ ≈ 1.0
    # end
    #
    # @testset "Test 2D-1-RK_H2 kernel" begin
    #     spl = prepare(p, RK_H2(0.001))
    #     c = get_cond(spl)
    #     @test c ≈ 1.0e15
    #
    #     spl = construct(spl, u)
    #     vt = [1.0, 3.0]
    #     σ = evaluate(spl, vt)
    #     @test σ ≈ 1.0
    #
    #     wt = [0.0, 3.0]
    #     σ1 = evaluate(spl, wt)
    #
    #     u2 = [0.0; 0.0; 0.0; 0.0; 2.0]
    #     spl = construct(spl, u2)
    #     σ2 = evaluate(spl, wt)
    #     @test σ2 ≈ 2.0 * σ1
    #
    #     spl = interpolate(p, u, RK_H2(0.001))
    #     σ = evaluate(spl, vt)
    #     @test σ ≈ 1.0
    # end


end

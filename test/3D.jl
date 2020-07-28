@testset "Test 3D" begin

    p = [-1.0 2.0; -1.0 4.0; 3.0 2.0; 3.0 4.0; 1.0 3.0]' # function knots
    u = [0.0; 0.0; 0.0; 0.0; 1.0]                        # function values in knots

    t = [-1.0 3.0; 0.0 3.0; 1.0 3.0; 2.0 3.0; 3.0 3.0]'  # evaluation points

    @testset "Test 3D-RK_H0 kernel" begin
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

    @testset "Test 3D-RK_H1 kernel" begin
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

    @testset "Test 3D-RK_H2 kernel" begin
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

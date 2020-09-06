export demo

using Gadfly

function demo(type_of_kernel::Int = 1)
# Setup
    x = collect(1.0:1.0:20)
    u = x.*0.0
    s = x
    v = x.*0.0

    for i in 6:10
        u[i] = 1.0
    end
    for i in 11:14
        u[i] = -0.2 * i + 3.0
    end
    for i in 11:14
        v[i] = -0.2
    end

    p = collect(1.0:0.2:20)
    r = p.*0.0
    tol = 1e-15
    i = 0
    for pi in p
        i += 1
        if pi >= (6.0 - tol) && pi <= (10.0 + tol)
            r[i] = 1.0
        end
        if pi > (10.0 + tol) && pi < (15.0 - tol)
            r[i] = -0.2 * p[i] + 3.0
        end
    end
####
    rk = RK_H1()
    if type_of_kernel == 0
        rk = RK_H0()
    elseif type_of_kernel == 2
        rk = RK_H2()
    end
    spline = prepare(x, rk)
    cond = get_cond(spline)
    ε = get_epsilon(spline)
    println("cond = $cond, ε = $ε")

    spline = construct(spline, u)
    σ = evaluate(spline, p)

    δ = σ .- r
    rmse = round(norm(δ)/sqrt(length(δ)); digits=2)
    println("rmse = $rmse")

    set_default_plot_size(13cm, 13cm)
    plt = Gadfly.plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
          layer(x = p, y = r, Geom.line, Theme(default_color=colorant"red")),
          layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
          Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
          Guide.manual_color_key("Legend", ["Points", "True", "Spline"], ["orange", "red", "blue"]),
          Guide.xlabel("nodes, points"), Guide.ylabel("f, σ"),
          Guide.title("Fig.1a"))

    if type_of_kernel != 0
        Gadfly.draw(SVG("c:/0/example-1a.svg", 13cm, 13cm), plt)

        dσ = similar(σ)
        for i=1:length(p)
            dσ[i] = evaluate_derivative(spline, p[i])
        end
        plt = Gadfly.plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
              layer(x = p, y = dσ, Geom.line, Theme(default_color=colorant"green")),
              Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
              Guide.manual_color_key("Legend", ["Points", "Spline 1st derivative"], ["orange", "green"]),
              Guide.xlabel("nodes, points"), Guide.ylabel(raw"σ'"),
              Guide.title("Fig.2a"))
        Gadfly.draw(SVG("c:/0/example-1a-der.svg", 13cm, 13cm), plt)
    else
        Gadfly.draw(SVG("c:/0/example-1c.svg", 13cm, 13cm), plt)
    end

    σ = evaluate(spline, [3.1, 8.1, 18.1])

    u2 = 2.0 .* u
    spline = construct(spline, u2)
    σ = evaluate(spline, p)
    σ = evaluate(spline, [3.1, 8.1, 18.1])


###
      spline = interpolate(x, u, s, v, RK_H1())

      cond = get_cond(spline)
      ε = get_epsilon(spline)
      println("cond = $cond, ε = $ε")

      σ = evaluate(spline, p)
      δ = σ .- r
      rmse = round(norm(δ)/sqrt(length(δ)); digits=2)
      println("rmse = $rmse")

      set_default_plot_size(13cm, 13cm)
      plt = Gadfly.plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
            layer(x = p, y = r, Geom.line, Theme(default_color=colorant"red")),
            layer(x = p, y = σ, Geom.line, Theme(default_color=colorant"blue")),
            Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
            Guide.manual_color_key("Legend", ["Points", "True", "Spline"], ["orange", "red", "blue"]),
            Guide.xlabel("nodes, points"), Guide.ylabel("f, σ"),
            Guide.title("Fig.1b"))
      Gadfly.draw(SVG("c:/0/example-1b.svg", 13cm, 13cm), plt)

      dσ = similar(σ)
      for i=1:length(p)
          dσ[i] = evaluate_derivative(spline, p[i])
      end
      plt = Gadfly.plot(layer(x = x, y = u, Geom.point, Theme(default_color=colorant"orange")),
            layer(x = p, y = dσ, Geom.line, Theme(default_color=colorant"green")),
            Scale.y_continuous(minvalue=-0.5, maxvalue=1.5),
            Guide.manual_color_key("Legend", ["Points", "Spline 1st derivative"], ["orange", "green"]),
            Guide.xlabel("nodes, points"), Guide.ylabel(raw"σ'"),
            Guide.title("Fig.2b"))
      Gadfly.draw(SVG("c:/0/example-1b-der.svg", 13cm, 13cm), plt)

end

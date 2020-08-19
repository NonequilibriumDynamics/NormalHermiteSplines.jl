function get_1D_model1(p::Float64)
    val = 0.0
    t = abs(3.0*(1.0 - p) + 1.0)
    if t > 0.0
        val = 1.0 + p^2 + log(t)/3.3
    end
    return val
end

function get_1D_model1_grad(p::Float64, tol::Float64=1.e-15)
    val = 0.0
    if abs(p - 4.0/3.0) > tol
        val = 2*p + 1.0/(3.3 * (p - 4.0/3.0))
    end
end

function get_2D_model1(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    val = 0.0
    x = p[1]
    y = p[2]
    d = (x - 0.75)^2 + (y - 0.75)^2
    r2 = 0.04
    if d <= r2
        val = 1.0
    end
    if x <= 0.25 && x >= 0.05 && y <= 0.25 && y >= 0.05
        val = 1.0
    end
    return val
end

function get_2D_model2(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    f = 0.0
    x = p[1]
    y = p[2]
    if x > 1.0 || x < 0.0 || y > 1.0 || y < 0.0
        return f
    end

    if y <= x && y <= (1.0 - x)
        f = 2.0 * y
    end
    if y <= x && y >= (1.0 - x)
        f = 1.0 - 2.0 * (x - 0.5)
    end
    if y >= x && y >= (1.0 - x)
        f = 1.0 - 2.0 * (y - 0.5)
    end
    if y >= x && y <= (1.0 - x)
        f = 2.0 * x
    end
    return f
end

function get_2D_model2_grad(p::Vector{Float64})
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]

    if y < x && y < (1.0 - x)
        grad[1] = 0.0
        grad[2] = 2.0
    end
    if y < x && y > (1.0 - x)
        grad[1] = -2.0
        grad[2] = 0.0
    end
    if y > x && y > (1.0 - x)
        grad[1] = 0.0
        grad[2] = -2.0
    end
    if y > x && y < (1.0 - x)
        grad[1] = 2.0
        grad[2] = 0.0
    end
    return grad
end

# Franke's CLIFF Function F2 (Cliff function) from
# R. Renka, R. Brown, Algorithm 792: accuracy test of ACM algorithms for interpolation of scattered data in the plane, ACM Transactions on Mathematical Software (TOMS), Vol.25, No.1,1999
function get_2D_model3(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    x = p[1]
    y = p[2]
    f = (tanh(9.0 * (y - x)) + 1.0) / (tanh(9.0)+ 1.0)
    return f
end

function get_2D_model3_grad(p::Vector{Float64})
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]

    grad[1] = (-9.0 / (tanh(9.0)+ 1.0)) * (sech(9.0 * (x - y)))^2
    grad[2] = (9.0 / (tanh(9.0)+ 1.0)) * (sech(9.0 * (x - y)))^2
    return grad
end

# Function F10 from
# R. Renka, R. Brown, Algorithm 792: accuracy test of ACM algorithms for interpolation of scattered data in the plane, ACM Transactions on Mathematical Software (TOMS), Vol.25, No.1,1999
function get_2D_model4(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    x = p[1]
    y = p[2]
    d = sqrt((80.0 * x - 40.0)^2 + (90.0 * y - 45.0)^2)
    f = exp(-0.04 * d) * cos(0.15 * d)
    return f
end

function get_2D_model4_grad(p::Vector{Float64})
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]
    d = sqrt((80.0 * x - 40.0)^2 + (90.0 * y - 45.0)^2)
    if d > 0.0
        grad[1] = exp(-0.04 * d) * (-(960.0 * x - 480.0) * sin(0.15 * d) - (256.0 * x - 128.0) * cos(0.15 * d)) / d
        grad[2] = exp(-0.04 * d) * (-(1215.0 * y - 607.5) * sin(0.15 * d) - (324.0 * x - 162.0) * cos(0.15 * d)) / d
    end
    return grad
end

# Function F2 from
# S. De Marchi, A. Martinez, E. Perracchione, Fast and stable rational RBF-based partition of unity interpolation. Journal of Computational and Applied Mathematics, 2018.
function get_2D_model5(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    x = p[1]
    y = p[2]
    f = (x + y - 1.0)^7
    return f
end

function get_2D_model5_grad(p::Vector{Float64})
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]
    grad[1] = grad[2] = 7.0 * (x + y - 1.0)^6
    return grad
end

function get_2D_model6(p::Vector{Float64})
#  (x,y) \in [-1;1][-1;1]
    x = p[1]
    y = p[2]
    f = sin(4.0 * sqrt(x^2 + y^2))
    return f
end

function get_2D_model6_grad(p::Vector{Float64})
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]
    d = sqrt(x^2 + y^2)
    if d > 0.0
        grad[1] = 4.0 * x * cos(4.0 * d) / d
        grad[2] = 4.0 * y * cos(4.0 * d) / d
    end
    return grad
end

function get_2D_model10(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    return 1.0
end

function get_2D_model11(p::Vector{Float64})
#  (x,y) \in [0;1][0;1]
    val = 0.0
    x = p[1]
    y = p[2]
    d = (x - 0.5)^2 + (y - 0.5)^2
    r2 = 0.16 * (1.0 + 100.0 * eps(1.0))
    if d <= r2
        val = sqrt(abs(0.16 - d))
    end
    return val
end

function get_2D_model12(p::Vector{Float64})
#  (x,y) \in [-4;4][-2;2]
    val = 0.0
    x = p[1]
    y = p[2]
    d = x^2 + y^2
    r2 = 1.0 + 100.0 * eps(1.0)
    if d <= r2
        val = sqrt(abs(1.0 - d))
    end
    d = (x - 2.5)^2 + y^2
    if d <= r2
        val = sqrt(abs(1.0 - d))
    end
    d = (x + 2.5)^2 + y^2
    if d <= r2
        val = sqrt(abs(1.0 - d))
    end
    return val
end

function get_2D_model12_grad(p::Vector{Float64})
#  (x,y) \in [-4;4][-2;2]
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]
    d = x^2 + y^2
    r2 = 1.0 + 100.0 * eps(1.0)
    if d <= r2 && abs(1.0 - d) > eps()
        grad[1] =  -x / sqrt(abs(1.0 - d))
        grad[2] =  -y / sqrt(abs(1.0 - d))
    end
    d = (x - 2.5)^2 + y^2
    if d <= r2 && abs(1.0 - d) > eps()
        grad[1] =  -(x - 2.5) / sqrt(abs(1.0 - d))
        grad[2] =  -y / sqrt(abs(1.0 - d))
    end
    d = (x + 2.5)^2 + y^2
    if d <= r2 && abs(1.0 - d) > eps()
        grad[1] =  -(x + 2.5) / sqrt(abs(1.0 - d))
        grad[2] =  -y / sqrt(abs(1.0 - d))
    end
    return grad
end

function get_2D_model13(p::Vector{Float64})
# Franke's EXPONENTIAL f1
# R. Renka, R. Brown, Algorithm 792: accuracy test of ACM algorithms for interpolation of scattered data in the plane, ACM Transactions on Mathematical Software (TOMS), Vol.25, No.1,1999
# https://dl.acm.org/doi/10.1145/305658.305745
#  (x,y) \in [0;1][0;1]
    x = p[1]
    y = p[2]
    val = 0.75*exp(-((9.0* x - 2.)^2 + (9.0*y - 2.)^2)/4.0) +
          0.75*exp(-((9.0*x + 1.0)^2)/49. - (9.0*y + 1.0)/10.0) +
          0.5*exp(-((9.0*x - 7.)^2 + (9.0*y - 3.)^2)/4.0) -
          0.2*exp(-(9.0*x - 4.)^2 - (9.0*y - 7.)^2)
    return val
end

function get_2D_model13_grad(p::Vector{Float64})
    #  (x,y) \in [0;1][0;1]
    grad = [0.0; 0.0]
    x = p[1]
    y = p[2]
    grad[1] = -0.27551*(9.0*x + 1.0)*exp(0.1*(-9.0*y - 1.0) - (1.0/49.0)*(9.0*x + 1.0)^2) -
               2.25*(9.0*x - 7.)*exp(0.25*(-(9.0*x - 7.)^2 - (9.0*y - 3.)^2)) +
               3.6*(9.0*x - 4.)*exp(-(9.0*x - 4.)^2 - (9.0*y - 7.)^2) -
               3.375*(9.0*x - 2.)*exp(0.25*(-(9.0*x - 2.)^2 - (9.0*y - 2.)^2))
    grad[2] = -0.675*exp(0.1*(-9.0*y - 1.0) - (1.0/49.0)*(9.0*x + 1.0)^2) +
               3.6*(9.0*y - 7.)*exp(-(9.0*x - 4.)^2 - (9.0*y - 7.)^2) -
               2.25*(9.0*y - 3.)*exp(0.25*(-(9.0*x - 7.)^2 - (9.0*y - 3.)^2)) -
               3.375*(9.0*y - 2.)*exp(0.25*(-(9.0*x - 2.)^2 - (9.0*y - 2.)^2))
    return grad
end

function get_3D_model1(p::Vector{Float64})
    r = p[1] + p[2] + p[3]
    val = -1.0
    if(p[1] < 0.0 || p[2] < 0.0 || p[3] < 0.0 || r > 1.0)
        return val
    end
    return r
end

function get_3D_model1_grad()
    return [1.0; 1.0; 1.0]
end

function get_3D_model2(p::Vector{Float64})
    return p[3]
end

function get_3D_model2_grad()
    return [0.0; 0.0; 1.0]
end


function get_3D_model3(p::Vector{Float64})
# T. Foley, Interpolation and approximation of 3-D and 4-D scattered data, Comput. Math. Appl., Vol.13, No.8, 1987.
# https://www.sciencedirect.com/science/article/pii/0898122187900435
# G_1(x,y,z), (x,y,z) \in [0;1][0;1][0;1]
    x = p[1]
    y = p[2]
    z = p[3]

    val = 0.75*exp(-16.0*((x - 0.25)^2 + (y - 0.25)^2 + (z - 0.25)^2)) +
          0.50*exp(-10.0*((x - 0.25)^2 + (y - 0.25)^2)) +
          0.50*exp(-10.0*((x - 0.75)^2 + (y - 0.125)^2 + (z - 0.5)^2)) -
          0.25*exp(-20*((x - 0.75)^2 + (y - 0.75)^2))
    return val
end

function get_3D_model3_grad(p::Vector{Float64})
    grad = [0.0; 0.0; 0.0]
    x = p[1]
    y = p[2]
    z = p[3]
    grad[1] = -24.0*exp(-16.0*((-0.25 + x)^2 + (-0.25 + y)^2 + (-0.25 + z)^2))*(-0.25 + x) -
              10.0*exp(-10.0*((-0.25 + x)^2 + (-0.25 + y)^2))*(-0.25 + x) -
              10.0*exp(-10.0*((-0.75 + x)^2 + (-0.125 + y)^2 + (-0.5 + z)^2))*(-0.75 + x) +
              10.0*exp(-20.0*((-0.75 + x)^2 + (-0.75 + y)^2))*(-0.75 + x)
    grad[2] = -24.0*exp(-16.0*((-0.25 + x)^2 + (-0.25 + y)^2 + (-0.25 + z)^2))*(-0.25 + y) -
              10.0*exp(-10.0*((-0.25 + x)^2 + (-0.25 + y)^2))*(-0.25 + y) -
              10.0*exp(-10.0*((-0.75 + x)^2 + (-0.125 + y)^2 + (-0.5 + z)^2))*(-0.125 + y) +
              10.0*exp(-20.0*((-0.75 + x)^2 + (-0.75 + y)^2))*(-0.75 + y)
    grad[3] = -24.0*exp(-16.0*((-0.25 + x)^2 + (-0.25 + y)^2 + (-0.25 + z)^2))*(-0.25 + z) -
              10.0*exp(-10.0*((-0.75 + x)^2 + (-0.125 + y)^2 + (-0.5 + z)^2))*(-0.5 + z)
    return grad
end

function get_3D_model4(p::Vector{Float64})
# T. Foley, Interpolation and approximation of 3-D and 4-D scattered data, Comput. Math. Appl., Vol.13, No.8, 1987.
# https://www.sciencedirect.com/science/article/pii/0898122187900435
# G_2(x,y,z), (x,y,z) \in [0;1][0;1][0;1]
    val = get_3D_model3(p)
    if val < 0.0
        val = 0.0
    end
    if val > 0.5
        val = 0.5
    end
    return val
end

function get_3D_model4_grad(p::Vector{Float64})
    grad = get_3D_model3_grad(p)
    val = get_3D_model3(p)
    if val < 0.0
        grad = [0.0; 0.0; 0.0]
    end
    if val > 0.5
        grad = [0.0; 0.0; 0.0]
    end
    return grad
end

# Function F2 (Cliff function) from
# R. Renka, Multivariate interpolation of large sets of scattered data. ACM Transactions on Mathematical Software, Vol.14, No.2, 1988.
function get_3D_model5(p::Vector{Float64})
# Franke's SADDLE F3
#  (x,y,z) \in [0;1][0;1][0;1]
    x = p[1]
    y = p[2]
    z = p[3]
    f = (1.25 + cos(5.4*y))*cos(6.0*z)/(6.0 + 6.0*(3.0*x - 1.0)^2)
    return f
end

function get_3D_model5_grad(p::Vector{Float64})
    grad = [0.0; 0.0; 0.0]
    x = p[1]
    y = p[2]
    z = p[3]
    grad[1] = -((108.0*x - 36.0)*(cos(5.4*y) + 1.25)*cos(6.0*z))/(6.0 + 6.0*(3.0*x - 1.0)^2)^2
    grad[2] = -(5.4*sin(5.4*y)*cos(6.0*z))/(6.0 + 6.0*(3.0*x - 1.0)^2)
    grad[3] = -(cos(5.4*y) + 1.25)*sin(6.0*z)/(1.0 + (3.0*x - 1.0)^2)
    return grad
end



# ROBERT E. BARNHILL and SARAH E. STEAD, Multistage trivariate surfaces, The Rocky Mountain Journal of Mathematics, Vol. 14, No. 1, 1984.
# https://www.jstor.org/stable/44236789?seq=1
# http://web.eecs.utk.edu/research/imp/tools/interp/

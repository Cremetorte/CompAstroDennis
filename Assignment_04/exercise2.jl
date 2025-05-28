# symmetric numeric differentiation
function d_dx(f)
    h = eps(Float64)^(1/3) # ³√ϵ as h
    res(x) = (f(x + h) - f(x - h)) / (2 * h)
    return res  # type: function
end

# i-th derivative using d_dx
function ith_derivative(i, f)
    g = f
    for i in 1:i    # i included!
        g = d_dx(g)
    end
    return g   # type: Function
end

# Error of Simpson Rule
function simpson_error(f, a, b)
    interval = LinRange(a, b, 100)
    d⁴f_dx⁴ = ith_derivative(4, f)  # cursed name
    ξ = maximum(d⁴f_dx⁴, interval)  # Worst case for error
    domain = b - a
    return d⁴f_dx⁴(ξ) * (domain/2)^5 / 90    
end

# Trapez rule
function trapezoidal_integral(f, a, b)
    domain = b - a
    # as in lecture
    return domain/2 * (f(a) + f(b))
end

#Simpson rule
function simpsons_integral(f, a, b)
    domain = b - a
    # as in lecture
    return domain/6 * (f(a) + 4*f((a+b)/2) + f(b))
end

# adaptive simpson rule
function adaptive_simpson_integral(f, a, b, epsilon=1e-8)
    if simpson_error(f, a, b) < epsilon
        # if error is small enough, return the actual numerical integral
        return simpsons_integral(f, a, b)
    else
        # cut interval in two equal halves
        a_new = a
        b_new = (a+b)/2
        c_new = b

        # sum over these two halves
        return adaptive_simpson_integral(f, a_new, b_new, epsilon) + adaptive_simpson_integral(f, b_new, c_new, epsilon)
    end
end



function si_simpson(x)
    f(t) = sinc(t)
    return adaptive_simpson_integral(f, 0, x)
end

function C_simpson(x)
    f(t) = cos(π/2 * t^2)
    return adaptive_simpson_integral(f, 0, x)
end

println("Using the adaptive Simpson's rule: \nSi(1) = $(si_simpson(1)) \nC(5) = $(C_simpson(5))")


# Trapez Method adaptive

function subdivide(a, b, N)
    if N == 1
        return [[a,b]]
    end
    interval = LinRange(a, b, N+1)
    res = []
    for i in 1:N
        push!(res, [interval[i], interval[i+1]])
    end
    return res
end

function adaptive_trapezoidal_rule(f, a, b, epsilon=1e-8, max_steps=1e10)
    res = trapezoidal_integral(f, a, b)

    N = 2 # one was already done
    for i in 1:max_steps-1
        subdivisions = subdivide(a,b,N)
        temp_res = 0
        for i in subdivisions
            temp_res += trapezoidal_integral(f, i[1], i[2])
        end

        if abs(temp_res - res) < epsilon
            return temp_res
        end

        N *= 2
    end
    error("Adaptive trapezoidal rule did not converge within the maximum number of steps.")
end



function si_trapez(x)
    f(t) = sinc(t)
    return adaptive_trapezoidal_rule(f, 0, x)
end

function C_trapez(x)
    f(t) = cos(π/2 * t^2)
    return adaptive_trapezoidal_rule(f, 0, x)
end

println("Using the adaptive Trapez rule: \nSi(1) = $(si_trapez(1)) \nC(5) = $(C_trapez(5))")

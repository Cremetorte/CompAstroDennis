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
    interval = LinRange(a, b, 400)
    d⁴f_dx⁴ = ith_derivative(4, f)  # cursed name
    max_d4dx = maximum(d⁴f_dx⁴.(interval))  # Worst case for error <=> ξ is chosen to maximize error
    domain = b - a
    return max_d4dx * (domain/2)^5 / 90    
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
    f(t) = sinc(t/π)   # sinc function (Julias definition uses sinc(x) = sin(πx)/πx)
    return adaptive_simpson_integral(f, 0, x)
end

function C_simpson(x)
    f(t) = cos(π/2 * t^2)
    return adaptive_simpson_integral(f, 0, x)
end

println("Using the adaptive Simpson's rule: \nSi(1) = $(si_simpson(1)) \nC(5) = $(C_simpson(5))")


# Trapez Method adaptive

# helper
function trapez_fixed_nr_of_steps(f, a, b, N)
    nodes = LinRange(a, b, N+1)
    res = 0.0
    for i in 1:N
        #use linrange as nodes
        res += trapezoidal_integral(f, nodes[i], nodes[i+1])
    end
    return res
end

function adaptive_trapezoidal_rule(f, a, b, epsilon=1e-8, max_steps=1e4)
    res = trapezoidal_integral(f, a, b)

    N = 2 # one was already done
    for i in 1:max_steps-1
        temp_res = trapez_fixed_nr_of_steps(f, a, b, N)

        if abs(temp_res - res) < epsilon
            return temp_res
        end
        res = temp_res

        N *= 2
    end
    error("Adaptive trapezoidal rule did not converge within the maximum number of steps.")
end






function si_trapez(x)
    f(t) = sinc(t/π)   # sinc function (Julias definition uses sinc(x) = sin(πx)/πx)
    return adaptive_trapezoidal_rule(f, 0, x)
end

function C_trapez(x)
    f(t) = cos(π/2 * t^2)
    return adaptive_trapezoidal_rule(f, 0, x)
end

println("Using the \"adaptive Trapez rule\": \nSi(1) = $(si_trapez(1)) \nC(5) = $(C_trapez(5))")


try
    # Compare to Scipys result
    using PyCall    # Has to be installed via REPL -> Pkg mode (press "]") -> add PyCall

    scipy = pyimport("scipy") # apparently this is possible 

    # Second funtion; sinc is already defined in Julia
    cos_sq(t) = cos(π/2 * t^2)

    # Compute Si(1) using scipy's simpson rule
    x_sinc = range(0, 1, length=1000) 
    y_sinc = sinc.(x_sinc./π)
    si_scipy = scipy.integrate.simpson(y_sinc, x_sinc)

    # Compute C(5) using scipy's simpson rule
    x_cos = range(0, 5, length=1001)
    y_cos = cos_sq.(x_cos)
    c_scipy = scipy.integrate.simpson(y_cos, x_cos)

    println("Using scipy's Simpson rule via PyCall: \nSi(1) = $si_scipy \nC(5) = $c_scipy")
catch ModuleNotFoundError
    println("PyCall is not installed. Please install it via Pkg mode (press \"]\" -> add PyCall) and try again.")
end



""" Output:
Using the adaptive Simpson's rule: 
Si(1) = 0.8414709848078965 
C(5) = 0.5636311907377783
Using the "adaptive Trapez rule": 
Si(1) = 0.9460830688712603 
C(5) = 0.5636311867991783
Using scipy's Simpson rule via PyCall: 
Si(1) = 0.9460830703671914 
C(5) = 0.5636312021714622
"""
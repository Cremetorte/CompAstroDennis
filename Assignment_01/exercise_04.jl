using Plots
using LaTeXStrings


function π_polygon(n)
    if n == 2
        return 2 * sqrt(2)
    end
    
    term_1 = sqrt(1 - (π_polygon(n-1)/2^(n-1))^2)   # There is a term ~π^2/2^(2n) subtracted from 1 
                                                    # -> after 15 it: ≈ 10/2^30 ≈ 2^(-27)
    term_2 = sqrt(2*(1- term_1))
    return 2^(n-1) * term_2
end


function π_kahan(n)
    # A_0 = 2√2

    if n == 0
        A_0 = 2*sqrt(2)
        Z_0 = 2*(A_0/2)^2 / (1 + sqrt(1 - A_0^2))
        return Z_0
    end

    A_n = 2^(n-1) * sqrt(4 * π_kahan(n-1))

    term_1 = 2 * (A_n/2^(n+1))^2
    term_2 = 1 + sqrt(1 - (A_n/2^n)^2)

    return term_1/term_2


end



function logarithmic_error(x)
    return log(abs(x-pi))
    
end


n_array = 2:25
pi_approx = π_polygon.(n_array)
# pi_approx_kahan = π_kahan.(n_array .- 2)

n_array = n_array * n_array
display(n_array)


plot(n_array, logarithmic_error.(pi_approx),
    marker=(:xcross, 4), label="Approximation using Polygons")
plot!(n_array, logarithmic_error.(pi_approx_kahan),
    marker=(:xcross, 4), label="Approximation using Kahan algorithm")
xlabel!("Iterations")
ylabel!(L"\log|\pi_{approx} - \pi|")
title!("Logarithmic deviation from π")
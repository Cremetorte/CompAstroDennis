using Plots
gr()
using LaTeXStrings
using BenchmarkTools  # added just now

relevant_ks = [10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9]


function run_sum_forwards()
    #define 1 in 32 and 64 bit floats
    one_32::Float32 = 1
    one_64::Float64 = 1

    sum_32::Float32 = 0
    sum_64::Float64 = 0

    # Prepare results
    results_32 = Float32[]
    results_64 = Float64[]

    # the k-s we want to plot 
    relevant_ks = [10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9]


    for k in 1:10^8
        # calculate sums
        sum_32 += one_32 / k^2
        sum_64 += one_64 / k^2

        if k ∈ relevant_ks 
            # push to back of results
            push!(results_32, sum_32)
            push!(results_64, sum_64)
        end

    end
    return results_32, results_64
end

function run_sum_backwards()
    #define 1 in 32 and 64 bit floats
    one_32::Float32 = 1
    one_64::Float64 = 1

    sum_32::Float32 = 0
    sum_64::Float64 = 0

    # Prepare result 
    results_32 = Float32[]
    results_64 = Float64[]

    
    # the k-s we want to plot
    relevant_ks = [10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9]

    for k in reverse(1:10^8)
        # calculate sums
        sum_32 += one_32 / k^2
        sum_64 += one_64 / k^2

        if k ∈ relevant_ks 
            # push to back of results
            push!(results_32, sum_32)
            push!(results_64, sum_64)
        end

    end
    return results_32, results_64
end


function rel_err(x)
    # π::Float64 = 3.141592653589793115997963468544185161590576171875  # Not neccessary; pi is integrated in julia

    x_64::Float64 = convert(Float64, x) # make sure to convert x to Float64 to have the most precise comparison possible

    true_limit = π^2/6
    # return relative error
    return abs(x_64-true_limit)/true_limit
end

# # Uncomment if forward sum is wanted
# results_32, results_64 = run_sum_forwards()
# display(results_32)

# Uncomment if backward sum is wanted
results_32, results_64 = run_sum_forwards()
display(results_32)

@time res_32_err = rel_err.(results_32)
@time res_64_err = rel_err.(results_64)





plot(relevant_ks, log.(res_32_err), 
    label="Logarithmic rel. error using 32-bit Floats",
    marker=(:xcross, 5))
plot!(relevant_ks, log.(res_64_err), 
    label="Logarithmic rel. error using 64-bit Floats",
    marker=(:xcross, 5))
xlabel!("Iteration steps, backward sum")    # Change backward ↔ forward
ylabel!(L"\log(|s_k - s_\infty|/s_\infty)")
title!("Relative error vs iteration steps")
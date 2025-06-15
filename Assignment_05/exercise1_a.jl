include("exercise_1_shared_code.jl")
using Plots


# -------------------------------------------------- to find dimensionless mass
function trapez_int(x,y)
    res = 0
    for i in 1:length(x)-1
        res += (y[i] + y[i+1]) * (x[i+1] - x[i]) / 2
    end    
    return res
end



n_from_gamma(γ) = 1/(γ-1)

γs = [3, 5/3, 7/5]
ns = n_from_gamma.(γs)

for i in eachindex(γs)
    n = ns[i]
    gamma = γs[i]
    println("Solving Lane-Emden with n=$(n)...")
    ξ, w = solve_Lane_Emden(1e-2, n)

    # total mass
    mass = trapez_int(ξ, w)

    println("Solved! \nRadius: $(ξ[end])\ntotal mass: $(mass)\n")


    plot(ξ, w, xlabel="ξ", ylabel="w", title="Lane-Emden Solution (n=$(n))", legend=false)
    xlabel!("ξ")
    ylabel!("w")
    title!("Lane-Emden Solution for n = $n")
    savefig("Plots/numeric_gamma=$(gamma).png")
end
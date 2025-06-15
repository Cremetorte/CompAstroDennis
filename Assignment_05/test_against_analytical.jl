include("exercise_1_shared_code.jl")

using Plots
gr()


w_0_analytic(ξ) = 1 - 1/6 * ξ^2
w_1_analytic(ξ) = sinc(ξ/π)
w_5_analytic(ξ) = 1/sqrt(1+ξ^2/3)

n_arr = [0,1,5]

analytic_funcs = (w_0_analytic, w_1_analytic, w_5_analytic)


for i in 1:3
    println("Solving Lane-Emden with n=$(n_arr[i])...")
    ξ, w = @time solve_Lane_Emden(1e-2, n_arr[i])
    println("Solved!\n")

    plot(ξ, w, label="n = $i", dpi=600)
    plot!(ξ, analytic_funcs[i].(ξ), label="analytic n = $i", linestyle=:dash)
    xlabel!("ξ")
    ylabel!("w")
    title!("Lane-Emden Solution for n = $i")
    savefig("Plots/analytic_lane_emden_n=$(n_arr[i]).png")
end


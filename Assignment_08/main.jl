include("methods.jl")
using Plots
gr()


# ---------------------
# Simulation Parameters
# ---------------------
xmin = -1
xmax = 1
gridpoints = 400
σ = 0.8
a = 1
end_t = 4   

# -----------
# Step Widths
# -----------
Δx = (xmax - xmin)/gridpoints
Δt = σ * Δx / a

# ----------------------
# Lattice_initialization
# ----------------------
u_0 = [(i*Δx > -1/3 + gridpoints*Δx/2) && (i*Δx < 1/3 + gridpoints*Δx/2) ? 1 : 0 for i in 1:gridpoints]
# "readability counts"

# -----------------------------------
# Analytical solution: shifted by a*t
# -----------------------------------
u_analytical = circshift(u_0, a*end_t)
# a*end_t can only be an integer because else the rectangle peak doesnt fit onto the grid perfectly
# But a=1, end_t=4 is the only thing needed


# ----------------------------------
# Helper function for scaling x-axis
# ----------------------------------
x_axis(grid) = [xmin + i*Δx for i in 1:length(grid)]



# --------------------------------
# solve linear advection ∀ methods
# --------------------------------

methods = Dict("Upwind" => upwind, "Centered Differences" => centered_diff, "Lax-Wendroff" => lax_wendroff)


for method_label in keys(methods)
    # Run solving algorithm and save data in u_sol
    println("Solving Linear Advection using the $method_label-Method")
    u_sol = solve_lin_adv((methods[method_label]), σ, u_0, Δt, end_t)

    # run solving algorithm after it compiled when it was called above to get pure run time
    println("Benchmark results on compiled functions:")
    @time solve_lin_adv((methods[method_label]), σ, u_0, Δt, end_t)

    # plot results
    plot(x_axis(u_sol), u_sol, label="Numeric Solution using the $method_label-Method", grid=true)
    plot!(x_axis(u_analytical), u_analytical, label="Analytical Solution")
    xaxis!("x")
    yaxis!("u(x,t=4)")
    title!("Numeric Solution using the $method_label-Method")
    savefig("Plots/$(method_label)_plot.png")

    println("\n")
end
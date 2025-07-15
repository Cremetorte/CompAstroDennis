using Plots
try
    using DataFrames
    using CSV
catch
    println("Could not import all necessary Packages. Do you wish to install DataFrames and CSV right now?")
    println("type y or n: ")
    usrinput = readline()
    if usrinput == "y"
        using Pkg
        Pkg.add("DataFrames")
        Pkg.add("CSV")
    else
        error("Packages could not be installed")
    end
end
using DataFrames
using CSV
include("lib.jl")

# ---------------------
# Simulation Parameters
# ---------------------
xmin = 0
xmax = 100
gridpoints = 500
t_end = 100
c_s = 1

# -----------
# Step Widths
# -----------
Δx = (xmax - xmin)/gridpoints

# --------------
# Initial Q1, Q2
# --------------
x(index) = (index-1) * Δx
Q_1 = [1 + 0.3 * exp(-(x(i)-50)^2/10) for i in 1:gridpoints]
Q_2 = zeros(gridpoints)

plot(x.(eachindex(Q_1)),Q_1)
savefig("Plots/Q1_initial.png")



output = Matrix{Float64}(undef, 0, 1 + gridpoints)
t = 0
while t < t_end
    global t, output
    t += calculate_step!(Q_1, Q_2, Δx)
    row = [t; Q_1]
    output = vcat(output, row')
end

plot(x.(eachindex(Q_1)),Q_1)
savefig("Plots/Q1_final.png")


# -----------
# File Output
# -----------
df = DataFrame(output, ["time", ["ρ_$i" for i in 1:gridpoints]...])
CSV.write("output.csv", df)

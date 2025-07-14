using Plots
try
    using DataFrames
    using CSV
catch
    println("Could not import all necessary Packages. Do you wish to install DataFrames and CSV?")
    println("type y or n: ")
    usrinput = readline()
    if usrinput == "y"
        using Pkg
        add("DataFrames")
        add("CSV")
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


output = []
temp_output = []
t = 0
while t <= t_end
    global t
    t += calculate_step!(Q_1, Q_2, Δx)
    push!(output, t, Q_1...)
end


display(output)

# -----------
# File Output
# -----------
df = DataFrame(output)
CSV.write("output.csv", df)

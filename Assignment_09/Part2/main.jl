using Plots
using Base.Threads

try
    using DataFrames
    using CSV
    using Dates
catch
    println("Could not import all necessary Packages. Do you wish to install Dates, DataFrames and CSV right now?")
    println("type y or n: ")
    usrinput = readline()
    if usrinput == "y"
        using Pkg
        Pkg.add("DataFrames")
        Pkg.add("CSV")
        Pkg.add("Dates")
    else
        error("Packages could not be installed")
    end
end
using DataFrames
using CSV
using Dates
include("lib.jl")

# ---------------------
# Simulation Parameters
# ---------------------
xmin = 0
xmax = 100
gridpoints = 500
t_end = 40
γ = 1.4

# -----------
# Step Widths
# -----------
Δx = (xmax - xmin)/gridpoints

# --------------
# Initial Q1, Q2
# --------------
x(index) = (index-1) * Δx
Q_1 = [i*Δx <=50 ? 2.0 : 1.0 for i in 1:gridpoints]
Q_2 = zeros(Float64, gridpoints)
Q_3 = ones(Float64, gridpoints) .* Q_1

println("Type of Q_3: ", typeof(Q_3))

plot(x.(eachindex(Q_1)),Q_1)
savefig("Plots/Q1_initial.png")



println("Starting simulation with $(Threads.nthreads()) threads...")

start_time = now()

output_rho = Matrix{Float64}(undef, 0, 1 + gridpoints)
output_momentum = Matrix{Float64}(undef, 0, 1 + gridpoints)
output_internal_energy = Matrix{Float64}(undef, 0, 1 + gridpoints)
t = 0
while t < t_end
    global t, output_rho, output_momentum, output_internal_energy
    t += calculate_step!(Q_1, Q_2, Q_3, Δx, γ)

    row_rho = [t; Q_1]
    output_rho = vcat(output_rho, row_rho')

    row_momentum = [t; Q_2]
    output_momentum = vcat(output_momentum, row_momentum')
    
    row_internal_energy = [t; Q_3]
    output_internal_energy = vcat(output_internal_energy, row_internal_energy')
end

end_time = now()
println("Simulation block elapsed time: ", end_time - start_time)

plot(x.(eachindex(Q_1)),Q_1)
savefig("Plots/Q1_final.png")


# -----------
# File output
# -----------
df = DataFrame(output_rho, ["time", ["ρ_$i" for i in 1:gridpoints]...])
CSV.write("Output/output_rho.csv", df)

df = DataFrame(output_momentum, ["time", ["ρ_$i" for i in 1:gridpoints]...])
CSV.write("Output/output_momentum.csv", df)

df = DataFrame(output_internal_energy, ["time", ["ρ_$i" for i in 1:gridpoints]...])
CSV.write("Output/output_internal_energy.csv", df)

println("Simulation complete! Data saved to CSV files.")

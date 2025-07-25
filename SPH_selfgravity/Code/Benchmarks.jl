using BenchmarkTools
include("DataHandling.jl")
include("Physics.jl")
using Profile
using PProf


particles, N = init_particles("../initial_distribution/random_distribution.dat")


n = 1.0                 # polytropic index
γ = 1.0 + 1/n             # polytropic exponent
K = 0.1                 # polytropic constant
ν = 0.1                 # viscosity coefficient
λ = 2.01203286081606    # linear gravity
h = 0.2                 # smoothing length
t_end = 10.0            # end time
Δt = 0.01               # time step

println("Benchmarking kernel function...")
display(@benchmark div_kernel(part_i, part_j, $h) setup=(part_i = particles.pos[rand(1:N)]; part_j = particles.pos[rand(1:N)]; h = $h))
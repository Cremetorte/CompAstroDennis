include("Code/DataHandling.jl")
include("Code/Physics.jl")
using Plots
using Dates
using Base.Threads
using StaticArrays


initial_file = "initial_distribution/random_distribution.dat"


# Initialize Particles using Data
println("Reading Data...")
particles, N = init_particles(initial_file)
println("Data read. N=$N Particles initialized.")

println("plotting initial positions...")

x = [p[1] for p in particles.pos]
y = [p[2] for p in particles.pos]
z = [p[3] for p in particles.pos]
scatter(x,y,z, title="Initial Particle Positions", xlabel="X", ylabel="Y", zlabel="Z", legend=false)
savefig("Plots/initial_particle_positions.png")


# ---------------------
# Simulation parameters
# ---------------------
n = 1.0                 # polytropic index
γ = 1.0 + 1/n             # polytropic exponent
K = 0.1                 # polytropic constant
ν = 0.1                 # viscosity coefficient
λ = 2.01203286081606    # linear gravity
h = 0.2                 # smoothing length
t_end = 20.0            # end time
Δt = 0.01               # time step
println("Simulation parameters set: n=$n, γ=$γ, K=$K, ν=$ν, h=$h, t_end=$t_end, Δt=$Δt")

t = 0.0         # initial time
nr_timesteps = Int(floor(t_end / Δt) + 1)


# ---------------
# Setup animation
# ---------------
frames_position = Vector{Plots.Plot}(undef, 0)
frames_density = Vector{Plots.Plot}(undef, 0)

# ----------------
# Start Simulation
# ----------------
println("Starting simulation with $(Threads.nthreads()) threads...\n")
start_time = now()
for step in 1:nr_timesteps
    # Update time
    global t += Δt
    
    # Calculate calculate_properties
    calculate_properties!(particles, K, γ, h, λ, ν)
    
    # Calculate accelerations
    calculate_accelerations_symmetric!(particles, h, λ, ν)
    
    # integrate a, vel
    integrate_vel_acc!(particles, Δt)
    
    if step % 10 == 0
        print("Calculating time $(rpad(round(t, digits=4), 4, '0')) ≡ step $(lpad(step, 4, '0'))/$nr_timesteps...    \r")
        global x = [p[1] for p in particles.pos]
        global y = [p[2] for p in particles.pos]
        global z = [p[3] for p in particles.pos]
        global r = [norm(p) for p in particles.pos]
        # Store frame for animation
        push!(frames_position, scatter(x, y, z, title="Particle Positions at time $(rpad(round(t, digits=4), 4, '0'))", xlabel="X", ylabel="Y", zlabel="Z", legend=false,))
        push!(frames_density, scatter(r, particles.ρ, title="Particle Densities at time $(rpad(round(t, digits=4), 4, '0'))", xlabel="Calculated Density", ylabel="rho", legend=true))
    end
end
end_time = now()
println("Simulation completed in $(end_time - start_time).")


# -----------------------------
# calculate denisty between 0,2R
# -----------------------------
r = range(0, stop=2*0.75, length=100)
densities = zeros(length(r))
for i in 1:length(r)
    densities[i] = 0
    for j in 1:N
        densities[i] += particles.mass[j] * kernel(SVector{3, Float64}(r[i],0,0), particles.pos[j], h)
    end
end

density_analytical(radius) = (λ/(2*K*(1+n)) * (0.75^2 - radius^2))^n
densities_analytical = density_analytical.(r)
# Plotting density profile
println("Plotting density profile...")
plot(r, densities, title="Density Profile", xlabel="Distance (r)", ylabel="Density", legend=true)
plot!(r, densities_analytical, label="Analytical Density", linestyle=:dash)
savefig("Plots/density_profile.png")




# ---------------------------------
# Plotting final particle positions
# ---------------------------------
println("Saving final particle positions...")
x = [p[1] for p in particles.pos]
y = [p[2] for p in particles.pos]
z = [p[3] for p in particles.pos]
scatter(x,y,z, title="Particle Positions", xlabel="X", ylabel="Y", zlabel="Z")
savefig("Plots/final_particle_positions.png")


# ------------
# Generate GIF
# ------------
println("Generating GIFs...")
anim_pos = @animate for i in 1:length(frames_position)
    plot!(frames_position[i], xlims=(-1, 1), ylims=(-1, 1), zlims=(-1, 1))
end
gif(anim_pos, "Animations/Toy_Star_pos.gif", fps=30)

anim_density = @animate for i in 1:length(frames_position)
    plot!(frames_density[i], xlims=(0, 0.8), ylims=(0, 3))
    plot!(r, densities_analytical, label="Analytical Density", linestyle=:dash)
end
gif(anim_density, "Animations/Toy_Star_density.gif", fps=30)
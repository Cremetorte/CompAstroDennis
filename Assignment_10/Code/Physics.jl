using StaticArrays
using LinearAlgebra
using Base.Threads
include("DataHandling.jl")


# Actually slower, dont use
function kernel_branchless(pos_i::SVector{3,Float64}, pos_j::SVector{3,Float64}, h=0.2::Float64)
    l = norm(pos_j .- pos_i)/h
    
    #branchless calculation
    w = (l < 0.5) * (6*l^3 - 6*l^2 +1) + (l >= 0.5 && l <= 1) * 2 * (1 - l)^3 + 0.0
    return 8/π / h^3 * w
end

# faster version
function kernel_ifelse(pos_i::SVector{3,Float64}, pos_j::SVector{3,Float64}, h=0.2::Float64)
    l = norm(pos_j .- pos_i)/h

    divisor = π*h^3
    
    if l < 0.5
        w = (6*l^3 - 6*l^2 +1)
    elseif l >= 0.5 && l <= 1
        w = 2 * (1 - l)^3
    else
        w = 0
    end

    return 8/divisor * w
end


# Fastest implementation
function kernel(pos_i::SVector{3,Float64}, pos_j::SVector{3,Float64}, h=0.2::Float64)
    l = norm(pos_j - pos_i)/h
    
    if l > 1
        return 0.0
    elseif l < 0.5
        w = (6*l^3 - 6*l^2 +1)
    else
        w = 2 * (1 - l)^3
    end

    return 8/π / h^3 * w
end

function div_kernel(pos_i::SVector{3,Float64}, pos_j::SVector{3,Float64}, h=0.2::Float64)
    q = pos_i - pos_j
    
    l = norm(q)/h
    


    if l > 1 || l == 0
        # Avoid division by zero and outside support radius
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif l < 0.5
        w = 3*l^2 - 2*l
    else
        w = -(1 - l)^2
    end
    scalar = 6 * 8 / (l * π * h^4) * w
    return q * scalar
end


function pressure(density, K::Float64, γ::Float64)
    return K * density^γ
end


function soundspeed(density, K::Float64, γ::Float64)
    return sqrt(K * γ * density^(γ-1))
end

function acceleration_grav_damp(position, velocity, λ::Float64, ν::Float64)
    return (-λ .* position) - (ν .* velocity)
end


function calculate_properties!(particles, K=0.1, γ=2, h=0.2, λ=2.01203286081606, ν=0.1)
    # Äußere For-loop in Threads
    N = length(particles.ρ)
    @threads for i in 1:N
        # lokale temporäre Variable für Thread-Safety
        rho_i = 0
        # Vektorisiert mit @simd, theoretisch mit cuda machbar?
        @inbounds for j in 1:N
            rho_i += kernel(particles.pos[i], particles.pos[j], h) * particles.mass[j]
        end

        # Speichern der berechneten Werte
        particles.ρ[i] = rho_i
        particles.pressure[i] = particles.ρ[i], K, γ
        particles.soundspeed[i] = soundspeed(particles.ρ[i], K, γ)
    end
end

function calculate_accelerations!(particles, h=0.2, λ=2.01203286081606, ν=0.1)
    # Berechnung der Beschleunigungen
    N = length(particles.ρ)
    @threads for i in 1:N
        particles.acc[i] = acceleration_grav_damp(particles.pos[i], particles.vel[i], λ, ν)

        pressure_dens_term_i = particles.pressure[i]/particles.ρ[i]^2

        @inbounds @simd for j in 1:N
            
            particles.acc[i] = particles.acc[i] - div_kernel(particles.pos[i], particles.pos[j], h) .* 
            (particles.mass[j] * (pressure_dens_term_i + 
            particles.pressure[j]/particles.ρ[j]^2) )
            
        end
    end
end


function calculate_accelerations_symmetric!(particles, h=0.2, λ=2.01203286081606, ν=0.1)
    N = length(particles.ρ)
    
    # Acc before other loop to get all particles
    @threads for i in 1:N
        @inbounds particles.acc[i] = acceleration_grav_damp(particles.pos[i], particles.vel[i], λ, ν)
    end
    
    # Calculate symmtric forces
    @threads for i in 1:N-1
        for j in i+1:N
            # i!=j but all i,j get caught
            @inbounds begin
                force_ij = div_kernel(particles.pos[i], particles.pos[j], h) * 
                          particles.mass[j] * (particles.pressure[i]/particles.ρ[i]^2 + 
                                              particles.pressure[j]/particles.ρ[j]^2)
                
                particles.acc[i] = particles.acc[i] - force_ij
                particles.acc[j] = particles.acc[j] + force_ij 
            end
        end
    end
end

function integrate_vel_acc!(particles, Δt=1e-2)
    N = length(particles.ρ)
    @threads for i in 1:N
        @inbounds begin
            # Update velocity
            particles.vel[i] = particles.vel[i] + particles.acc[i] .* Δt

            # Update position
            particles.pos[i] = particles.pos[i] + particles.vel[i] .* Δt
        end
    end
end

function integrate_leapfrog!(particles, Δt=1e-2, K=0.1, γ=2, h=0.2, λ=2.01203286081606, ν=0.1)
    N = length(particles.ρ)
    # calculate t_1/2 for all particles
    @threads for i in 1:N
        @inbounds begin
            particles.vel[i] = particles.vel[i] + particles.acc[i] * 0.5 * Δt
            particles.pos[i] = particles.pos[i] + particles.vel[i] * 0.5 * Δt
        end
    end

    calculate_properties!(particles, K, γ, h, λ, ν)
    calculate_accelerations!(particles, h, λ, ν)
    #calculate new accelerations, pos, vel get_starting_values
    @threads for i in 1:N
        @inbounds begin
            particles.pos[i] = particles.pos[i] + particles.vel[i] * 0.5 * Δt
            particles.vel[i] = particles.vel[i] + particles.acc[i] * 0.5 * Δt
        end
    end
    
end





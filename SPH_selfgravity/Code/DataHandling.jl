using StaticArrays
using CSV
using DataFrames

# Struct definition, irgendwas ist komisch
struct Particles
    pos::Vector{SVector{3, Float64}}
    vel::Vector{SVector{3, Float64}}
    acc::Vector{SVector{3, Float64}}
    mass::Vector{Float64}
    ρ::Vector{Float64}
    pressure::Vector{Float64}
    soundspeed::Vector{Float64}
end

function init_particles(file)
    # Read file
    df = CSV.read(file, DataFrame)
    
    #extract Data
    N = size(df)[1]
    r = [SVector{3, Float64}(df[i, 1], df[i, 2], df[i, 3]) for i in 1:N]
    v = [SVector{3, Float64}(df[i, 4], df[i, 5], df[i, 6]) for i in 1:N]
    a = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:N]
    m = [df[i, 7] for i in 1:N]
    ρ = [0.0 for _ in 1:N]
    pressure = [0.0 for _ in 1:N]
    soundspeed = [0.0 for _ in 1:N]

    return Particles(r, v, a, m, ρ, pressure, soundspeed), N    
end


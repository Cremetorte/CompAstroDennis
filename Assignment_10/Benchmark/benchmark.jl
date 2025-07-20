using BenchmarkTools
using StaticArrays
include("../lib/Physics.jl")

N = 100000
println("Benchmarking Kernel with $N Particles")
println("Branchless:")
@btime kernel_branchless.(positions_i, positions_j, h) setup=(
    positions_i = [SVector{3}(rand(3)) for _ in 1:N];
    positions_j = [SVector{3}(rand(3)) for _ in 1:N];
    h = rand()
)

# Benchmark für kernel_ifelse
println("\nIf/else:")
@btime kernel_ifelse.(positions_i, positions_j, h) setup=(
    positions_i = [SVector{3}(rand(3)) for _ in 1:N];
    positions_j = [SVector{3}(rand(3)) for _ in 1:N];
    h = rand(1) * 0.5
)

# Benchmark für kernel:kernel_optimized
println("\nOptimized If/else")
@btime kernel_optimized.(positions_i, positions_j, h) setup=(
    positions_i = [SVector{3}(rand(3)) for _ in 1:N];
    positions_j = [SVector{3}(rand(3)) for _ in 1:N];
    h = rand(1) * 0.5
)

include("Physics.jl")
using StaticArrays

h = 0.2


function kernel_scalar(r_tilde)
    global h
    r_0 = SVector{3, Float64}(0,0,0)
    r_tilde_vec = SVector{3, Float64}(r_tilde,0,0)
    return kernel(r_0, r_tilde_vec, h)
end

using Plots

r_vals = range(0, 2h, length=100)
k_vals = [kernel_scalar(r) for r in r_vals]

area = sum(k_vals) * (2h / 100)
println("Integral of kernel: ", area*2π)
println("Kernel at r=0: $(kernel_scalar(0))")
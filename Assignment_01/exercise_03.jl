function sqrt_sum_infinity()
    a::Float32 = 1e30
    b::Float32 = 1.0

    res::Float32 = sqrt(a*a + b*b)  # res = inf because a*a = 1e60 > Float32 limit (~10^38)
    return res
end


function sqrt_sum_rescaled()
    a::Float32 = 1e30
    b::Float32 = 1.0

    res::Float32 = sqrt(a) * sqrt(a + b^2/a)    # sqrt(a) = 1e15, sqrt(a+b^2/a) = sqrt(10^30+10^-30) -> inside Float32-Limit
    return res
end

function sqrt_sum_Float64()
    # Float64 - Version
    a::Float64 = 1e30
    b::Float64 = 1.0

    res::Float64 =  sqrt(a^2 + b^2)     # a^2 will swallow b^2 -> a^2 + b^2 = a^2
    return res
end

function sqrt_sum_rescaled_Float64()
    # Float64 - Version
    a::Float64 = 1e30
    b::Float64 = 1.0

    res::Float64 =  sqrt(a)*sqrt(a + b^2/a)     # a^2 will swallow b^2 -> a^2 + b^2 = a^2
    return res
end

standard = sqrt_sum_infinity()
rescaled = sqrt_sum_rescaled()
sum_64 = sqrt_sum_Float64()
sum_64_rescaled = sqrt_sum_rescaled_Float64()

println("Standard result: $standard")
println("Rescaled result: $rescaled") # a little bit off due to Floating point rounding errors
println("Float64 result: $sum_64")
println("Rescaled Float64 result: $sum_64_rescaled")

# Finds machine precision ϵ depending on var_type (eg. Float64, Float32, Float16,...)
function find_epsilon(var_type)
    one = convert(var_type, 1.0)    # Convert 1 to var_type

    epsilon = one # start value
    while true
        if one + epsilon == one
            break
        else
            epsilon /=2
        end
    end
    return epsilon*2    # *2 because epsilon itself is the value where 1 + ϵ = 1
end

function machine_prec(var_type)
    one = convert(var_type, 1.0)
    return nextfloat(one) - one     # nextfloat: smallest float bigger than 1 ⟺ machine precision
end


# calculate ϵ-s
ϵ_32 = find_epsilon(Float32)
true_ϵ_32 = machine_prec(Float32)
ϵ_64 = find_epsilon(Float64)
true_ϵ_64 = machine_prec(Float64)

println("Float32: $ϵ_32, true value: $true_ϵ_32")
println("Float64: $ϵ_64, true value: $true_ϵ_64")
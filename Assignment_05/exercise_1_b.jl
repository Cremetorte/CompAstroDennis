include("exercise_1_shared_code.jl")


# data
m_sun = 1.989e33 # g
R_sun = 6.96e10 # cm
mean_density_sun = 1.41 # g/cm^3
mean_molecular_weight_sun = 0.62 # g/mol
G = 6.67430e-8 # cm^3/(g*s^2) # gravitational constant
k_B = 1.380649e-16 # erg/K # Boltzmann constant
m_H = 1.6735575e-24 # g # mass of hydrogen

# calculate the sun
n=3

ξ, w, z = solve_Lane_Emden_with_derivative(1e-2, n)

# at the root
ξ_1 = ξ[end]
z_1 = z[end] # = dw/dξ at the root

ρ_c_over_ρ_bar = ξ_1^3 / (-3 * ξ_1^2 * z_1)

critical_density = ρ_c_over_ρ_bar * mean_density_sun    # multipied by mean density of the sun

α = R_sun / ξ_1

K = 4π * G *α^2 / (n + 1) /critical_density^(1/n - 1)

p_c = K * critical_density^(1 + 1/n)

T_c = p_c * mean_molecular_weight_sun * m_H / (k_B * critical_density)

println("Critical density: $(critical_density) g/cm^3") #
println("α: $(α) cm")
println("K: $(K) cm^2/s^2")
println("p_c: $(p_c) erg/cm^3")
println("T_c: $(T_c) K")
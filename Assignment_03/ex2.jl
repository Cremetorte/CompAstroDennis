using Roots


# M = 5.974× 1024kg, m = 7.348× 1022kg, R = 3.844× 108m, ω = 2.662× 10−6s−1
M = 5.974e24
m = 7.348e22
R = 3.844e8
ω = 2.662e-6
G = 6.674e-11

f(r) = G*M/r^2 - G*m/(R - r)^2 - ω^2*r

# r finden sodass f(r) = 0
r_sol = find_zero(f, (0, R))

println("Lagrange Point 1: $r_sol")
# Lagrange Point 1: 3.2604e8 m
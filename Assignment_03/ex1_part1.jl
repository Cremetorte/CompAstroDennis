using Plots

# Variable um Iterationsschritte zu zählen
nr_it = 0

# Fixpunkt-Schema zur Lösung der Keplergleichung
function fixed_point_E(mean_anomaly, eccentricity, epsilon=1e-9, max_steps=1e10)
    global nr_it

    # Startwert finden
    if eccentricity <= 0.8
        E = mean_anomaly
    else
        E = pi
    end

    # Iterationen durchführen
    for i in 1:max_steps
        nr_it += 1
        E_new = mean_anomaly + eccentricity * sin(E)
        # convergence criterion
        if abs(E_new - E) < epsilon
            return E_new
        end
        #set new E
        E = E_new
    end
    println("Fixed-Point-Scheme did not converge in $max_steps steps")
end

# Newton-Raphson Methode für Keplergleichung
function newton_E(mean_anomaly, eccentricity, epsilon=1e-9, max_steps=1e10)
    global nr_it

    # Startwert finden
    if eccentricity <= 0.8
        E = mean_anomaly
    else
        E = pi
    end

    # Functions g, g'
    g(x) = x - eccentricity * sin(x) - mean_anomaly
    g_prime(x) = 1- eccentricity * cos(x)

    # Iteratioen durchführen
    for i in 1:max_steps
        nr_it += 1
        E_new = E - g(E)/g_prime(E)

        # convergence criterion
        if abs(E_new - E) < epsilon
            return E_new
        end

        # set E to the new E
        E = E_new
    end
    println("Newton-Raphson-method did not converge in $max_steps steps")
end

# true anomaly f
function true_anomaly(eccentricity, eccentric_anomaly)
    root_stuff = sqrt((1+eccentricity)/(1-eccentricity))
    return 2 * atan(root_stuff * tan(eccentric_anomaly/2))
end

# Abstand r zwischen Körper und Fokus
function distance(semi_major_axis, eccentricity, true_anomaly)
    return semi_major_axis * (1-eccentricity^2) / (1 + eccentricity*cos(true_anomaly))
end

# Position des Körpers
function position(M, e, a, solver=newton_E)
    E = solver(M, e)
    f = true_anomaly(e, E)
    r = distance(a, e, f)

    return [r*cos(f), r*sin(f)]
end

# lineare Range als Anomalie
mean_anomalies = LinRange(0, 2π, 256)

#Konstanten
e_mercury = 0.205
a_mercury = 0.39

e_halley = 0.967
a_halley = 17.8



positions_mercury = position.(mean_anomalies, e_mercury, a_mercury)
println("Newton-Raphson-method iterations for Mercury: $nr_it")
# Gives 962 Iterations

nr_it = 0  # Reset iteration count for Halley's Comet
positions_mercury = position.(mean_anomalies, e_mercury, a_mercury, fixed_point_E)
println("Fixed-Point-Scheme iterations for Mercury: $nr_it")
# Gives 2640 Iterations

nr_it = 0  # Reset iteration count for Halley's Comet
positions_halley = position.(mean_anomalies, e_halley, a_halley)
println("Newton-Raphson-method iterations for Halley's Comet: $nr_it")
# 1266 Iterations

nr_it = 0  # Reset iteration count for Halley's Comet
positions_halley = position.(mean_anomalies, e_halley, a_halley, fixed_point_E)
println("Fixed-Point-Scheme iterations for Halley's Comet: $nr_it")
# 33102 Iterations

positions_mercury = reduce(hcat, positions_mercury)'
positions_halley = reduce(hcat, positions_halley)'



# Plots
plot(positions_mercury[:, 1], positions_mercury[:, 2], label="Mercury", color=:blue, aspect_ratio=:equal)
xlabel!("x-position (AU)")
ylabel!("y-position (AU)")
title!("Orbits of Mercury")
savefig("orbit_mercury.png")


plot(positions_halley[:, 1], positions_halley[:, 2], label="Halley's Comet", color=:red, aspect_ratio=:equal)
xlabel!("x-position (AU)")
ylabel!("y-position (AU)")
title!("Orbits of Halley's Comet")
savefig("orbit_halley.png")
# Display the plot
# display(plot())
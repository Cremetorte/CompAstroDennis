using Plots
using Dates


function newton_E(mean_anomaly, eccentricity, epsilon=1e-9, max_steps=1e10)
    if eccentricity <= 0.8
        E = mean_anomaly
    else
        E = pi
    end

    # Functions g, g'
    g(x) = x - eccentricity * sin(x)
    g_prime(x) = 1- eccentricity * cos(x)

    for i in 1:max_steps
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



function newton_E(mean_anomaly, eccentricity, epsilon=1e-9, max_steps=1e10)
    global nr_it
    if eccentricity <= 0.8
        E = mean_anomaly
    else
        E = pi
    end

    # Functions g, g'
    g(x) = x - eccentricity * sin(x) - mean_anomaly
    g_prime(x) = 1- eccentricity * cos(x)

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


function true_anomaly(eccentricity, eccentric_anomaly)
    root_stuff = sqrt((1+eccentricity)/(1-eccentricity))
    return 2 * atan(root_stuff * tan(eccentric_anomaly/2))
end


function distance(semi_major_axis, eccentricity, true_anomaly)
    return semi_major_axis * (1-eccentricity^2) / (1 + eccentricity*cos(true_anomaly))
end


function position(M, e, a, solver=newton_E)
    E = solver(M, e)
    f = true_anomaly(e, E)
    r = distance(a, e, f)

    return [r*cos(f), r*sin(f)]
end


function distance(r_1, r_2)
    return sqrt((r_1[1] - r_2[1])^2 + (r_1[2] - r_2[2])^2) 
end


# duration
t_1 = Date(1985, 1, 1)
t_2 = Date(2025, 5, 21)
duration_days = Dates.value(t_2 - t_1)

# time vector
t = collect(0:duration_days)
t = t ./ 365.25  # convert to years


# Mean anomalies as functions
M_0(lambda, phi_0, period) = lambda - phi_0 - 15 * 2π / period
M(t, M0, period) = M0 + 2π * t / period 



# Constants of earth
period_earth = 1
a_earth = 1
λ_earth = 100.46 / 360 * 2π
Φ_0_earth = 102.95 / 360 * 2π
M_0_earth = M_0(λ_earth, Φ_0_earth, period_earth)
e_earth = 0.0167

# earth's mean anomalies
M_earth = M.(t, M_0_earth, period_earth)


# Constants of mars
a_mars = 1.524
period_mars = (a_mars^3/a_earth^3 * period_earth^3)^0.5     # using Kepler 3
λ_mars = 355.45 / 360 * 2π
Φ_0_mars = 336.04 / 360 * 2π
M_0_mars = M_0(λ_mars, Φ_0_mars, period_mars)
e_mars = 0.0935

# mars' mean anomalies
M_mars = M.(t, M_0_mars, period_mars)




@time r_earth = position.(M_earth, e_earth, a_earth)
@time r_mars = position.(M_mars, e_mars, a_mars)

distances = distance.(r_earth, r_mars)

# Plotting
plot(t, distances, label="Distance between Earth and Mars", 
    xlabel="Time (years)", ylabel="Distance (AU)", 
    title="Distance between Earth and Mars over time",
    legend=:topright, grid=true
    )

# Save the plot
savefig("distance_earth_mars.png")
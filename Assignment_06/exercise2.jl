using Plots
gr()

# from assignment 5
function RK4_integrator_step(f, t_n, u_n, Ω, h)
    K_1 = f(u_n[1], u_n[2], u_n[3], u_n[4], Ω, 1 - 1e-3, 1e-3)
    K_2 = f(u_n[1] + h/2*K_1[1], u_n[2] + h/2*K_1[2], u_n[3] + h/2*K_1[3], u_n[4] + h/2*K_1[4], Ω, 1 - 1e-3, 1e-3)
    K_3 = f(u_n[1] + h/2*K_2[1], u_n[2] + h/2*K_2[2], u_n[3] + h/2*K_2[3], u_n[4] + h/2*K_2[4], Ω, 1 - 1e-3, 1e-3)
    K_4 = f(u_n[1] + h*K_3[1], u_n[2] + h*K_3[2], u_n[3] + h*K_3[3], u_n[4] + h*K_3[4], Ω, 1 - 1e-3, 1e-3)

    return u_n .+ h/6 .* (K_1 .+ 2 .* K_2 .+ 2 .* K_3 .+ K_4)
end


function f(x, y, x_dot, y_dot, Ω, μ_1, μ_2)
    r_1 = sqrt((x + μ_2)^2 + y^2)
    r_2 = sqrt((x - μ_1)^2 + y^2)
    x_ddot = 2 * Ω * y_dot  +  Ω^2 * x  -  (μ_1 * (x + μ_2)/r_1^3 + μ_2 * (x - μ_1)/r_2^3) 
    y_ddot = -2 * Ω * x_dot  -  Ω^2 * y  -  (μ_1/r_1^3 + μ_2/r_2^3)*y
    return [x_dot, y_dot, x_ddot, y_ddot]
end

function get_starting_values(C_J, x_0, Ω, μ_1, μ_2)
    r_1 = sqrt((x_0 + μ_2)^2)
    r_2 = sqrt((x_0 - μ_1)^2)
    y_dot = sqrt(Ω^2 * x_0^2 + 2 * (μ_1/r_1 + μ_2/r_2) - C_J)
    return [x_0, 0, 0, y_dot]
end


function solve_RC3BP(C_J, x_0, step_size=1e-3, nr_steps=5e6)
    Ω = 1
    μ_2 = 1e-3
    μ_1 = 1 - μ_2

    u_0 = get_starting_values(C_J, x_0, Ω, μ_1, μ_2)
    u_n = [u_0]

    x_n = [u_0[1]]
    y_n = [u_0[2]]
    x_dot_n = [u_0[3]]
    y_dot_n = [u_0[4]]

    t_0 = 1e-4
    t_n = [t_0]

    for i in 1:nr_steps
        u_n_i = RK4_integrator_step(f, t_n[end], u_n[end], Ω, step_size)
        push!(u_n, u_n_i)

        x_n_i = u_n_i[1]
        y_n_i = u_n_i[2]
        x_dot_n_i = u_n_i[3]
        y_dot_n_i = u_n_i[4]

        push!(x_n, x_n_i)
        push!(y_n, y_n_i)
        push!(x_dot_n, x_dot_n_i)
        push!(y_dot_n, y_dot_n_i)

        t_0 += step_size
        push!(t_n, t_0)
    end
    return (t_n, x_n, y_n, x_dot_n, y_dot_n)
        
end


starting_xs = [0.21, 0.24, 0.26, 0.27, 0.4, 0.5, 0.6, 0.8]

for x_0 in starting_xs
    C_J = 3.03
    println("Solving RC3BP for C_J=$(C_J) and x_0=$(x_0)...")
    t_n, x_n, y_n, x_dot_n, y_dot_n = @time solve_RC3BP(C_J, x_0)
    println("Solved!\n")

    plot(x_n[1:1000:end], y_n[1:1000:end], label="C_J=$(C_J), x_0=$(x_0)", dpi=600, markersize=2)
    xlabel!("x")
    ylabel!("y")
    title!("RC3BP Solution for C_J=$(C_J) and x_0=$(x_0)")
    savefig("Plots/RC3BP_C_J=$(C_J)_x_0=$(x_0).png")

end


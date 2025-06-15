using Plots
gr()


function RK4_integrator_step(f, ξ_n, y_n, n, h)
    K_1 = f(ξ_n, y_n, n)
    K_2 = f(ξ_n + h/2, y_n + h*K_1/2, n)
    K_3 = f(ξ_n + h/2, y_n + h*K_2/2, n)
    K_4 = f(ξ_n + h, y_n + h*K_3, n)

    return y_n .+ h ./ 6 .* (K_1 .+ 2 .* K_2 + 2 .* K_3 .+ K_4)
end


function f(ξ, y, n)
    w, z = y
    return [z, -(2/ξ * z + w^n)]
end


function solve_Lane_Emden(step_size, n_cur, max_steps=1e4)

    y_0 = [1.0, 0.0]
    y_n = [y_0]

    w_n = [y_0[1]]
    z_n = [y_0[2]]

    ξ_0 = 1e-4

    ξ_n = [ξ_0]

    found_root = false
    for i in 1:max_steps

        new_y = RK4_integrator_step(f, ξ_n[end], y_n[end], n_cur, step_size)
        # detect root
        if new_y[1] < 0
            found_root = true
            break
        end

        #add to output
        push!(y_n, new_y)
        push!(ξ_n, ξ_n[end] + step_size)
        push!(w_n, y_n[end][1])
        push!(z_n, y_n[end][2])
    end

    # poor mans for-else
    if found_root == false
        print("Did not find the root after $max_steps steps.")
    end

    return (ξ_n, w_n)
end



# step_size = 1e-4

# n_cur = 0

# y_0 = [1.0, 0.0]

# y_n = [y_0]

# w_n = [y_0[1]]
# z_n = [y_0[2]]

# ξ_0 = 1e-4

# ξ_n = [ξ_0]


# for i in 1:800840
#     push!(y_n, RK4_integrator_step(f, ξ_n[end], y_n[end], n_cur, step_size))
#     push!(ξ_n, ξ_n[end] + step_size)
#     push!(w_n, y_n[end][1])
#     push!(z_n, y_n[end][2])
# end

# display(w_n)

# w_0_analytic(ξ) = 1 - 1/6 * ξ^2
# w_1_analytic(ξ) = sinc(ξ/π)
# w_5_analytic(ξ) = 1/sqrt(1+ξ^2/3)

# n_arr = [0,1,5]

# analytic_funcs = (w_0_analytic, w_1_analytic, w_5_analytic)


# for i in 1:3
#     ξ, w = solve_Lane_Emden(1e-4, n_arr[i])
#     plot(ξ, w, label="n = $i", dpi=600)
#     plot!(ξ, analytic_funcs[i].(ξ), label="analytic n = $i", linestyle=:dash)
#     xlabel!("ξ")
#     ylabel!("w")
#     title!("Lane-Emden Solution for n = $i")
#     savefig("analytic_lane_emden_n=$(n_arr[i]).png")
# end

function upwind(u_n, σ)
    u_n_left = circshift(u_n, 1)
    return u_n .- σ .* (u_n .- u_n_left)
end

function centered_diff(u_n, σ)
    u_n_left = circshift(u_n, 1)
    u_n_right = circshift(u_n, -1)
    return u_n .- σ/2 .* (u_n_right .- u_n_left)
end

function lax_wendroff(u_n, σ)
    u_n_left = circshift(u_n, 1)
    u_n_right = circshift(u_n, -1)
    return (u_n 
            .- σ/2 .* (u_n_right .- u_n_left) 
            .+ σ^2/2 .* (u_n_right .- 2 .* u_n .+ u_n_left)
            )
end


function solve_lin_adv(method, σ, u_0, Δt, t_end)
    nr_timesteps = floor(Int64, t_end/Δt)
    # println("Running $nr_timesteps timesteps")
    
    u_cur = u_0
    for i in 1:nr_timesteps
        print("\rComputing step $i/$nr_timesteps    ")
        u_cur = method(u_cur, σ)
    end
    print("\n")
    
    return u_cur
end
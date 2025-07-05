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
    nr_timesteps = round(Int64, t_end/Δt)
    
    u_cur = u_0
    for i in 1:nr_timesteps
        u_cur = method(u_cur, σ)
    end
    
    return u_cur
end
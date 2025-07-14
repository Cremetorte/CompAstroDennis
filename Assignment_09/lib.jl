using Base.Threads


c_s = 1

function calculate_step!(Q1_array, Q2_array, Δx)
    N = length(Q1_array)
    Q1_new = similar(Q1_array)
    Q2_new = similar(Q2_array)

    v_max = maximum(c_s .+ abs.(Q2_array ./ Q1_array))
    Δt = 0.5 * Δx / v_max

    Q1_half = zeros(N)
    Q2_half = zeros(N)


    # First step of op-splitting
    @threads for i in 1:N
        Q1 = Q1_array[i]
        Q1_left = (i==1) ? Q1_array[1] : Q1_array[i-1]
        Q1_right = (i==N) ? Q1_array[N] : Q1_array[i+1]

        Q2 = Q2_array[i]
        Q2_left = (i==1) ? -Q2_array[1] : Q2_array[i-1]
        Q2_right = (i==N) ? -Q2_array[N] : Q2_array[i+1]


        u_half_right = 0.5 * (Q2 / Q1 + Q2_right / Q1_right)
        u_half_left = 0.5 * (Q2 / Q1 + Q2_left / Q1_left)

        flux1_right = u_half_right > 0 ? Q1 * u_half_right : Q1_right * u_half_right
        flux2_right= u_half_right < 0 ? Q2 * u_half_right : Q2_right * u_half_right

        flux1_left = u_half_left < 0 ? Q1 * u_half_left : Q1_left * u_half_left
        flux2_left = u_half_left > 0 ? Q2 * u_half_left : Q2_left * u_half_left

        Q1_half[i] = Q1 - Δt/Δx * (flux1_right - flux1_left)
        Q2_half[i] = Q2 - Δt/Δx * (flux2_right - flux2_left)
        
    end



    @threads for i in 1:N
        Q1_left = (i==1) ? Q1_half[1] : Q1_half[i-1]
        Q1_right = (i==N) ? Q1_half[N] : Q1_half[i+1]

        Q1_new[i] = Q1_half[i]
        Q2_new[i] = Q2_half[i] - Δt * c_s^2 * (Q1_right - Q1_left) / 2*Δx 
    end

    Q1_array .= Q1_new
    Q2_array .= Q2_new

    return Δt # for timekeeping
    
end
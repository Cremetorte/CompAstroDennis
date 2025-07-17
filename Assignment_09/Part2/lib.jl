using Base.Threads


c_s(ϵ, γ) = sqrt(γ*(γ-1)*ϵ)

p(ρ, γ, ϵ) = (γ-1) * ρ * ϵ

function calculate_step!(Q1_array, Q2_array, Q3_array, Δx, γ)
    N = length(Q1_array)

    Q1_new = similar(Q1_array)
    Q2_new = similar(Q2_array)
    Q3_new = similar(Q3_array)


    v_max = maximum(c_s.(Q3_array ./ Q1_array, γ) .+ abs.(Q2_array ./ Q1_array))
    Δt = 0.5 * Δx / v_max

    Q1_half = zeros(N)
    Q2_half = zeros(N)
    Q3_half = zeros(N)


    # First step of op-splitting
    @threads for i in 1:N
        Q1 = Q1_array[i]
        Q1_left = (i==1) ? Q1_array[1] : Q1_array[i-1]
        Q1_right = (i==N) ? Q1_array[N] : Q1_array[i+1]

        Q2 = Q2_array[i]
        Q2_left = (i==1) ? -Q2_array[1] : Q2_array[i-1]
        Q2_right = (i==N) ? -Q2_array[N] : Q2_array[i+1]

        Q3 = Q3_array[i]
        Q3_left = (i==1) ? Q3_array[1] : Q3_array[i-1]
        Q3_right = (i==N) ? Q3_array[N] : Q3_array[i+1]


        u_half_right = 0.5 * (Q2 / Q1 + Q2_right / Q1_right)
        u_half_left = 0.5 * (Q2 / Q1 + Q2_left / Q1_left)

        flux1_right = u_half_right > 0 ? Q1 * u_half_right : Q1_right * u_half_right
        flux2_right = u_half_right > 0 ? Q2 * u_half_right : Q2_right * u_half_right
        flux3_right = u_half_right > 0 ? Q3 * u_half_right : Q3_right * u_half_right
        flux1_left = u_half_left > 0 ? Q1_left * u_half_left : Q1 * u_half_left
        flux2_left = u_half_left > 0 ? Q2_left * u_half_left : Q2 * u_half_left
        flux3_left = u_half_left > 0 ? Q3_left * u_half_left : Q3 * u_half_left

        Q1_half[i] = Q1 - Δt/Δx * (flux1_right - flux1_left)
        Q2_half[i] = Q2 - Δt/Δx * (flux2_right - flux2_left)
        Q3_half[i] = Q3 - Δt/Δx * (flux3_right - flux3_left)

    end


    # 2nd step of op-splitting
    @threads for i in 1:N
        Q1_left = (i==1) ? Q1_half[1] : Q1_half[i-1]
        Q1_right = (i==N) ? Q1_half[N] : Q1_half[i+1]

        Q2_left = (i==1) ? -Q2_half[1] : Q2_half[i-1]
        Q2_right = (i==N) ? -Q2_half[N] : Q2_half[i+1]

        Q3_left = (i==1) ? Q3_half[1] : Q3_half[i-1]
        Q3_right = (i==N) ? Q3_half[N] : Q3_half[i+1]

        energy_left = Q3_left/Q1_left
        energy_right = Q3_right/Q1_right

        u_left = Q2_left / Q1_left
        u_right = Q2_right / Q1_right

        p_left = p(Q1_left, γ, energy_left)
        p_right = p(Q1_right, γ, energy_right)

        Q1_new[i] = Q1_half[i]
        Q2_new[i] = Q2_half[i] - Δt / (2*Δx) * (p_right - p_left)
        Q3_new[i] = Q3_half[i] - Δt / (2*Δx) * (p_right * u_right - p_left * u_left) 
    end

    Q1_array .= Q1_new
    Q2_array .= Q2_new
    Q3_array .= Q3_new

    return Δt # for timekeeping
    
end
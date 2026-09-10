include("imp_data_ex.jl")
# ============================================================
# Estimate approximate Kc from the temporal slope of Y
# ============================================================
function estimate_Kc_from_slopes(K_vals, X, Y)

    M, N = size(Y)

    slopes = zeros(M)

    # Linear fit Y = a + b X for every K
    for i in 1:M
        x = vec(X)
        y = vec(Y[i, :])

        x̄ = mean(x)
        ȳ = mean(y)

        slopes[i] = sum((x .- x̄) .* (y .- ȳ)) /
                    sum((x .- x̄).^2)
    end

    # Look for a sign change in the slope
    crossings = Float64[]

    for i in 1:M-1
        s1 = slopes[i]
        s2 = slopes[i+1]

        if s1 * s2 <= 0 && s1 != s2
            # Linear interpolation of slope(K) = 0
            Kcross = K_vals[i] +
                     (K_vals[i+1] - K_vals[i]) *
                     (-s1) / (s2 - s1)

            push!(crossings, Kcross)
        end
    end

    if !isempty(crossings)
        # If several crossings exist, choose the one closest
        # to the K value having the smallest |slope|.
        Kref = K_vals[argmin(abs.(slopes))]
        Kc_est = crossings[argmin(abs.(crossings .- Kref))]
    else
        # No sign change: use the K with the smallest slope magnitude
        Kc_est = K_vals[argmin(abs.(slopes))]
    end

    return Kc_est, slopes
end


# ============================================================
# Generate synthetic curves around Kc
# ============================================================
function generate_sampling_K(Kc, ΔK, n_new, K_min, K_max)

    if n_new == 0
        return Float64[]
    end

    # Symmetric positions around Kc:
    #
    # n_new = 4, ΔK = 0.02
    # -> Kc - 0.04, Kc - 0.02, Kc + 0.02, Kc + 0.04
    #
    # n_new = 5
    # -> Kc - 0.04, Kc - 0.02, Kc,
    #    Kc + 0.02, Kc + 0.04
    offsets = (collect(1:n_new) .- (n_new + 1)/2) .* ΔK

    K_new = Kc .+ offsets

    # Do not generate points outside the existing experimental range
    K_new = K_new[(K_new .>= K_min) .&
                  (K_new .<= K_max)]

    return K_new
end


# ============================================================
# Main sampling function
# ============================================================
function finite_time_scaling_sampling(
    K_vals,
    t_vals,
    p2_mat,
    p2_err_mat;
    critic_estimate=true,
    Kc_input=nothing,
    ΔK=0.02,
    N_new=10,
    rng=Random.default_rng()
)

    # --------------------------------------------------------
    # Calculate X and Y from the experimental data
    # --------------------------------------------------------
    X, Y, Y_err = finite_time_scaling_data(
        t_vals,
        p2_mat,
        p2_err_mat
    )

    # --------------------------------------------------------
    # Estimate Kc or use supplied Kc
    # --------------------------------------------------------
    if critic_estimate

        Kc = estimate_Kc_from_slopes(
            K_vals,
            X,
            Y
        )

    else

        if Kc_input === nothing
            error("Kc_input must be provided when critic_estimate=false")
        end

        Kc = Kc_input
    end

    # --------------------------------------------------------
    # Generate new K values around Kc
    # --------------------------------------------------------
    K_new = generate_sampling_K(Kc, ΔK, N_new, minimum(K_vals), maximum(K_vals))

    # Keep only K values inside the experimental range
    valid = (
        (K_new .>= minimum(K_vals)) .&
        (K_new .<= maximum(K_vals))
    )

    K_new = K_new[valid]

    N_actual = length(K_new)

    if N_actual == 0
        error("No new K values fall inside the existing K range.")
    end

    # --------------------------------------------------------
    # Existing collapse
    #
    # Used ONLY to obtain the existing horizontal shifts,
    # which are needed to construct the synthetic master curve.
    # --------------------------------------------------------
    shifts, _ = finite_time_scaling(X, Y)

    # --------------------------------------------------------
    # Shift as a function of K
    # --------------------------------------------------------
    order = sortperm(K_vals)

    K_sorted = K_vals[order]
    shifts_sorted = shifts[order]

    function shift_from_K(K)

        K_clamped = clamp(
            K,
            minimum(K_sorted),
            maximum(K_sorted)
        )

        j = searchsortedlast(
            K_sorted,
            K_clamped
        )

        if j <= 1
            return shifts_sorted[1]

        elseif j >= length(K_sorted)
            return shifts_sorted[end]
        end

        K1 = K_sorted[j]
        K2 = K_sorted[j+1]

        a1 = shifts_sorted[j]
        a2 = shifts_sorted[j+1]

        α = (K_clamped - K1) / (K2 - K1)

        return (1 - α) * a1 + α * a2
    end

    # --------------------------------------------------------
    # Construct empirical master curve
    # --------------------------------------------------------
    X_shifted = X .+ shifts .* ones(1, length(t_vals))

    x_master = vec(X_shifted)
    y_master = vec(Y)
    err_master = vec(Y_err)

    order_master = sortperm(x_master)

    x_master = x_master[order_master]
    y_master = y_master[order_master]
    err_master = err_master[order_master]

    # --------------------------------------------------------
    # Master curve interpolation
    # --------------------------------------------------------
    function master_curve(x)

        x_clamped = clamp(
            x,
            minimum(x_master),
            maximum(x_master)
        )

        j = searchsortedlast(
            x_master,
            x_clamped
        )

        if j <= 1
            return y_master[1], err_master[1]

        elseif j >= length(x_master)
            return y_master[end], err_master[end]
        end

        x1 = x_master[j]
        x2 = x_master[j+1]

        y1 = y_master[j]
        y2 = y_master[j+1]

        e1 = err_master[j]
        e2 = err_master[j+1]

        α = (x_clamped - x1) / (x2 - x1)

        y = (1 - α) * y1 + α * y2
        e = (1 - α) * e1 + α * e2

        return y, e
    end

    # --------------------------------------------------------
    # Generate new synthetic data
    # --------------------------------------------------------
    p2_new = zeros(N_actual, length(t_vals))
    p2_err_new = zeros(N_actual, length(t_vals))

    for i in 1:N_actual

        K = K_new[i]

        # Shift corresponding to this K
        aK = shift_from_K(K)

        for j in eachindex(t_vals)

            # Scaling coordinate
            x = X[1, j] + aK

            # Synthetic master-curve value
            Y_mean, Y_error = master_curve(x)

            # Convert back from log(Λ)
            Λ_mean = exp(Y_mean)

            # Convert Λ back to <p²>
            p2_mean = Λ_mean * t_vals[j]^(2/3)

            # Typical relative experimental uncertainty
            relative_error = median(Y_err[:, j])

            σp2 = p2_mean * relative_error

            # Statistical realization
            p2_new[i, j] =
                max(
                    p2_mean + randn(rng) * σp2,
                    eps()
                )

            p2_err_new[i, j] = σp2
        end
    end

    # --------------------------------------------------------
    # Combine old and new data
    # --------------------------------------------------------

    K_total = vcat(K_vals, K_new)

    p2_total = vcat(
        p2_mat,
        p2_new
    )

    p2_err_total = vcat(
        p2_err_mat,
        p2_err_new
    )


    # --------------------------------------------------------
    # Sort everything according to increasing K
    # --------------------------------------------------------

    order = sortperm(K_total)

    K_total = K_total[order]

    p2_total = p2_total[order, :]

    p2_err_total = p2_err_total[order, :]


    # --------------------------------------------------------
    # Return complete dataset
    # --------------------------------------------------------

    return K_total, p2_total, p2_err_total
end

#Monte-Carlo sampling to estimate the quality of collapse: if repeated N_mc times, what is the mean and std of sX_rel?
function monte_carlo_sampling(
    K_vals,
    t_vals,
    p2_mat,
    p2_err_mat;
    critic_estimate=true,
    Kc_input=nothing,
    ΔK=0.02,
    N_new=10,
    N_mc=500,
    rng=MersenneTwister(1234)
)

    # Store sX_rel from each realization
    sX_rel_values = zeros(N_mc)

    # --------------------------------------------------------
    # Repeat the experiment N_mc times
    # --------------------------------------------------------
    for mc in 1:N_mc

        # Generate a new statistical realization
        K_new, p2_new, p2_err_new =
            finite_time_scaling_sampling(
                K_vals,
                t_vals,
                p2_mat,
                p2_err_mat;
                critic_estimate=critic_estimate,
                Kc_input=Kc_input,
                ΔK=ΔK,
                N_new=N_new,
                rng=rng
            )

        # ----------------------------------------------------
        # Add the new curves to the experimental data
        # ----------------------------------------------------
        K_total = vcat(
            K_vals,
            K_new
        )

        p2_total = vcat(
            p2_mat,
            p2_new
        )

        p2_err_total = vcat(
            p2_err_mat,
            p2_err_new
        )

        # ----------------------------------------------------
        # Your existing preprocessing
        # ----------------------------------------------------
        X_total, Y_total, Y_err_total =
            finite_time_scaling_data(
                t_vals,
                p2_total,
                p2_err_total
            )

        # ----------------------------------------------------
        # Your existing collapse
        # ----------------------------------------------------
        shifts, sX_rel =
            finite_time_scaling(
                X_total,
                Y_total
            )

        # Store collapse quality
        sX_rel_values[mc] = sX_rel
    end

    # --------------------------------------------------------
    # Monte-Carlo statistics
    # --------------------------------------------------------
    sX_rel_mean = mean(sX_rel_values)
    sX_rel_std = std(sX_rel_values)

    return (
        sX_rel_values = sX_rel_values,
        sX_rel_mean = sX_rel_mean,
        sX_rel_std = sX_rel_std
    )
end
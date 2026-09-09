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
    d=3,
    n_kicks_i=1,
    n_kicks_f=0,
    nbins=30,

    # Critical point
    critic_estimate=true,
    Kc_input=nothing,

    # Sampling
    ΔK=0.02,
    n_new=10,

    # Monte Carlo
    nmc=500,
    rng=MersenneTwister(1234))

    # --------------------------------------------------------
    # 1. Calculate X,Y,Yerr from the original data
    # --------------------------------------------------------
    X, Y, Yerr = finite_time_scaling_data(
        t_vals,
        p2_mat,
        p2_err_mat;
        d=d,
        n_kicks_i=n_kicks_i,
        n_kicks_f=n_kicks_f
    )

    # --------------------------------------------------------
    # 2. Estimate Kc, or use supplied value
    # --------------------------------------------------------
    if critic_estimate

        Kc_est, slopes = estimate_Kc_from_slopes(
            K_vals,
            X,
            Y
        )

    else

        if Kc_input === nothing
            error("Kc_input must be supplied when critic_estimate=false")
        end

        Kc_est = Kc_input

        # Still calculate slopes for diagnostics
        _, slopes = estimate_Kc_from_slopes(
            K_vals,
            X,
            Y
        )
    end

    # --------------------------------------------------------
    # 3. Generate the additional K values
    # --------------------------------------------------------
    K_new = generate_sampling_K(
        Kc_est,
        ΔK,
        n_new,
        minimum(K_vals),
        maximum(K_vals)
    )

    n_actual = length(K_new)

    if n_actual == 0
        error("No new K values were generated. Check Kc, ΔK and the K range.")
    end

    # --------------------------------------------------------
    # 4. Obtain the original collapse
    # --------------------------------------------------------
    shifts, sX_rel_original = finite_time_scaling(
        X,
        Y;
        nbins=nbins
    )

    # --------------------------------------------------------
    # 5. Build empirical master curve
    #
    #    Xshift = X + shift
    #    Y      = log(Λ)
    # --------------------------------------------------------
    Xshift = X .+ shifts

    x_master = vec(Xshift)
    y_master = vec(Y)
    err_master = vec(Yerr)

    # Sort according to shifted X
    order = sortperm(x_master)

    x_master = x_master[order]
    y_master = y_master[order]
    err_master = err_master[order]

    # --------------------------------------------------------
    # 6. Interpolation function for the empirical master curve
    # --------------------------------------------------------
    function master_curve(x)

        # Keep interpolation inside the experimentally sampled range
        x_clamped = clamp(
            x,
            minimum(x_master),
            maximum(x_master)
        )

        # Locate neighboring points
        j = searchsortedlast(x_master, x_clamped)

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

        y_interp = (1-α)*y1 + α*y2
        e_interp = (1-α)*e1 + α*e2

        return y_interp, e_interp
    end

    # --------------------------------------------------------
    # 7. Estimate shift for each new K
    #
    #    We interpolate shift(K) from the existing data.
    # --------------------------------------------------------
    shift_order = sortperm(K_vals)

    K_sorted = K_vals[shift_order]
    shift_sorted = shifts[shift_order]

    function shift_from_K(K)

        K_clamped = clamp(
            K,
            minimum(K_sorted),
            maximum(K_sorted)
        )

        j = searchsortedlast(K_sorted, K_clamped)

        if j <= 1
            return shift_sorted[1]
        elseif j >= length(K_sorted)
            return shift_sorted[end]
        end

        K1 = K_sorted[j]
        K2 = K_sorted[j+1]

        a1 = shift_sorted[j]
        a2 = shift_sorted[j+1]

        α = (K_clamped - K1) / (K2 - K1)

        return (1-α)*a1 + α*a2
    end

    # --------------------------------------------------------
    # 8. Generate ONE synthetic dataset
    # --------------------------------------------------------
    function generate_synthetic_dataset()

        p2_new = zeros(n_actual, length(t_vals))
        p2err_new = zeros(n_actual, length(t_vals))

        for i in 1:n_actual

            K = K_new[i]
            aK = shift_from_K(K)

            for j in eachindex(t_vals)

                # Original X coordinate:
                # X = -log(t)/d
                x = X[1,j] + aK

                # Empirical master curve
                y_mean, yerr = master_curve(x)

                # Λ = exp(Y)
                Λ_mean = exp(y_mean)

                # Convert back to <p²>
                p2_mean = Λ_mean * t_vals[j]^(2.0/d)

                # Relative error in p².
                #
                # We use a typical relative experimental error
                # from the existing data at this time.
                relerr_existing =
                    median(Yerr[:,j])

                σp2 = p2_mean * relerr_existing

                # Gaussian experimental noise
                p2_sample =
                    p2_mean + randn(rng) * σp2

                # Protect logarithm downstream
                p2_new[i,j] = max(p2_sample, eps())

                p2err_new[i,j] = σp2
            end
        end

        return p2_new, p2err_new #new generated datasets
    end

    # --------------------------------------------------------
    # 9. Monte-Carlo sampling of collapse quality
    # --------------------------------------------------------
    sX_rel_mc = zeros(nmc)

    for mc in 1:nmc

        p2_new, p2err_new =
            generate_synthetic_dataset()

        # Add synthetic curves to experimental data
        K_combined = vcat(K_vals, K_new)

        p2_combined = vcat(
            p2_mat,
            p2_new
        )

        p2err_combined = vcat(
            p2_err_mat,
            p2err_new
        )

        # Recalculate X,Y using exactly the same
        # finite-time-scaling preprocessing
        X_combined, Y_combined, _ =
            finite_time_scaling_data(
                t_vals,
                p2_combined,
                p2err_combined;
                d=d,
                n_kicks_i=n_kicks_i,
                n_kicks_f=n_kicks_f
            )

        # Perform the actual collapse
        _, sX_rel =
            finite_time_scaling(
                X_combined,
                Y_combined;
                nbins=nbins
            )
    
        sX_rel_mc[mc] = sX_rel
    end

    # --------------------------------------------------------
    # 10. Statistics of collapse quality
    # --------------------------------------------------------
    sX_rel_mean = mean(sX_rel_mc)
    sX_rel_std  = std(sX_rel_mc)

    return (
        Kc_est = Kc_est,
        K_new = K_new,

        original_sX_rel = sX_rel_original,

        sX_rel_mc = sX_rel_mc,
        sX_rel_mean = sX_rel_mean,
        sX_rel_std = sX_rel_std,

        slopes = slopes,
        shifts = shifts
    )
end
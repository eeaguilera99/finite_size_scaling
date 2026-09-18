include("fss_qkr3_timescaling.jl")
include("fss_qkr3_timescaling_analysis.jl")

using Interpolations
using Statistics
using Random


# ============================================================
# Estimate approximate Kc from temporal slope of Y
# ============================================================
function estimate_Kc_from_slopes(K_vals, X, Y)

    M, N = size(Y)

    slopes = zeros(M)

    # Linear fit Y = a + bX for every K curve
    for i in 1:M

        x = vec(X[1, :])
        y = vec(Y[i, :])

        x̄ = mean(x)
        ȳ = mean(y)

        slopes[i] =
            sum((x .- x̄) .* (y .- ȳ)) /
            sum((x .- x̄).^2)
    end


    # Look for sign changes in the slope
    crossings = Float64[]

    for i in 1:(M-1)

        s1 = slopes[i]
        s2 = slopes[i+1]

        if s1 * s2 <= 0 && s1 != s2

            Kcross =
                K_vals[i] +
                (K_vals[i+1] - K_vals[i]) *
                (-s1) / (s2 - s1)

            push!(crossings, Kcross)
        end
    end


    if !isempty(crossings)

        # Select the crossing closest to the K value
        # having the smallest absolute slope
        Kref = K_vals[argmin(abs.(slopes))]

        Kc_est =
            crossings[
                argmin(abs.(crossings .- Kref))
            ]

    else

        # No sign change:
        # use the curve with slope closest to zero
        Kc_est =
            K_vals[
                argmin(abs.(slopes))
            ]

    end

    return Kc_est, slopes
end

function estimate_delta_k(K_Vals; Kc_input=nothing)
    #estimate interval to sample around Kc
    if Kc_input === nothing
        Kc, _ = estimate_Kc_from_slopes(K_Vals, X, Y)
    else
        Kc = Kc_input
    end
    #closes value in K_vals to Kc
    Kc_index = argmin(abs.(K_Vals .- Kc))
    ΔK = K_Vals[Kc_index + 1] - K_Vals[Kc_index-1]

    return ΔK
end
# ============================================================
# Generate new K values around Kc
# ============================================================

function generate_sampling_K(Kc, ΔK, n_new, K_min, K_max)
    if n_new == 0
        return Float64[]
    end

    # Symmetric positions around Kc.
    #
    # For even n_new:
    #
    # n_new = 4
    # -> Kc-1.5ΔK, Kc-0.5ΔK,
    #    Kc+0.5ΔK, Kc+1.5ΔK
    #
    # For odd n_new the central point is exactly Kc.
    offsets =
        (collect(1:n_new) .- (n_new + 1)/2) .* ΔK

    K_new = Kc .+ offsets


    # Keep only values inside the experimental K range
    mask =
        (K_new .>= K_min) .&
        (K_new .<= K_max)

    K_new = K_new[mask]

    return K_new
end

# ============================================================
# Build one branch of the empirical master curve
#
# Localized and diffusive branches are treated independently.
# ============================================================
function build_master_curve_branch(
    X,
    Y,
    shifts,
    branch_mask;
    nbins=30,
    min_bin_points=2
)

    N = size(Y, 2)

    # Select only curves belonging to this branch
    Y_branch = Y[branch_mask, :]
    shifts_branch = shifts[branch_mask]

    if length(shifts_branch) < 2
        error(
            "Too few K curves on one side of Kc " *
            "to construct a master branch."
        )
    end


    # --------------------------------------------------------
    # Construct shifted X coordinates
    #
    # X is normally common to every K curve.
    # --------------------------------------------------------

    X_base = reshape(vec(X[1, :]), 1, N)

    X_branch =
        repeat(
            X_base,
            length(shifts_branch),
            1
        )

    X_shifted =
        X_branch .+
        reshape(shifts_branch, :, 1)


    x_raw = vec(X_shifted)
    y_raw = vec(Y_branch)


    # Remove invalid points
    valid =
        isfinite.(x_raw) .&
        isfinite.(y_raw)

    x_raw = x_raw[valid]
    y_raw = y_raw[valid]


    # --------------------------------------------------------
    # Divide shifted X into bins
    # --------------------------------------------------------

    x_min = minimum(x_raw)
    x_max = maximum(x_raw)

    bin_edges =
        range(
            x_min,
            x_max;
            length=nbins + 1
        )


    x_bin = Float64[]
    y_bin = Float64[]


    # --------------------------------------------------------
    # Average points inside every bin
    # --------------------------------------------------------

    for b in 1:nbins

        if b < nbins

            mask =
                (x_raw .>= bin_edges[b]) .&
                (x_raw .< bin_edges[b+1])

        else

            # Include upper boundary in final bin
            mask =
                (x_raw .>= bin_edges[b]) .&
                (x_raw .<= bin_edges[b+1])

        end


        if count(mask) >= min_bin_points

            push!(
                x_bin,
                mean(x_raw[mask])
            )

            push!(
                y_bin,
                mean(y_raw[mask])
            )
        end
    end


    # --------------------------------------------------------
    # Sort bins according to increasing X
    # --------------------------------------------------------

    order = sortperm(x_bin)

    x_bin = x_bin[order]
    y_bin = y_bin[order]


    if length(x_bin) < 4
        error(
            "Too few populated bins to construct " *
            "this master-curve branch."
        )
    end


    if any(diff(x_bin) .<= 0)
        error(
            "Binned X coordinates must be strictly increasing."
        )
    end


    # --------------------------------------------------------
    # Shape-preserving interpolation of this branch
    # --------------------------------------------------------

    master_itp = extrapolate(
        interpolate(
            x_bin,
            y_bin,
            SteffenMonotonicInterpolation()
        ),
        Interpolations.Line()
    )


    return master_itp, x_bin, y_bin
end



# ============================================================
# Build interpolation of 1/xi(K) for ONE side of Kc
#
# q(K) = 1/xi(K) = 1/exp(shift)
#
# This follows your current implementation.
# ============================================================
function build_inverse_xi_interpolation(
    K_branch,
    shifts_branch
)

    if length(K_branch) < 2
        error(
            "At least two K values are required " *
            "to interpolate 1/xi on each side of Kc."
        )
    end


    order = sortperm(K_branch)

    K_sorted =
        K_branch[order]

    shifts_sorted =
        shifts_branch[order]


    # q = 1/xi = exp(-shift)
    q_sorted =
        1.0 ./ exp.(shifts_sorted)


    # K must be strictly increasing
    if any(diff(K_sorted) .<= 0)
        error(
            "K values used for interpolation " *
            "must be strictly increasing."
        )
    end


    q_itp = extrapolate(
        interpolate(
            K_sorted,
            q_sorted,
            SteffenMonotonicInterpolation()
        ),
        Interpolations.Line()
    )

    return q_itp
end



# ============================================================
# Main synthetic-data sampling function
# ============================================================
function finite_time_scaling_sampling(
    K_vals,
    t_vals,
    p2_mat,
    p2_err_mat;
    critic_estimate=false,
    Kc_input=nothing,
    ΔK=nothing,
    N_new=10,
    master_nbins=30,
    rng=Random.default_rng()
)

    # --------------------------------------------------------
    # Original scaling data
    # --------------------------------------------------------

    X, Y, Y_err =
        finite_time_scaling_data(
            t_vals,
            p2_mat,
            p2_err_mat
        )


    if ΔK === nothing
        ΔK = estimate_delta_k(K_vals; Kc_input=Kc_input)
    end
    # --------------------------------------------------------
    # Estimate or specify Kc
    # --------------------------------------------------------

    if critic_estimate

        Kc, _ =
            estimate_Kc_from_slopes(
                K_vals,
                X,
                Y
            )

    else

        if Kc_input === nothing
            error(
                "Kc_input must be supplied when " *
                "critic_estimate=false."
            )
        end

        Kc = Kc_input
    end


    # --------------------------------------------------------
    # Generate new K values
    # --------------------------------------------------------

    K_new =
        generate_sampling_K(
            Kc,
            ΔK,
            N_new,
            minimum(K_vals),
            maximum(K_vals)
        )


    # Do not generate a curve exactly at Kc.
    #
    # The current two-branch model does not define
    # a separate critical master curve.
    tol_Kc =
        100 * eps(Float64) *
        max(1.0, abs(Kc))

    K_new =
        K_new[
            abs.(K_new .- Kc) .> tol_Kc
        ]


    N_actual = length(K_new)


    if N_actual == 0
        error(
            "No valid new K values were generated. " *
            "Use an even N_new or modify ΔK."
        )
    end


    # --------------------------------------------------------
    # Collapse ORIGINAL dataset
    #
    # Used only to obtain the original shifts.
    # --------------------------------------------------------

    shifts, _ =
        finite_time_scaling(
            X,
            Y,
            nbins=30
        )


    # --------------------------------------------------------
    # Split original curves into two physical branches
    # --------------------------------------------------------

    loc_mask =
        K_vals .<= Kc

    diff_mask =
        K_vals .>= Kc


    if count(loc_mask) < 2
        error(
            "Too few original curves below Kc."
        )
    end

    if count(diff_mask) < 2
        error(
            "Too few original curves above Kc."
        )
    end


    # ========================================================
    # Interpolation of 1/xi(K) separately on the two sides
    # ========================================================

    q_loc_itp =
        build_inverse_xi_interpolation(
            K_vals[loc_mask],
            shifts[loc_mask]
        )


    q_diff_itp =
        build_inverse_xi_interpolation(
            K_vals[diff_mask],
            shifts[diff_mask]
        )
    #=    
    plt1 = plot(title=latexstring("ξ(K) interpolation, \$d=$(D)\$"), xlabel=L"κ", ylabel=L"ξ")
    plot!(plt1, K_vals[loc_mask], log.(1 ./ q_loc_itp.(K_vals[loc_mask])), label="Original data (localized)", color=:blue)
    plot!(plt1, K_vals[diff_mask], log.(1 ./ q_diff_itp.(K_vals[diff_mask])), label="Original data (diffusive)", color=:green)
    display(plt1)
    =#

    # --------------------------------------------------------
    # Convert interpolated q(K)=1/xi back into shift:
    #
    # shift = ln(xi) = ln(1/q)
    # --------------------------------------------------------

    function shift_from_K(K)

        if K < Kc

            qK = q_loc_itp(K)

        elseif K > Kc

            qK = q_diff_itp(K)

        else

            error(
                "Shift interpolation exactly at Kc " *
                "is not defined."
            )
        end


        # 1/xi must remain positive.
        #
        # Line extrapolation close to Kc can in principle
        # cross zero. Do not silently accept an unphysical value.
        if qK <= 0

            error(
                "Interpolation of 1/xi produced a non-positive " *
                "value at K = $K. " *
                "Reduce the extrapolation distance or inspect " *
                "the 1/xi(K) interpolation."
            )
        end


        return log(1.0 / qK)
    end



    # ========================================================
    # Construct TWO empirical master curves
    # ========================================================

    master_loc,
    x_loc_bin,
    y_loc_bin =
        build_master_curve_branch(
            X,
            Y,
            shifts,
            loc_mask;
            nbins=master_nbins
        )


    master_diff,
    x_diff_bin,
    y_diff_bin =
        build_master_curve_branch(
            X,
            Y,
            shifts,
            diff_mask;
            nbins=master_nbins
        )



    # ========================================================
    # Generate synthetic curves
    # ========================================================

    p2_new =
        zeros(
            N_actual,
            length(t_vals)
        )

    p2_err_new =
        zeros(
            N_actual,
            length(t_vals)
        )


    for i in 1:N_actual

        K = K_new[i]


        # ----------------------------------------------------
        # Select branch and corresponding shift interpolation
        # ----------------------------------------------------

        if K < Kc

            aK =
                shift_from_K(K)

            master_itp =
                master_loc

        elseif K > Kc

            aK =
                shift_from_K(K)

            master_itp =
                master_diff

        else

            error(
                "Synthetic data exactly at Kc are not defined " *
                "by the present two-branch construction."
            )
        end


        # ----------------------------------------------------
        # Generate all time points for this K
        # ----------------------------------------------------

        for j in eachindex(t_vals)

            # Original unshifted X
            x0 = X[1, j]

            # Shifted scaling coordinate
            x =
                x0 + aK


            # Correct branch of empirical master curve
            Y_mean =
                master_itp(x)


            # Y = ln Λ
            Λ_mean =
                exp(Y_mean)


            # Λ = <p²>/t^(2/3)
            p2_mean =
                Λ_mean *
                t_vals[j]^(2/3)


            # ------------------------------------------------
            # Experimental uncertainty
            #
            # Since Y_err ≈ σ_p2 / p2,
            # use the median experimental relative error
            # at this time.
            # ------------------------------------------------

            relative_error =
                median(
                    Y_err[:, j]
                )


            σp2 =
                p2_mean *
                relative_error


            # ------------------------------------------------
            # Synthetic experimental measurement
            # ------------------------------------------------

            sampled_value =
                p2_mean +
                randn(rng) * σp2


            p2_new[i, j] =
                max(
                    sampled_value,
                    eps(Float64)
                )


            p2_err_new[i, j] =
                σp2
        end
    end



    # ========================================================
    # Combine original and synthetic data
    # ========================================================

    K_total =
        vcat(
            K_vals,
            K_new
        )


    p2_total =
        vcat(
            p2_mat,
            p2_new
        )


    p2_err_total =
        vcat(
            p2_err_mat,
            p2_err_new
        )



    # --------------------------------------------------------
    # Sort all rows according to increasing K
    # --------------------------------------------------------

    order_total =
        sortperm(K_total)


    K_total =
        K_total[order_total]


    p2_total =
        p2_total[
            order_total,
            :
        ]


    p2_err_total =
        p2_err_total[
            order_total,
            :
        ]


    println(
        "Data sampling generated around Kc = $Kc " *
        "with ΔK = $ΔK and N_new = $N_actual."
    )


    # --------------------------------------------------------
    # FINAL OUTPUT
    # --------------------------------------------------------

    return (
        K_total,
        p2_total,
        p2_err_total
    )
end



# ============================================================
# Monte-Carlo estimate of collapse-quality statistics
#
# IMPORTANT:
# finite_time_scaling_sampling now ALREADY returns the
# complete old + new dataset.
# ============================================================
function monte_carlo_sampling(
    K_vals,
    t_vals,
    p2_mat,
    p2_err_mat;
    critic_estimate=true,
    Kc_input=nothing,
    ΔK=0.02,
    N_new=10,
    master_nbins=30,
    N_mc=500,
    rng=MersenneTwister(1234)
)

    sX_rel_values =
        zeros(N_mc)


    for mc in 1:N_mc

        # ----------------------------------------------------
        # This already returns ORIGINAL + SYNTHETIC data
        # ----------------------------------------------------

        K_total,
        p2_total,
        p2_err_total =
            finite_time_scaling_sampling(
                K_vals,
                t_vals,
                p2_mat,
                p2_err_mat;
                critic_estimate=critic_estimate,
                Kc_input=Kc_input,
                ΔK=ΔK,
                N_new=N_new,
                master_nbins=master_nbins,
                rng=rng
            )


        # ----------------------------------------------------
        # Standard FSS preprocessing
        # ----------------------------------------------------

        X_total,
        Y_total,
        Y_err_total =
            finite_time_scaling_data(
                t_vals,
                p2_total,
                p2_err_total
            )


        # ----------------------------------------------------
        # Standard collapse
        # ----------------------------------------------------

        shifts_total,
        sX_rel =
            finite_time_scaling(
                X_total,
                Y_total
            )


        sX_rel_values[mc] =
            sX_rel
    end


    sX_rel_mean =
        mean(sX_rel_values)

    sX_rel_std =
        std(sX_rel_values)


    return (
        sX_rel_values=sX_rel_values,
        sX_rel_mean=sX_rel_mean,
        sX_rel_std=sX_rel_std
    )
end
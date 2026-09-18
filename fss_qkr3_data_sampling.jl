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
    x = vec(X[1, :])
    x̄ = mean(x)
    denominator = sum((x .- x̄).^2)
    denominator > 0 || error("At least two distinct times are required.")

    for i in 1:M
        y = vec(Y[i, :])
        ȳ = mean(y)
        slopes[i] = sum((x .- x̄) .* (y .- ȳ)) / denominator
    end

    # Search neighboring K values even if the input is unsorted.
    order = sortperm(K_vals)
    crossings = Float64[]
    for j in 1:(M - 1)
        i1, i2 = order[j], order[j + 1]
        s1, s2 = slopes[i1], slopes[i2]
        if s1 * s2 <= 0 && s1 != s2
            Kcross = K_vals[i1] +
                     (K_vals[i2] - K_vals[i1]) * (-s1) / (s2 - s1)
            push!(crossings, Kcross)
        end
    end

    Kref = K_vals[argmin(abs.(slopes))]
    Kc_est = isempty(crossings) ? Kref :
             crossings[argmin(abs.(crossings .- Kref))]
    return Kc_est, slopes
end

function estimate_delta_k(K_vals; Kc_input=nothing, X=nothing, Y=nothing)
    if Kc_input === nothing
        (X === nothing || Y === nothing) &&
            error("Supply Kc_input, or supply both X and Y.")
        Kc, _ = estimate_Kc_from_slopes(K_vals, X, Y)
    else
        Kc = Kc_input
    end

    K_sorted = sort(collect(K_vals))
    length(K_sorted) >= 2 || error("At least two K values are required.")
    any(diff(K_sorted) .<= 0) && error("K values must be distinct.")
    i = argmin(abs.(K_sorted .- Kc))

    # Preserve the original two-neighbor interval in the interior.
    # At an endpoint, use the available neighboring interval.
    return K_sorted[min(i + 1, end)] - K_sorted[max(i - 1, 1)]
end

# ============================================================
# Generate new K values around Kc
# ============================================================
function generate_sampling_K(Kc, ΔK, n_new, K_min, K_max)
    (n_new isa Integer && n_new >= 0) ||
        error("n_new must be a non-negative integer.")
    n_new == 0 && return Float64[]
    (ΔK isa Real && isfinite(ΔK) && ΔK > 0) ||
        error("ΔK must be finite and positive.")

    offsets = (collect(1:n_new) .- (n_new + 1)/2) .* ΔK
    K_new = Kc .+ offsets
    mask = (K_new .>= K_min) .& (K_new .<= K_max)
    return K_new[mask]
end

# ============================================================
# Generate final time grid: original times + earlier/later times
#
# t_step controls ONLY the new points; measured times are kept.
# Earlier points start at t_early (at t_step if t_early == 0).
# Later points start at maximum(t_vals) + t_step.
# A bound is included only if it lies on the generated grid.
# t = 0 is excluded because the scaling coordinate contains log(t).
# ============================================================
function generate_sampling_t(
    t_vals;
    t_early=nothing,
    t_late=nothing,
    Δt=50
)
    isempty(t_vals) && error("t_vals must not be empty.")
    all(t -> isfinite(t) && t > 0, t_vals) ||
        error("Original times must be finite and strictly positive.")
    (isfinite(Δt) && Δt > 0) ||
        error("Δt must be finite and positive.")

    t_min, t_max = extrema(t_vals)
    new_early = Float64[]
    new_late = Float64[]

    if t_early !== nothing
        (isfinite(t_early) && t_early >= 0) ||
            error("t_early must be finite and non-negative.")

        if t_early < t_min
            t_start = iszero(t_early) ? Δt : t_early
            new_early = [
                t for t in t_start:Δt:t_min
                if 0 < t < t_min
            ]
        end
    end

    if t_late !== nothing
        (isfinite(t_late) && t_late > 0) ||
            error("t_late must be finite and positive.")

        if t_late > t_max
            new_late = collect((t_max + Δt):Δt:t_late)
        end
    end

    return sort!(unique(vcat(
        collect(t_vals),
        new_early,
        new_late
    )))
end

# ============================================================
# Build one branch of the empirical master curve
# Localized and diffusive branches are treated independently.
# ============================================================
function build_master_curve_branch(
    X, Y, shifts, branch_mask;
    nbins=30,
    min_bin_points=2
)
    N = size(Y, 2)
    Y_branch = Y[branch_mask, :]
    shifts_branch = shifts[branch_mask]

    length(shifts_branch) >= 2 ||
        error(
            "Too few K curves on one side of Kc " *
            "to construct a master branch."
        )

    X_base = reshape(vec(X[1, :]), 1, N)
    X_branch = repeat(X_base, length(shifts_branch), 1)
    X_shifted = X_branch .+ reshape(shifts_branch, :, 1)

    x_raw, y_raw = vec(X_shifted), vec(Y_branch)
    valid = isfinite.(x_raw) .& isfinite.(y_raw)
    x_raw, y_raw = x_raw[valid], y_raw[valid]

    isempty(x_raw) &&
        error("No valid points in this master-curve branch.")

    x_min, x_max = extrema(x_raw)
    bin_edges = range(x_min, x_max; length=nbins + 1)
    x_bin, y_bin = Float64[], Float64[]

    for b in 1:nbins
        if b < nbins
            mask =
                (x_raw .>= bin_edges[b]) .&
                (x_raw .< bin_edges[b + 1])
        else
            mask =
                (x_raw .>= bin_edges[b]) .&
                (x_raw .<= bin_edges[b + 1])
        end

        if count(mask) >= min_bin_points
            push!(x_bin, mean(x_raw[mask]))
            push!(y_bin, mean(y_raw[mask]))
        end
    end

    order = sortperm(x_bin)
    x_bin, y_bin = x_bin[order], y_bin[order]

    length(x_bin) >= 4 ||
        error("Too few populated bins to construct this master-curve branch.")

    any(diff(x_bin) .<= 0) &&
        error("Binned X coordinates must be strictly increasing.")

    # Preserve the original Steffen interpolation + Line extrapolation.
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
# q(K) = 1/xi(K) = 1/exp(shift)
# ============================================================
function build_inverse_xi_interpolation(K_branch, shifts_branch)
    length(K_branch) >= 2 ||
        error(
            "At least two K values are required " *
            "to interpolate 1/xi on each side of Kc."
        )

    order = sortperm(K_branch)
    K_sorted = K_branch[order]
    shifts_sorted = shifts_branch[order]

    q_sorted = 1.0 ./ exp.(shifts_sorted)

    any(diff(K_sorted) .<= 0) &&
        error("K values used for interpolation must be strictly increasing.")

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
#
# N_new = 0: no additional K curves.
# t_new = nothing: return the original t_vals, in original order.
# Otherwise: sorted union of t_vals and t_new.
# Original measurements and their errors are copied unchanged.
# Always returns K_total, t_total, p2_total, p2_err_total.
# ============================================================
function finite_time_scaling_sampling(
    K_vals, t_vals, p2_mat, p2_err_mat;
    critic_estimate=false,
    Kc_input=nothing,
    ΔK=nothing,
    N_new=10,
    t_new=nothing,
    master_nbins=30,
    rng=Random.default_rng()
)
    M, N = length(K_vals), length(t_vals)

    M > 0 && N > 0 ||
        error("K_vals and t_vals must not be empty.")

    size(p2_mat) == (M, N) ||
        error("p2_mat must have size (length(K_vals), length(t_vals)).")

    size(p2_err_mat) == (M, N) ||
        error("p2_err_mat must have the same size as p2_mat.")

    all(isfinite, K_vals) ||
        error("K values must be finite.")

    length(unique(K_vals)) == M ||
        error("K values must be distinct.")

    all(t -> isfinite(t) && t > 0, t_vals) ||
        error("Times must be finite and positive.")

    length(unique(t_vals)) == N ||
        error("Original times must be distinct.")

    all(p -> isfinite(p) && p > 0, p2_mat) ||
        error("p2 values must be finite and positive.")

    all(e -> isfinite(e) && e >= 0, p2_err_mat) ||
        error("Errors must be finite and non-negative.")

    (N_new isa Integer && N_new >= 0) ||
        error("N_new must be a non-negative integer.")

    # --------------------------------------------------------
    # Final time grid and mapping back to measured columns
    # --------------------------------------------------------
    if t_new === nothing
        t_total = t_vals
    else
        t_new isa AbstractVector ||
            error("t_new must be a vector or nothing.")

        all(t -> isfinite(t) && t > 0, t_new) ||
            error(
                "New times must be finite and positive; " *
                "t = 0 is not allowed."
            )

        t_total = sort!(unique(vcat(
            collect(t_vals),
            collect(t_new)
        )))
    end

    original_column = Dict(
        t => j for (j, t) in enumerate(t_vals)
    )

    old_columns = [
        get(original_column, t, 0) for t in t_total
    ]

    added_time_indices = findall(iszero, old_columns)

    if N_new == 0 && isempty(added_time_indices)
        order = sortperm(K_vals)
        return (
            K_vals[order],
            t_total,
            p2_mat[order, old_columns],
            p2_err_mat[order, old_columns]
        )
    end

    N >= 2 ||
        error("At least two original times are required for sampling.")

    # --------------------------------------------------------
    # Original scaling data; estimate or specify Kc
    # --------------------------------------------------------
    X, Y, Y_err = finite_time_scaling_data(
        t_vals,
        p2_mat,
        p2_err_mat
    )

    if critic_estimate
        Kc, _ = estimate_Kc_from_slopes(K_vals, X, Y)
    else
        Kc_input === nothing &&
            error("Kc_input must be supplied when critic_estimate=false.")
        Kc = Kc_input
    end

    isfinite(Kc) || error("Kc must be finite.")

    # --------------------------------------------------------
    # Generate new K values; omit Kc and existing curves
    # --------------------------------------------------------
    if N_new > 0
        if ΔK === nothing
            ΔK = estimate_delta_k(K_vals; Kc_input=Kc)
        end

        K_new = generate_sampling_K(
            Kc,
            ΔK,
            N_new,
            minimum(K_vals),
            maximum(K_vals)
        )
    else
        K_new = Float64[]
    end

    tol_Kc = 100 * eps(Float64) * max(1.0, abs(Kc))

    K_new = unique(filter(
        K ->
            abs(K - Kc) > tol_Kc &&
            !any(Kold -> abs(K - Kold) <= tol_Kc, K_vals),
        K_new
    ))

    N_actual = length(K_new)

    if N_new > 0 && N_actual == 0
        @warn "No distinct noncritical K values were generated; continuing with the requested time grid."
    end

    if N_actual == 0 && isempty(added_time_indices)
        order = sortperm(K_vals)
        return (
            K_vals[order],
            t_total,
            p2_mat[order, old_columns],
            p2_err_mat[order, old_columns]
        )
    end

    # --------------------------------------------------------
    # Collapse ONLY the original dataset to obtain its shifts
    # --------------------------------------------------------
    shifts, _ = finite_time_scaling(X, Y; nbins=30)
    shifts = vec(shifts)

    all(isfinite, shifts) ||
        error("The original collapse produced invalid shifts.")

    loc_mask = K_vals .<= Kc
    diff_mask = K_vals .>= Kc

    count(loc_mask) >= 2 ||
        error("Too few original curves below Kc.")

    count(diff_mask) >= 2 ||
        error("Too few original curves above Kc.")

    # Interpolate inverse xi only for newly generated K curves.
    q_loc_itp = q_diff_itp = nothing

    if N_actual > 0
        q_loc_itp = build_inverse_xi_interpolation(
            K_vals[loc_mask],
            shifts[loc_mask]
        )

        q_diff_itp = build_inverse_xi_interpolation(
            K_vals[diff_mask],
            shifts[diff_mask]
        )
    end

    function shift_from_K(K)
        if K < Kc
            qK = q_loc_itp(K)
        elseif K > Kc
            qK = q_diff_itp(K)
        else
            error("Shift interpolation exactly at Kc is not defined.")
        end

        (isfinite(qK) && qK > 0) || error(
            "Interpolation of 1/xi produced an invalid value at K = $K. " *
            "Reduce the extrapolation distance or inspect the interpolation."
        )

        return log(1.0 / qK)
    end

    # --------------------------------------------------------
    # Construct the same TWO empirical master curves
    # --------------------------------------------------------
    master_loc, x_loc_bin, y_loc_bin = build_master_curve_branch(
        X,
        Y,
        shifts,
        loc_mask;
        nbins=master_nbins
    )

    master_diff, x_diff_bin, y_diff_bin = build_master_curve_branch(
        X,
        Y,
        shifts,
        diff_mask;
        nbins=master_nbins
    )

    # --------------------------------------------------------
    # Experimental relative errors on the final time grid
    # Preserve the original median across K at measured times.
    # For new times use Steffen + Line, without Flat or clipping.
    # --------------------------------------------------------
    rel_err_original = [
        median(Y_err[:, j]) for j in 1:N
    ]

    all(e -> isfinite(e) && e >= 0, rel_err_original) ||
        error("The original relative errors must be finite and non-negative.")

    rel_err_total = zeros(length(t_total))

    for j in eachindex(t_total)
        j_old = old_columns[j]
        if j_old != 0
            rel_err_total[j] = rel_err_original[j_old]
        end
    end

    if !isempty(added_time_indices)
        order_t = sortperm(t_vals)

        err_itp = extrapolate(
            interpolate(
                collect(t_vals[order_t]),
                rel_err_original[order_t],
                SteffenMonotonicInterpolation()
            ),
            Interpolations.Line()
        )

        for j in added_time_indices
            rel_err_total[j] = err_itp(t_total[j])
        end
    end

    all(e -> isfinite(e) && e >= 0, rel_err_total) || error(
        "Relative-error extrapolation produced a negative or non-finite value. " *
        "Reduce the requested time range or revise the uncertainty model."
    )

    # --------------------------------------------------------
    # Combine original and synthetic rows on the final grid
    # --------------------------------------------------------
    K_total = vcat(K_vals, K_new)

    p2_total = zeros(
        M + N_actual,
        length(t_total)
    )

    p2_err_total = zeros(
        M + N_actual,
        length(t_total)
    )

    for j in eachindex(t_total)
        j_old = old_columns[j]

        if j_old != 0
            p2_total[1:M, j] .= p2_mat[:, j_old]
            p2_err_total[1:M, j] .= p2_err_mat[:, j_old]
        end
    end

    for i in eachindex(K_total)
        is_original = i <= M

        time_indices =
            is_original ? added_time_indices : eachindex(t_total)

        isempty(time_indices) && continue

        K = K_total[i]

        # Exact fitted shift for measured K; interpolated shift for new K.
        aK = is_original ? shifts[i] : shift_from_K(K)

        if K < Kc
            master_itp = master_loc

        elseif K > Kc
            master_itp = master_diff

        else
            # An existing measured curve can lie exactly at Kc.
            # Extend its own Y(X+aK) with the same Steffen + Line choice,
            # rather than arbitrarily assigning it to either master branch.
            # No additional K curve at Kc is generated.
            x_critical = vec(X[1, :]) .+ aK
            order_x = sortperm(x_critical)

            master_itp = extrapolate(
                interpolate(
                    x_critical[order_x],
                    vec(Y[i, :])[order_x],
                    SteffenMonotonicInterpolation()
                ),
                Interpolations.Line()
            )
        end

        for j in time_indices
            t = t_total[j]
            j_old = old_columns[j]

            x0 = j_old == 0 ? -log(t)/3 : X[1, j_old]
            Y_mean = master_itp(x0 + aK)

            p2_mean = exp(Y_mean) * t^(2/3)

            (isfinite(p2_mean) && p2_mean > 0) ||
                error("Invalid synthetic mean at K = $K, t = $t.")

            σp2 = p2_mean * rel_err_total[j]

            isfinite(σp2) ||
                error("Invalid uncertainty at K = $K, t = $t.")

            sampled_value = p2_mean + randn(rng) * σp2

            isfinite(sampled_value) ||
                error("Invalid sampled value at K = $K, t = $t.")

            p2_total[i, j] = max(
                sampled_value,
                eps(Float64)
            )

            p2_err_total[i, j] = σp2
        end
    end

    # --------------------------------------------------------
    # Sort rows by K; time columns already match t_total
    # --------------------------------------------------------
    order_total = sortperm(K_total)

    K_total = K_total[order_total]
    p2_total = p2_total[order_total, :]
    p2_err_total = p2_err_total[order_total, :]

    println(
        "Data sampling generated around Kc = $Kc " *
        "with ΔK = $ΔK, N_new = $N_actual, " *
        "and $(length(added_time_indices)) new time points."
    )

    return K_total, t_total, p2_total, p2_err_total
end

# ============================================================
# Monte-Carlo estimate of collapse-quality statistics
# The sampling function returns the COMPLETE old + new dataset.
# ============================================================
function monte_carlo_sampling(
    K_vals, t_vals, p2_mat, p2_err_mat;
    critic_estimate=true,
    Kc_input=nothing,
    ΔK=0.02,
    N_new=10,
    t_new=nothing,
    master_nbins=30,
    N_mc=500,
    rng=MersenneTwister(1234)
)
    (N_mc isa Integer && N_mc >= 2) ||
        error("N_mc must be an integer of at least 2.")

    sX_rel_values = zeros(N_mc)

    for mc in 1:N_mc
        K_total, t_total, p2_total, p2_err_total =
            finite_time_scaling_sampling(
                K_vals,
                t_vals,
                p2_mat,
                p2_err_mat;
                critic_estimate=critic_estimate,
                Kc_input=Kc_input,
                ΔK=ΔK,
                N_new=N_new,
                t_new=t_new,
                master_nbins=master_nbins,
                rng=rng
            )

        X_total, Y_total, Y_err_total =
            finite_time_scaling_data(
                t_total,
                p2_total,
                p2_err_total
            )

        shifts_total, sX_rel =
            finite_time_scaling(
                X_total,
                Y_total
            )

        sX_rel_values[mc] = sX_rel
    end

    return (
        sX_rel_values=sX_rel_values,
        sX_rel_mean=mean(sX_rel_values),
        sX_rel_std=std(sX_rel_values)
    )
end
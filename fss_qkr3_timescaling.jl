using LinearAlgebra
using Statistics
using Optim
using Plots
using LsqFit
using CSV, DataFrames
using LaTeXStrings

K_vals = vec(Matrix(CSV.read("data2/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("data2/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("data2/nc_matrix.csv", DataFrame; header=false))             # ⟨p²⟩ values
p2_err_mat = Matrix(CSV.read("data2/nc_err_matrix.csv", DataFrame; header=false))     # Errors



function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=100, d, n_kicks_i=1, n_kicks_f=0)
    
    #filter Nkicks range
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]

    M, N = size(p2_mat)

    # Observable: Λ = <p^2>/t^(2/3)
    Λ = p2_mat ./ (t_vals' .^ (2/d))
    Λ_err = p2_err_mat ./ (t_vals' .^ (2/d))

    # Log variables
    X = -log.(t_vals' .^ (1/d))     # 1×N
    Y = log.(Λ)                # M×N
    # Propagate errors: Δ(ln Λ) ≈ ΔΛ / Λ
    Yerr = Λ_err ./ Λ

    # Flatten for binning
    allY = vec(Y)
    y_min, y_max = minimum(allY), maximum(allY) #range of data Yaxis
    bins = range(y_min, y_max; length=nbins+1)
    bin_ids = [searchsortedlast(bins, y) for y in allY] #assigns each y-value to a bin number from 1 to nbins

    # Cost function: variance of shifted X within Y-bins
    function cost(a_full::Vector)
        # Shift X values by the amount specified in a_full
        shiftedX = vec(X .+ a_full .* ones(1,N))
        
        total_var = 0.0
        for b in 1:nbins
            # Find all points that fall in bin b
            mask = (bin_ids .== b)
            
            # Only consider bins with more than 1 point
            if count(mask) > 1
                # Get X values for points in this bin
                xb = shiftedX[mask]
                # Add weighted variance of X values in this bin
                total_var += var(xb) * count(mask)
            end
        end
        return total_var
    end

    # Fix gauge: a₁ = 0
    function constrained_cost(a_free::Vector)
        a_full = vcat(0.0, a_free)
        return cost(a_full)
    end

    # Initial guess
    a0 = zeros(M-1)

    # Minimize
    res = optimize(constrained_cost, a0, NelderMead())
    shifts = vcat(0.0, Optim.minimizer(res))

    # === Compute normalized scatter directly from res.minimum ===
    total_points = M * N
    sX = sqrt(res.minimum / total_points)
    Xp = X .+ shifts                     # shifted X matrix
    sX_rel = sX / (maximum(Xp) - minimum(Xp) + eps())

    # === Return everything
    return res, shifts, X, Y, Yerr, sX_rel
end

function tot_variance(a_full::Vector, X::Matrix, Y::Matrix; nbins=100)# calculates rel var for arbitrary shifts
    Xp = vec(X .+ a_full .* ones(1,size(X,2)))
    
    # Flatten for binning
    allY = vec(Y)
    y_min, y_max = minimum(allY), maximum(allY) #range of data Yaxis
    bins = range(y_min, y_max; length=nbins+1)
    bin_ids = [searchsortedlast(bins, y) for y in allY]

    totw, totvar = 0.0, 0.0
    for b in 1:nbins
        mask = (bin_ids .== b)
        nb = count(mask)
        if nb > 1
            xb = vec(Xp)[mask]
            totvar += var(xb) * nb
            totw   += nb
        end
    end
    sX = sqrt(totvar / max(totw, 1.0))
    sX_rel = sX / (maximum(vec(Xp)) - minimum(vec(Xp)) + eps())

    return totvar, sX_rel
end

# === corrected collapse with irrelevant scaling variable =====================

# Lightweight weighted linear regression: y = a + b*x
# Returns (a, b, sse, ok)
function _wlinfit(x::AbstractVector, y::AbstractVector; w::AbstractVector=ones(length(x)))
    n = length(x)
    if n < 2 || any(!isfinite, x) || any(!isfinite, y) || any(!isfinite, w)
        return (NaN, NaN, Inf, false)
    end
    W   = sum(w)
    Wx  = sum(w .* x)
    Wy  = sum(w .* y)
    Wxx = sum(w .* x .* x)
    Wxy = sum(w .* x .* y)
    denom = W*Wxx - Wx*Wx
    if abs(denom) ≤ 1e-14
        return (NaN, NaN, Inf, false)
    end
    b = (W*Wxy - Wx*Wy) / denom
    a = (Wy - b*Wx) / W
    yhat = a .+ b .* x
    sse  = sum(w .* (y .- yhat).^2)
    return (a, b, sse, true)
end

# Minimal linear interpolator for (u, value); extrapolates linearly at the ends.
struct _LinInterp
    u::Vector{Float64}
    v::Vector{Float64}
end
function _LinInterp(u::AbstractVector, v::AbstractVector)
    us, vs = collect(u), collect(v)
    p = sortperm(us)
    _LinInterp(us[p], vs[p])
end
function (li::_LinInterp)(x::AbstractVector)
    u, v = li.u, li.v
    n = length(u)
    out = similar(x, Float64)
    for (k,xx) in enumerate(x)
        if xx ≤ u[1]
            out[k] = n==1 ? v[1] : (v[1] + (xx-u[1])*(v[2]-v[1])/(u[2]-u[1]))
        elseif xx ≥ u[end]
            out[k] = n==1 ? v[end] : (v[end] + (xx-u[end])*(v[end]-v[end-1])/(u[end]-u[end-1]))
        else
            j = searchsortedfirst(u, xx)
            u1,u2 = u[j-1], u[j]
            v1,v2 = v[j-1], v[j]
            t = (xx-u1)/(u2-u1)
            out[k] = (1-t)*v1 + t*v2
        end
    end
    out
end

"""
    corrected_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat;
                           d,
                           nbins::Int=40,
                           min_points_per_bin::Int=6,
                           min_w_span::Real=0.10,
                           t_min_cut = nothing,
                           bounds = (Kc_min=minimum(K_vals), Kc_max=maximum(K_vals),
                                     nu_min=0.3, nu_max=5.0,
                                     y_min=0.2,  y_max=5.0),
                           verbose::Bool=true)

Corrected finite-time scaling collapse with one irrelevant variable:

    Λ(K,t) = F0(u) + w F1(u),
    u = (K - Kc) * t^(1/(d*ν)),  w = t^(-y/d),

estimated by binning in `u` and weighted linear fits in `w`, while minimizing a global χ²
over (Kc, ν, y).

Returns:
  - best::NamedTuple(Kc, ν, y, χ2)
  - F0, F1 :: interpolants in u (callable on vectors)
  - u, w   :: matrices (size M×N) for best params
  - Y, Yerr:: ln Λ and its errors (M×N)
  - Ycorr  :: corrected data: ln Λ - w·F1(u) (M×N)
"""
function corrected_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat;
                                d,
                                nbins::Int=40,
                                min_points_per_bin::Int=6,
                                min_w_span::Real=0.10,
                                t_min_cut = nothing,
                                bounds = (Kc_min=minimum(K_vals), Kc_max=maximum(K_vals),
                                          nu_min=0.3, nu_max=5.0,
                                          y_min=0.2,  y_max=5.0),
                                verbose::Bool=true)

    M, N = size(p2_mat)
    # Observable and log form (consistent with your other functions)
    Λ   = p2_mat ./ (t_vals' .^ (2/d))
    Y   = log.(Λ)
    σΛ  = p2_err_mat ./ (t_vals' .^ (2/d))
    Yerr = σΛ ./ Λ           # Δ lnΛ ≈ ΔΛ / Λ
    W    = 1.0 ./ max.(Yerr.^2, 1e-16)

    # Optional time cut to suppress very short times
    tidx_keep = t_min_cut === nothing ? collect(1:N) : findall(t_vals .≥ t_min_cut)

    # Helpers: (Kc, ν, y) -> u, w
    u_w_from = function(p::AbstractVector{<:Real})
        Kc, ν, y = p
        u = [(K - Kc) * (t^(1/(d*ν))) for K in K_vals, t in t_vals]    # M×N
        w = [t^(-y/d) for t in t_vals]                                 # 1×N (we'll broadcast)
        return u, w
    end

    # Given u, w, estimate F0(u), F1(u) binwise via weighted LS in w
    function estimate_F0F1(u::Matrix{Float64}, w_vec::Vector{Float64})
        keep = vec([j in tidx_keep for i in 1:M, j in 1:N])
        u_flat = vec(u)
        Y_flat = vec(Y)
        W_flat = vec(W)
        w_flat = vcat([fill(w_vec[j], M) for j in 1:N]...)

        u_use = u_flat[keep]
        Y_use = Y_flat[keep]
        w_use = w_flat[keep]
        WT_use = W_flat[keep]

        order = sortperm(u_use)
        u_sorted  = u_use[order]
        Y_sorted  = Y_use[order]
        w_sorted  = w_use[order]
        WT_sorted = WT_use[order]

        # adaptive number of bins based on available points
        nb = min(nbins, max(5, Int(floor(length(u_sorted) / min_points_per_bin))))
        edges = [round(Int, (k-1)*length(u_sorted)/nb)+1 for k in 1:nb]
        push!(edges, length(u_sorted))

        u_cent, A, B = Float64[], Float64[], Float64[]
        for b in 1:nb
            lo, hi = edges[b], edges[b+1]
            if hi - lo + 1 < min_points_per_bin; continue; end
            ub = u_sorted[lo:hi]
            yb = Y_sorted[lo:hi]
            wb = w_sorted[lo:hi]
            Wb = WT_sorted[lo:hi]

            # require spread in w to identify the slope robustly
            span = (maximum(wb) - minimum(wb)) / max(abs(mean(wb)), 1e-12)
            if span < min_w_span; continue; end

            a, b̂, _, ok = _wlinfit(wb, yb; w=Wb)   # y = a + b*w
            if ok && isfinite(a) && isfinite(b̂)
                push!(u_cent, mean(ub))
                push!(A, a)   # F0(u_b)
                push!(B, b̂)  # F1(u_b)
            end
        end

        if length(u_cent) < 6
            return nothing, nothing
        end
        return _LinInterp(u_cent, A), _LinInterp(u_cent, B)
    end

    # χ² objective for (Kc, ν, y)
    function chi2_with_functions(p)
        Kc, ν, y = p
        if !(isfinite(Kc) && ν>0 && y>0); return (Inf, nothing, nothing); end
        u, w = u_w_from(p)
        F0, F1 = estimate_F0F1(u, [t^(-y/d) for t in t_vals])
        if F0 === nothing; return (Inf, nothing, nothing); end
        u_flat = vec(u)
        w_flat = vcat([fill(t_vals[j]^(-y/d), M) for j in 1:N]...)
        Y_model = F0(u_flat) .+ w_flat .* F1(u_flat)
        resid = vec(Y) .- Y_model
        χ2 = sum(W .* reshape(resid.^2, M, N))
        return (χ2, F0, F1)
    end

    # Initial guesses & bounds
    Kc0 = median(K_vals); ν0 = 1.6; y0 = 1.0
    p0 = [Kc0, ν0, y0]
    lower = [bounds.Kc_min, bounds.nu_min, bounds.y_min]
    upper = [bounds.Kc_max, bounds.nu_max, bounds.y_max]

    obj(p) = chi2_with_functions(p)[1]
    res = optimize(obj, lower, upper, p0, Fminbox(NelderMead()), Optim.Options(iterations=400))
    p⋆ = Optim.minimizer(res)
    χ2⋆, F0⋆, F1⋆ = chi2_with_functions(p⋆)

    if verbose
        println("Corrected collapse:")
        println("  Kc = $(round(p⋆[1], digits=6))")
        println("   ν = $(round(p⋆[2], digits=4))")
        println("   y = $(round(p⋆[3], digits=3))")
        println("  χ² = $(round(χ2⋆, digits=3))")
    end

    # Build outputs at optimum
    u⋆, w⋆ = u_w_from(p⋆)
    u_flat = vec(u⋆)
    F1vals = reshape(F1⋆(u_flat), M, N)
    w_mat  = reshape([t^(-p⋆[3]/d) for _ in 1:M for t in t_vals], M, N)
    Ycorr  = Y .- w_mat .* F1vals

    best = (Kc=p⋆[1], ν=p⋆[2], y=p⋆[3], χ2=χ2⋆)
    return best, F0⋆, F1⋆, u⋆, w⋆, Y, Yerr, Ycorr
end

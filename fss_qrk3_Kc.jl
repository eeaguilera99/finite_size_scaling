include("fss_qkr3_timescaling.jl")  # for finite_time_scaling

using LsqFit
using Plots

"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi; exclude_tol_frac=0.02, plotshow=true)
    # Model function
    model(p, K) = p[1] .+ p[2] .* abs.(K .- p[4]).^(-p[3])   # ξ₀ + A|K−Kc|^{−ν}
    names = ["ξ₀", "A", "ν", "Kc"]

    # Initial guess
    ξ0₀ = minimum(xi)*0.5
    A₀  = maximum(xi)
    ν₀  = 1.5
    Kc₀ = K_vals[argmax(xi)]  # where ξ is largest
    p0 = [ξ0₀, A₀, ν₀, Kc₀]

    # Mask out values too close to trial Kc₀
    ΔK = maximum(K_vals) - minimum(K_vals)
    mask = abs.(K_vals .- Kc₀) .> exclude_tol_frac*ΔK
    Kfit, ξfit = K_vals[mask], xi[mask]

    # Perform nonlinear least squares fit
    fit = curve_fit(model, Kfit, ξfit, p0)
    pbest = coef(fit)
    covar = estimate_covar(fit)
    perr  = sqrt.(diag(covar))

    # Unpack results
    ξ0, A, ν, Kc = pbest
    err_ξ0, err_A, err_ν, err_Kc = perr

    if plotshow
        # Plot raw data
        plt_raw = plot(K_vals, xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)", title="Raw ξ(K) data", label="data")
        display(plt_raw)

        # Plot fit
        Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
        ξfit_model = model(pbest, Kgrid)
        plt_fit = plot(K_vals, xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="ξ(K) with offset fit (LsqFit.jl)", label="data")
        plot!(plt_fit, Kgrid, ξfit_model, lw=2, label="fit (ν ≈ $(round(ν,digits=3)))")
        vline!(plt_fit, [Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt_fit)
    end

    return (ξ0=ξ0, A=A, ν=ν, Kc=Kc,
            err_ξ0=err_ξ0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc)
end





dim = 3  # spatial dimension

# Perform collapse
res, shifts, X, Y, Yerr = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)

# ===== Example usage =====
 xi = exp.(shifts)
 results = fit_xi_offset_LsqFit(K_vals, xi; plotshow=true)
 println(results)



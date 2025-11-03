include("fss_qkr3_timescaling.jl")  # for finite_time_scaling

using LsqFit
using Plots

"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi; exclude_tol_frac=0.02, plotshow=true)
    xi=(1)./xi
    # Model function
    model(K, p) = p[1] .+ p[2] .* abs.(K .- p[4]).^(p[3])   #1/ξ(K) = β₀ + A|K−Kc|^{−ν}
    names = ["β₀", "A", "ν", "Kc"]

    # Initial guess
    β₀₀ = minimum(xi)*0.5
    A₀  = maximum(xi)
    ν₀  = 1.5
    Kc₀ = K_vals[argmax(xi)]  # where ξ is largest
    p0 = [β₀₀, A₀, ν₀, Kc₀]

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
    β0, A, ν, Kc = pbest
    err_β0, err_A, err_ν, err_Kc = perr

    if plotshow
        # Plot raw data
        plt_raw = plot(K_vals, (1)./xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)", title="Raw ξ(K) data", label="data")
        display(plt_raw)
        
        # Plot fit
        Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
        ξfit_model = (1)./model(Kgrid, pbest)
        plt_fit = plot(K_vals, (1)./xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="ξ(K) with offset, Kc ≈ $(round(best.Kc,digits=5))",
                       label="data")
        Kgrid = range(Kmin, Kmax, length=400)
        ξfit = best.ξ0 .+ best.A .* abs.(Kgrid .- best.Kc).^(-best.ν)
        plot!(plt_fit, Kgrid, ξfit, lw=2,
              label="fit (ν ≈ $(round(best.ν,digits=3)))", ylims=(minimum(xi)*0.8, 4))
        vline!(plt_fit, [best.Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt_fit)
    end

    return (β0=β0, A=A, ν=ν, Kc=Kc,
            err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc)
end


dim = 1  # spatial dimension

# Perform collapse
res, shifts, X, Y, Yerr = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)
# ===== Example usage =====
xi = exp.(shifts)
results = fit_xi_offset_LsqFit(K_vals, xi; plotshow=true)
#println(results)

println("\n===== Critical fit results with offset and error bars =====")
println("Kc  ≈ $(results.Kc)  ± $(results.err_Kc)")
println("ν   ≈ $(results.ν)   ± $(results.err_ν)")
println("A   ≈ $(results.A)   ± $(results.err_A)")
println("β_0  ≈ $(results.β0)  ± $(results.err_β0)")
println("Fit quality: ", s_rel)


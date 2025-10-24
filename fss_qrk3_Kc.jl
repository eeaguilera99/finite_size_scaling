include("fss_qkr3_collapse.jl")  # for finite_time_scaling

using Statistics
using Optim
using Plots
using Random

"""
    fit_xi_offset_vsK(K_vals, xi; ngrid=300, exclude_tol_frac=0.02,
                      nboot=200, rng=Random.GLOBAL_RNG,
                      plotshow=true)

Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} and plot ξ vs K.

- Uses grid search on Kc, nonlinear optimization for A, ν, ξ0.
- Returns best-fit parameters + bootstrap error estimates.

Arguments
---------
- K_vals :: Vector of K values
- xi     :: Vector of ξ(K) values
- ngrid  :: Number of Kc grid points (default = 300)
- exclude_tol_frac :: Fraction of K-range to exclude near Kc
- nboot  :: Number of bootstrap resamples for error bars
- rng    :: Random number generator
- plotshow :: whether to make plots

Returns NamedTuple with:
  Kc, ν, A, ξ0, sse,
  err_Kc, err_ν, err_A, err_ξ0
"""

# Perform collapse
res, shifts, X, Y, Yerr = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim)

# Plot before collapse
plt1 = plot(title="Raw data (before shifts)",
    xlabel="ln(t^(-1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label="K=$K")
end
display(plt1)

# Plot after collapse
plt2 = plot(title="Data collapse (after optimal shifts)",
    xlabel="ln(ξ/t^(1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label="K=$K")
end
display(plt2)

#println("Optimal shifts ln ξ(K): ", shifts)
println("Fit quality: ", res.minimum)

function fit_xi_offset_vsK(K_vals, xi; ngrid=300, exclude_tol_frac=0.02,
                           nboot=200, rng=Random.GLOBAL_RNG,
                           plotshow=true)

    Kmin, Kmax = minimum(K_vals), maximum(K_vals)
    ΔK = Kmax - Kmin
    Kc_grid = range(Kmin + 0.05ΔK, Kmax - 0.05ΔK, length=ngrid)

    # SSE for given params
    function sse_offset(Kc, A, ν, ξ0)
        dK = abs.(K_vals .- Kc)
        mask = (dK .> exclude_tol_frac*ΔK) .& (xi .> ξ0)
        if count(mask) < 3
            return 1e12
        end
        pred = ξ0 .+ A .* dK.^(-ν)
        return sum((xi[mask] .- pred[mask]).^2)
    end

    # Best-fit container
    best = (sse = Inf, Kc = NaN, A = NaN, ν = NaN, ξ0 = NaN)

    # Grid search over Kc
    for Kc in Kc_grid
        # crude initial guesses
        A0 = maximum(xi)
        ν0 = 1.5
        ξ00 = minimum(xi) * 0.5

        res = optimize(p -> sse_offset(Kc, exp(p[1]), exp(p[2]), exp(p[3])),
                       [log(A0), log(ν0), log(ξ00)],
                            NelderMead(), Optim.Options(iterations=2000))
        A = exp(Optim.minimizer(res)[1])
        ν = exp(Optim.minimizer(res)[2])
        ξ0 = exp(Optim.minimizer(res)[3])
        sse = Optim.minimum(res)

        if sse < best.sse && ν > 0
            best = (sse=sse, Kc=Kc, A=A, ν=ν, ξ0=ξ0)
        end
    end

    # === Bootstrap error estimates ===
    boot_params = zeros(nboot, 4)
    for b in 1:nboot
        idx = rand(rng, 1:length(K_vals), length(K_vals))  # resample with replacement
        Kb, xib = K_vals[idx], xi[idx]

        # inner fit
        local_best = (sse = Inf, Kc = NaN, A = NaN, ν = NaN, ξ0 = NaN)
        for Kc in Kc_grid
            A0, ν0, ξ00 = maximum(xib), 1.5, minimum(xib)*0.5
            res = optimize(p -> begin
                                dK = abs.(Kb .- Kc)
                                mask = (dK .> exclude_tol_frac*ΔK) .& (xib .> exp(p[3]))
                                if count(mask) < 3
                                    return 1e12
                                end
                                pred = exp(p[3]) .+ exp(p[1]) .* dK.^(-exp(p[2]))
                                return sum((xib[mask] .- pred[mask]).^2)
                            end,
                            [log(A0), log(ν0), log(ξ00)],
                            NelderMead(), Optim.Options(iterations=1000))
            A = exp(Optim.minimizer(res)[1])
            ν = exp(Optim.minimizer(res)[2])
            ξ0 = exp(Optim.minimizer(res)[3])
            sse = Optim.minimum(res)
            if sse < local_best.sse && ν > 0
                local_best = (sse=sse, Kc=Kc, A=A, ν=ν, ξ0=ξ0)
            end
        end
        boot_params[b,:] .= [local_best.Kc, local_best.ν, local_best.A, local_best.ξ0]
    end

    # Compute bootstrap std devs
    err_Kc = std(boot_params[:,1])
    err_ν  = std(boot_params[:,2])
    err_A  = std(boot_params[:,3])
    err_ξ0 = std(boot_params[:,4])

    if plotshow && isfinite(best.Kc)
        # 1) Plot raw ξ(K) vs K (no fit)
        plt_raw = plot(K_vals, xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="Raw ξ(K) data", label="data")
        display(plt_raw)

        # 2) Plot fit with divergence
        plt_fit = plot(K_vals, xi, seriestype=:scatter, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="ξ(K) with offset, Kc ≈ $(round(best.Kc,digits=5))",
                       label="data")
        Kgrid = range(Kmin, Kmax, length=400)
        ξfit = best.ξ0 .+ best.A .* abs.(Kgrid .- best.Kc).^(-best.ν)
        plot!(plt_fit, Kgrid, ξfit, lw=2,
              label="fit (ν ≈ $(round(best.ν,digits=3)))")
        vline!(plt_fit, [best.Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt_fit)
    end

    return (Kc = best.Kc, ν = best.ν, A = best.A, ξ0 = best.ξ0,
            err_Kc=err_Kc, err_ν=err_ν, err_A=err_A, err_ξ0=err_ξ0,
            sse = best.sse)
end

# ===== Example usage (after you have shifts) =====
xi = exp.(shifts)
res = fit_xi_offset_vsK(K_vals, xi; ngrid=400, exclude_tol_frac=0.02, nboot=200, plotshow=true)

println("\n===== Critical fit results with offset and error bars =====")
println("Kc  ≈ $(res.Kc)  ± $(res.err_Kc)")
println("ν   ≈ $(res.ν)   ± $(res.err_ν)")
println("A   ≈ $(res.A)   ± $(res.err_A)")
println("ξ_0  ≈ $(res.ξ0)  ± $(res.err_ξ0)")



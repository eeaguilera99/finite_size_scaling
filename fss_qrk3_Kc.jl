include("fss_qkr3_collapse.jl")  # for finite_time_scaling

using Statistics
using Optim
using Plots

"""
    fit_xi_offset_vsK(K_vals, xi; ngrid=300, exclude_tol_frac=0.02, plotshow=true)

Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} and plot ξ vs K.

- K_vals :: Vector of kick strengths
- xi     :: Vector of extracted ξ(K) = exp.(shifts)
- ngrid  :: number of trial Kc values (grid search)
- exclude_tol_frac :: fraction of K-range excluded near Kc
- plotshow :: whether to plot results

Returns: NamedTuple with best-fit parameters (Kc, ν, A, ξ0)
"""
function fit_xi_offset_vsK(K_vals, xi; ngrid=300, exclude_tol_frac=0.02, plotshow=true)
    Kmin, Kmax = minimum(K_vals), maximum(K_vals)
    ΔK = Kmax - Kmin
    Kc_grid = range(Kmin + 0.05ΔK, Kmax - 0.05ΔK, length=ngrid)

    best = (sse = Inf, Kc = NaN, A = NaN, ν = NaN, ξ0 = NaN) #Initialize the “best fit” record with dummy values (infinite error).

    # SSE for given params
    function sse_offset(Kc, A, ν, ξ0)#Define an inner function to compute the sum of squared errors (SSE).
        dK = abs.(K_vals .- Kc)
        mask = (dK .> exclude_tol_frac*ΔK) .& (xi .> ξ0)
        if count(mask) < 3
            return 1e12
        end
        pred = ξ0 .+ A .* dK.^(-ν)
        return sum((xi[mask] .- pred[mask]).^2)#Return squared error between data and prediction
    end

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

    if plotshow && isfinite(best.Kc)
        plt = plot(K_vals, xi, seriestype=:scatter, ms=6,
                   xlabel="K", ylabel="ξ(K)",
                   title="ξ(K) with offset, diverging at Kc ≈ $(round(best.Kc,digits=5))",
                   label="data")
        Kgrid = range(Kmin, Kmax, length=400)
        ξfit = best.ξ0 .+ best.A .* abs.(Kgrid .- best.Kc).^(-best.ν)
        plot!(plt, Kgrid, ξfit, lw=2, label="fit (ν ≈ $(round(best.ν,digits=3)))")
        vline!(plt, [best.Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt)
    end

    return best
end

# ===== Example usage (after collapse) =====
xi = exp.(shifts)
fitres = fit_xi_offset_vsK(K_vals, xi; ngrid=400, exclude_tol_frac=0.02, plotshow=true)

println("\n===== Critical fit results with offset =====")
println("Kc  ≈ ", fitres.Kc)
println("ν   ≈ ", fitres.ν)
println("A   ≈ ", fitres.A)
println("ξ0  ≈ ", fitres.ξ0)

#=
# ===================== Fit ξ(K) ≈ A |K - Kc|^{-ν} =====================
function fit_xi_powerlaw_vsK(K_vals, xi; ngrid=400, exclude_tol_frac=0.02, plotshow=true)
    Kmin, Kmax = minimum(K_vals), maximum(K_vals)
    ΔK = Kmax - Kmin
    Kc_grid = range(Kmin + 0.05ΔK, Kmax - 0.05ΔK, length=ngrid)#define grid of Kc values to search over

    best = (sse = Inf, Kc = NaN, A = NaN, ν = NaN, used_idx = falses(length(K_vals)))

    # Define simple linear regression in log–log: log(ξ) = log(A) - ν log|K - Kc|
    linfit(x, y) = begin
        mx, my = mean(x), mean(y)
        Sxx = sum((x .- mx).^2)
        Sxy = sum((x .- mx) .* (y .- my))
        slope = Sxy / Sxx
        intercept = my - slope*mx
        yhat = intercept .+ slope .* x
        sse = sum((y .- yhat).^2)
        (intercept, slope, sse) # y = intercept + slope * x
    end
    #Perform linear regression for each Kc in grid
    for Kc in Kc_grid
        tol = exclude_tol_frac * ΔK
        mask = abs.(K_vals .- Kc) .> tol
        if count(mask) < 3
            continue
        end
        x = log.(abs.(K_vals[mask] .- Kc))
        y = log.(xi[mask])
        if std(x) < 1e-12
            continue
        end
        c0, c1, sse = linfit(x, y)
        A, ν = exp(c0), -c1
        if isfinite(sse) && sse < best.sse && ν > 0
            best = (sse=sse, Kc=Kc, A=A, ν=ν, used_idx=mask)
        end
    end

    if plotshow && isfinite(best.Kc)
        # Plot ξ(K) vs K with divergence at Kc
        plt = plot(K_vals, xi, seriestype=:scatter, ms=6,
                   xlabel="K", ylabel="ξ(K)",
                   title="ξ(K) diverging at Kc ≈ $(round(best.Kc,digits=4))")
        # Overlay fit curve
        Kfit = range(Kmin, Kmax, length=300)
        ξfit = best.A .* abs.(Kfit .- best.Kc).^(-best.ν)
        plot!(plt, Kfit, ξfit, lw=2, label="fit (ν ≈ $(round(best.ν,digits=3)))")
        vline!(plt, [best.Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt)
    end

    return (Kc = best.Kc, ν = best.ν, A = best.A, sse = best.sse, used_idx=best.used_idx)
end

# ===== Example usage (after you have `shifts`) =====
xi = exp.(shifts)# Convert back from log ξ to ξ
fitres = fit_xi_powerlaw_vsK(K_vals, xi; ngrid=500, exclude_tol_frac=0.02, plotshow=true)

println("\n===== Critical fit results =====")
println("Kc ≈ ", fitres.Kc)
println("ν  ≈ ", fitres.ν)
println("A  ≈ ", fitres.A)
println("SSE (log-space) = ", fitres.sse)=#
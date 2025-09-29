include("fss_qkr3_collapse.jl")  # for finite_time_scaling

using Statistics
using Optim
using Plots
using Random

"""
    fit_xi_offset_vsK(K_vals, xi, xi_err; ngrid=300, exclude_tol_frac=0.02,
                      nboot=200, plotshow=true)

Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using grid search + nonlinear LSQ,
with bootstrap error bars and optional plots.

Inputs
------
- K_vals :: Vector{Float64}  : kick strengths
- xi     :: Vector{Float64}  : scaling lengths exp.(shifts)
- xi_err :: Vector{Float64}  : error bars for xi
- ngrid  :: Int (default=300): number of trial Kc values
- exclude_tol_frac :: Float64 (default=0.02): fraction of K-range excluded near Kc
- nboot  :: Int (default=200): number of bootstrap resamples
- plotshow :: Bool : whether to make plots

Returns
-------
NamedTuple with fields:
  Kc, ν, A, ξ0, sse, Kc_err, ν_err, A_err, ξ0_err
"""
function fit_xi_offset_vsK(K_vals, xi, xi_err;
                           ngrid=300, exclude_tol_frac=0.02,
                           nboot=200, plotshow=true)

    Kmin, Kmax = minimum(K_vals), maximum(K_vals)
    ΔK = Kmax - Kmin
    Kc_grid = range(Kmin + 0.05ΔK, Kmax - 0.05ΔK, length=ngrid)

    # --- Helper: SSE for given params ---
    function sse_offset(Kc, A, ν, ξ0; mask=nothing)
        dK = abs.(K_vals .- Kc)
        local_mask = (dK .> exclude_tol_frac*ΔK) .& (xi .> ξ0)
        if mask !== nothing
            local_mask .&= mask
        end
        if count(local_mask) < 3
            return 1e12
        end
        pred = ξ0 .+ A .* dK.^(-ν)
        return sum(((xi[local_mask] .- pred[local_mask])./xi_err[local_mask]).^2)
    end

    # --- Inner optimizer for given Kc ---
    function fit_for_Kc(Kc)
        A0, ν0, ξ00 = maximum(xi), 1.5, minimum(xi)*0.5
        res = optimize(p -> sse_offset(Kc, exp(p[1]), exp(p[2]), exp(p[3])),
                       [log(A0), log(ν0), log(ξ00)],
                       NelderMead(), Optim.Options(iterations=2000))
        A = exp(Optim.minimizer(res)[1])
        ν = exp(Optim.minimizer(res)[2])
        ξ0 = exp(Optim.minimizer(res)[3])
        sse = Optim.minimum(res)
        return (A, ν, ξ0, sse)
    end

    # --- Grid search over Kc ---
    best = (sse=Inf, Kc=NaN, A=NaN, ν=NaN, ξ0=NaN)
    for Kc in Kc_grid
        A, ν, ξ0, sse = fit_for_Kc(Kc)
        if sse < best.sse && ν > 0
            best = (sse=sse, Kc=Kc, A=A, ν=ν, ξ0=ξ0)
        end
    end

    # --- Bootstrap to get error bars ---
    boot_params = zeros(nboot,4)  # cols = [Kc, ν, A, ξ0]
    rng = MersenneTwister(1234)
    for b in 1:nboot
        # resample indices with replacement
        idx = rand(rng, 1:length(K_vals), length(K_vals))
        Kb, ξb, errb = K_vals[idx], xi[idx], xi_err[idx]
        # quick refit
        best_b = (sse=Inf, Kc=best.Kc, A=best.A, ν=best.ν, ξ0=best.ξ0)
        for Kc in Kc_grid
            function sse_b(Kc, A, ν, ξ0)
                dK = abs.(Kb .- Kc)
                mask = (dK .> exclude_tol_frac*ΔK) .& (ξb .> ξ0)
                if count(mask) < 3
                    return 1e12
                end
                pred = ξ0 .+ A .* dK.^(-ν)
                return sum(((ξb[mask] .- pred[mask])./errb[mask]).^2)
            end
            res = optimize(p -> sse_b(Kc, exp(p[1]), exp(p[2]), exp(p[3])),
                           [log(best.A), log(best.ν), log(best.ξ0)],
                           NelderMead(), Optim.Options(iterations=1000))
            A = exp(Optim.minimizer(res)[1])
            ν = exp(Optim.minimizer(res)[2])
            ξ0 = exp(Optim.minimizer(res)[3])
            sse = Optim.minimum(res)
            if sse < best_b.sse && ν > 0
                best_b = (sse=sse, Kc=Kc, A=A, ν=ν, ξ0=ξ0)
            end
        end
        boot_params[b,:] = [best_b.Kc, best_b.ν, best_b.A, best_b.ξ0]
    end

    Kc_err, ν_err, A_err, ξ0_err = std(boot_params, dims=1)[:]

    # --- Plots ---
    if plotshow
        # Scatter-only with error bars
        plt1 = scatter(K_vals, xi, yerr=xi_err, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="Raw ξ(K) data with error bars",
                       label="data")
        display(plt1)

        # Fit overlay
        plt2 = scatter(K_vals, xi, yerr=xi_err, ms=6,
                       xlabel="K", ylabel="ξ(K)",
                       title="Fit with offset: Kc ≈ $(round(best.Kc,digits=4))",
                       label="data")
        Kgrid = range(Kmin, Kmax, length=400)
        ξfit = best.ξ0 .+ best.A .* abs.(Kgrid .- best.Kc).^(-best.ν)
        plot!(plt2, Kgrid, ξfit, lw=2,
              label="fit (ν ≈ $(round(best.ν,digits=3)))")
        vline!(plt2, [best.Kc], linestyle=:dash, color=:red, label="Kc")
        display(plt2)
    end

    return (Kc=best.Kc, ν=best.ν, A=best.A, ξ0=best.ξ0, sse=best.sse,
            Kc_err=Kc_err, ν_err=ν_err, A_err=A_err, ξ0_err=ξ0_err)
end


# ===== Example usage (after collapse) =====
xi = exp.(shifts)
xi_err = xi .* 0.05  # <-- put your actual errors here!

fitres = fit_xi_offset_vsK(K_vals, xi, xi_err; ngrid=400, nboot=200, plotshow=true)

println("\n===== Critical fit with offset + error bars =====")
println("Kc   = $(fitres.Kc) ± $(fitres.Kc_err)")
println("ν    = $(fitres.ν) ± $(fitres.ν_err)")
println("A    = $(fitres.A) ± $(fitres.A_err)")
println("ξ0   = $(fitres.ξ0) ± $(fitres.ξ0_err)")

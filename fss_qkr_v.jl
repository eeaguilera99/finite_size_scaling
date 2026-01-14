include("fss_qkr3_timescaling.jl")  # for finite_time_scaling

"""
finite_time_linear_scaling(K_vals, t_vals, p2_mat; Kc, ΔKfit=0.3, plotshow=true)

Implements the Lemarié finite-time-scaling method:

1. Compute Λ = <p²>/t^(2/3).
2. For each t, fit ln Λ ≈ ln Λc + s(t)*(K−Kc)
   using only points |K−Kc| < ΔKfit.
3. Plot ln Λ(K) vs K (figure 13-like).
4. Plot ln|s(t)| vs ln t and extract ν from the slope (figure 14-like).

Returns a NamedTuple with ν, slope, intercept, and the vectors of s(t).
"""
function finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim; Kc, ΔKfit=0.1, plotshow=true, n_kicks_i=1, n_kicks_f=0)

    
    #filter Nkicks range
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]

    K_vals = collect(K_vals)
    t_vals = collect(t_vals)
    M, N = size(p2_mat)

    # 1️⃣ Build Λ(K,t)
    Λ = p2_mat ./ (t_vals' .^ (2/dim))
    Λ_err = p2_err_mat ./ (t_vals' .^ (2/dim))
    lnΛ = log.(Λ)
    ln_Λ_err = Λ_err ./ Λ   # propagate errors: Δ(ln Λ) ≈ ΔΛ / Λ

    # 2️⃣ Fit lnΛ ≈ lnΛc + s(t)*(K−Kc) near Kc
    #create similar dimension arrays
    s_vals = similar(t_vals)
    s_errs = similar(t_vals)
    lnΛc_vals = similar(t_vals)
    fit_lines = Vector{Tuple{Float64,Float64}}(undef, length(t_vals))  # (intercept,slope)    


    mask_global = abs.(K_vals .- Kc) .<= ΔKfit
    Kfit = K_vals[mask_global]
    if length(Kfit) < 3
        error("Not enough K points near Kc. Increase ΔKfit.")
    end
    X = Kfit .- Kc
    Kfit_interval = indexin(Kfit, K_vals)
    # allocate chi2 storage per time and compute per-time linear fits
    chi2_per_time = zeros(length(t_vals))
    redchi2_per_time = zeros(length(t_vals))
    for (j,t) in enumerate(t_vals)
        y = lnΛ[mask_global, j]
        # linear regression y = a + b*X
        A = hcat(ones(length(X)), X)
        coeffs = A \ y #solve least squares
        yfit = A * coeffs
        resid = y - yfit
        σ2 = sum(resid.^2) / (length(y)-2)
        cov = σ2 * inv(A'A)
        lnΛc_vals[j] = coeffs[1]
        s_vals[j] = coeffs[2]
        s_errs[j] = sqrt(cov[2,2])
        fit_lines[j] = (coeffs[1], coeffs[2])

        # chi-square for this linear fit using lnΛ measurement errors (guard zeros)
        σ_y = copy(ln_Λ_err[mask_global, j])
        pos = σ_y .> 0
        if any(pos)
            σ_y[.!pos] .= maximum(σ_y[pos]) + eps()
        else
            σ_y .= maximum(abs.(y)) + eps()
        end
        chi2_j = sum(((y .- yfit) ./ σ_y).^2)
        dof_j = length(y) - 2
        chi2_per_time[j] = chi2_j
        redchi2_per_time[j] = chi2_j / max(dof_j, 1)
    end

    # 3️⃣ log–log fit of |s(t)| vs ln t (weighted)
    logt = log.(t_vals)
    logs = log.(abs.(s_vals))
    logs_err = s_errs ./ abs.(s_vals)   # σ(ln s) = σ_s / |s|
    A = hcat(ones(length(logt)), logt)

# Weighted linear regression using weights = 1/σ²
    # guard logs_err zeros before building W
    logs_err_safe = copy(logs_err)
    posl = logs_err_safe .> 0
    if any(posl)
        logs_err_safe[.!posl] .= maximum(logs_err_safe[posl]) + eps()
    else
        logs_err_safe .= eps()
    end
    W = Diagonal(1.0 ./ (logs_err_safe.^2))
    covmat = inv(A' * W * A)
    coeff = covmat * (A' * W * logs) # weighted least squares solution
    logs_fit = A * coeff
    slope = coeff[2]
    slope_err = sqrt(covmat[2,2])
    ν = 1 / (dim * slope)
    # propagate error: dν/ds = -1/(dim * s^2)
    ν_err = (1.0 / (dim * slope^2)) * slope_err

    intercept = coeff[1]

    # 4️⃣ Plots
    if plotshow
        #=
        # Fig. 13: lnΛ vs K for several t
        plt1 = plot(title=L"\ln{Λ(K)}"*" for various "*L"t",
                    xlabel=L"K", ylabel=L"\ln{Λ(K)}", legend=:topleft)
        for j in 1:N
            plot!(plt1, K_vals, lnΛ[:,j], label="t=$(round(t_vals[j],digits=3))", lw=1.8)
        end
        
        vline!(plt1, [Kc], color=:red, linestyle=:dash, label="Kc")
        display(plt1)=#

        plt2 = plot(title=latexstring("Linear fits near \$a_s=$(a_s)a_0\$ \$Δκ_{fit}\$=$(ΔKfit)"),
                    xlabel=L"κ", ylabel=L"\ln{Λ(κ_c)}", legend=:topleft)
        for j in 1:N
            a, b = fit_lines[j]
            Kloc = Kfit
            yloc = a .+ b .* (Kloc .- Kc)
            scatter!(plt2, Kloc, lnΛ[Kfit_interval,j], yerror=ln_Λ_err[Kfit_interval,j], label="", lw=1.8)#t=$(round(t_vals[j],digits=3))
            plot!(plt2, Kloc, yloc, lw=2, ls=:dash, label="")
        end
        vline!(plt2, [Kc], color=:red, linestyle=:dash, label=L"\kappa_c="*"$(round(Kc,digits=3))")
        display(plt2)
        #savefig(plt2, "fss_linear_fits_d$(dim)_Kc$(round(Kc,digits=3)).png")
         
        # Fig. 14: ln|s| vs ln t
        plt3 = plot(xlabel=L"\ln{t}", ylabel=L"(\ln{Λ})'(κ_c)",
                    title="Scaling of slopes \$a_s=$(a_s)a_0\$", label="data")
        scatter!(plt3, logt, logs; yerr=logs_err, label="data", ms=6)
        plot!(plt3, logt, logs_fit, lw=2, label="fit ν≈$(round(ν,digits=3))"*" ± "*"$(round(ν_err,digits=3))")
        display(plt3)
        #savefig(plt3, "fss_slope_scaling_d$(dim)_κc$(round(Kc,digits=3)).png")
    end

        return (ν=ν, err_ν=ν_err, slope=slope, slope_err=slope_err, intercept=intercept,
            s_vals=s_vals, s_errs=s_errs,
            lnΛc_vals=lnΛc_vals, fit_lines=fit_lines,
            chi2_loglog=chi2, redchi2_loglog=reduced_chi2,
            chi2_per_time=chi2_per_time, redchi2_per_time=redchi2_per_time)
end

Kc = 0.926
dim = 3

res = finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim;
                                 Kc = Kc, ΔKfit = 0.5, n_kicks_i=2, n_kicks_f=0)

println("\n===== Linear finite-time-scaling results =====")
println("ν  = $(round(res.ν,digits=4)) ± $(round(res.err_ν,digits=4))")
#println("slope (1/3ν) = $(round(res.slope,digits=5)) ± $(round(res.slope_err,digits=5))")
println("Goodness of the fit")
println("χ2 loglog = $(round(res.chi2_loglog ,digits=4))")
println("χ2 log log red = $(round(res.redchi2_loglog, digits=4))")
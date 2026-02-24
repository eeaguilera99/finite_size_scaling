using Plots
using LaTeXStrings

#plotting function for raw and collapsed data
function perform_collapse(K_vals, X, Y, Yerr, shifts, s_rel; raw=false, d=dim, data_type="", plot_label_b=false)
    if plot_label_b == true
        plot_label = "K="*string(round(K, digits=3))
    else
        plot_label = ""
    end
    if raw == true
        plt1 = plot(title=latexstring("Raw data $(data_type) \$d=$(d)\$, \$a_s=$(a_s)a_0\$"),
            xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
        for (i,K) in enumerate(K_vals)
            plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label=plot_label)
        end
        display(plt1)
        #savefig(plt1, "fss_raw_data_d$(dim).png")
    end

    # Plot after collapse
    plt2 = plot(title=latexstring("Data collapse $(data_type) \$d=$(d)\$, \$a_s=$(a_s)a_0\$"),
        xlabel=latexstring("\$\\ln(\\xi/N^{1/d})\$"), ylabel=latexstring("\$\\ln(\\Lambda)\$"))
    for (i,K) in enumerate(K_vals)
        plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label=plot_label)
    end
    display(plt2)   
    #savefig(plt2, "fss_collapsed_data_d$(dim).png")
    println("Fit quality $(data_type): ", s_rel)
end

#critical Kc analysis from collapse shifts
function perform_Kc_anal(shifts, shifts_err, K_vals; data_type="", d=3, n_k_filter=0, K_guess_index=0)
        xi = exp.(shifts)
        xierr = xi.*shifts_err
        results = fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=n_k_filter, K_val_g=K_guess_index)
        # Plot fit
        Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
        plt_fit1 = plot(K_vals, xi, seriestype=:scatter, ms=6,
                        xlabel=L"κ", ylabel=L"ξ(κ)",
                        title=latexstring("\$$data_type\$, \$κ_c≈ $(round(results.Kc,digits=3))\$, \$d=$(d)\$"),
                        label="data")
        ξfit1 = (1)./(results.β0 .+ results.A .* abs.(Kgrid .- results.Kc).^(abs(results.ν)))
        plot!(plt_fit1, Kgrid, ξfit1, lw=2,
                label="fit (ν ≈ $(round(abs(results.ν),digits=3)))")
        vline!(plt_fit1, [results.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
        display(plt_fit1)
        println("\n===== Critical fit $(data_type) with offset and error bars =====")
        println("κc1  ≈ $(results.Kc)  ± $(results.err_Kc)")
        println("ν1   ≈ $(results.ν)   ± $(results.err_ν)")
        println("A1   ≈ $(results.A)   ± $(results.err_A)")
        println("β_01  ≈ $(results.β0)  ± $(results.err_β0)")
        println("χ²1  = $(results.χ2),  χ²_red1 = $(results.χ2_red)")
end

#slopes analysis
function finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim; data_type ="", Kc, ΔKfit=0.1, plotshow=true, n_kicks_i=2, n_kicks_f=0)

    
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

    # goodness of fit for the ln|s| vs ln t weighted fit
    resid_logs = logs .- logs_fit
    chi2 = sum((resid_logs ./ logs_err_safe).^2)
    dof = length(logs) - 2
    reduced_chi2 = chi2 / max(dof, 1)

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

        plt2 = plot(title=latexstring("Linear fits near $data_type \$a_s=$(a_s)a_0\$ \$Δκ_{fit}\$=$(ΔKfit)"),
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
                    title=latexstring("Scaling of slopes $data_type \$a_s=$(a_s)a_0\$"), label="data")
        scatter!(plt3, logt, logs; yerr=logs_err, label="data", ms=6)
        plot!(plt3, logt, logs_fit, lw=2, label="fit ν≈$(round(ν,digits=3))"*" ± "*"$(round(ν_err,digits=3))")
        display(plt3)
        #savefig(plt3, "fss_slope_scaling_d$(dim)_κc$(round(Kc,digits=3)).png")
    end
    println("\n===== Linear finite-time-scaling results $data_type=====")
    println("ν  = $(round(ν,digits=4)) ± $(round(ν_err,digits=4))")
    println("Goodness of the fit")
    println("χ2 loglog = $(round(chi2 ,digits=4))")
    println("χ2 log log red = $(round(reduced_chi2, digits=4))")
end

#Time evolution at criticality for dimension
function perform_tevol(d, shifts, t_vals, data_mat, data_mat_err; Kc_offset_index=0, data_type="")

    xi1 = exp.(shifts)
    K_c1_i = argmax(xi1) + Kc_offset_index

    #Plot straight lines 
    model(t, p) = p[1] .* t .+ p[2]  # y = m*x + b
    guess1 = [2/d, 0.0]  # Initial guess for [m, b]
    fit1 = curve_fit(model, log.(t_vals), log.(data_mat[K_c1_i, :]), guess1)
    fit_params1 = coef(fit1)
    fit_errs1 = standard_errors(fit1)
    α_err = fit_errs1[1]

    #χ^2 goodness of fit
    residuals1 = log.(data_mat[K_c1_i, :]) .- model(log.(t_vals), fit_params1)
    χ2_1 = sum((residuals1 ./ (data_mat_err[K_c1_i, :]./data_mat[K_c1_i, :])).^2)
    dof1 = length(t_vals) - length(fit_params1)
    println("χ2 =", χ2_1)
    println("χ2 rel =", χ2_1/dof1)

    #Plot
    plt1 = plot(title=latexstring("Time evolution of $(data_type) near \$κ_c\$"),
        xlabel="Time (kicks)", ylabel=latexstring("$(data_type)"))
    plot!(plt1, t_vals, data_mat[K_c1_i, :], xscale=:log10, yscale=:log10, marker=:o, label="κ_c=$(round(K_vals[K_c1_i], digits=3))", ms=3)
    plot!(plt1, t_vals, exp.(model(log.(t_vals),fit_params1)), xscale=:log10, yscale=:log10, label=latexstring("\$α≈$(round(fit_params1[1], digits=3)) ± $(round(α_err, digits=3))\$"))
    display(plt1)
    println("Approximate dimension from fit: d ≈ $(round(2/fit_params1[1], digits=3)) ± $(round(2*α_err/fit_params1[1]^2, digits=3))")
end
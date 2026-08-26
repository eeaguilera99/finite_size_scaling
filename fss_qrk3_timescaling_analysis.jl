
using Plots
using LaTeXStrings
using LsqFit

#plotting function for raw and collapsed data
function perform_collapse(K_vals, X, Y, Yerr, shifts; raw=false, d=dim, data_type="", plot_label_b=false, ploterr=true, save=false)
    if plot_label_b == true
        plot_label = "K="*string(round(K, digits=3))
    else
        plot_label = ""
    end
    if raw == true
        plt1 = plot(title=latexstring("Raw data $(data_type) \$d=$(d)\$, \$a_s=$(a_s)a_0\$"),
            xlabel=latexstring("\$\\ln(N_p^{-1/d})\$"), ylabel=latexstring("\$ \\ln (Λ)\$"))
        for (i,K) in enumerate(K_vals)
            plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label=plot_label)
        end
        if save == true
            savefig(plt1, "C:\\Users\\c7041417\\Documents\\Julia\\QKR\\finite_size_scaling\\plots\\Ex\\fss_raw_data_d=$(d)_a=$(a_s).pdf") 
        end 
        display(plt1)
    end

    # Plot after collapse
    plt2 = plot(title=latexstring("Data collapse $(data_type) \$d=$(d)\$, \$a_s=$(a_s)a_0\$"),
        xlabel=latexstring("\$\\ln(\\xi/N_p^{1/d})\$"), ylabel=latexstring("\$\\ln(\\Lambda)\$"))
    if ploterr == true
        for (i,K) in enumerate(K_vals)
            plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label=plot_label)
        end
    else
        for (i,K) in enumerate(K_vals)
            plot!(plt2, X[:] .+ shifts[i], Y[i,:], marker=:o, label=plot_label)
        end
    end
    if save == true
        savefig(plt2, "C:\\Users\\c7041417\\Documents\\Julia\\QKR\\finite_size_scaling\\plots\\Ex\\fss_collapsed_data_d=$(d)_a=$(a_s).pdf")
    end
    display(plt2)
end

function perform_collapse2(K_vals, X_data, Y_data, Yerr_data, shifts, χ2, χ2_red, pbest, perr, D, a_s; showxi=true)
    Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
    plt5 = plot(title="Scaling by Taylor fit \$d=$(D)\$, \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
    for i in eachindex(K_vals)

        #fixed K
        K = K_vals[i]

        # Horizontal shift for each K
        X_collapse = X_data' .+ shifts[i] 

        # ACTUAL DATA
        Y_actual = Y_data[i, :]
        Y_err_actual = Yerr_data[i, :]

        #fit data
        # Use the fitted model to calculate ln Λ
        #t = tt_vals[j]

        #=logΛ_fit = model(
            [K_vals'; fill(t, length(K_vals))'],
            pbest
        )

        valid_fit = isfinite.(X_collapse) .&
                isfinite.(logΛ_fit)=#

        valid_data = isfinite.(X_collapse) .&
                isfinite.(Y_actual)

        plot!(
            plt5,
            X_collapse[valid_data],
            Y_actual[valid_data],#logΛ_fit[valid_fit],
            yerror=Y_err_actual[valid_data],
            seriestype = :scatter,
            ms = 3,
            label = ""
        )
    end
    display(plt5)
    if showxi
        #plot localiation length ξ(k)
        plt6 = plot(title="Localization length ξ(k) \$d=$(D)\$, \$a_s=$(a_s)a_0\$", xlabel=latexstring("κ"), ylabel=latexstring("ξ(k)"))
        plot!(plt6, K_vals, exp.(shifts), seriestype=:scatter, ms=3, label="ξ(k) data")
        ξplot =  exp.(-pbest[4]*log.(abs.(pbest[1].*(Kgrid .- pbest[3]) .+ pbest[2].*(Kgrid .- pbest[3]).^2)))
        plot!(plt6, Kgrid, ξplot, lw=2, label="ξ(k) fit (ν=$(round(pbest[4], digits=3)))")
        vline!(plt6, [pbest[3]], lw=2, ls=:dash, color=:red, label="Kc=$(round(pbest[3], digits=3))")
        display(plt6)
    end
    

    println(" Collapse \$d=$(D)\$, \$a_s=$(a_s)a_0\$ fitted parameters with errors:")
    println("b1 = $(round(pbest[1], digits=3)) ± $(round(perr[1], digits=3))")
    println("κ_c = $(round(pbest[3], digits=3)) ± $(round(perr[3], digits=3))")
    println("ν = $(round(pbest[4], digits=3)) ± $(round(perr[4], digits=3))")
    println("ξsat = $(round(pbest[6], digits=3)) ± $(round(perr[6], digits=3))")
    println("χ2 = $(round(χ2, digits=3)), χ2_red = $(round(χ2_red, digits=3))")
    
end

function filter_data(X, Y, Y_err)
    Xf = Float64[]
    Yf = Float64[]
    Ef = Float64[]

    @inbounds for i in eachindex(X)
        xi = X[i]
        yi = Y[i]
        ei = Y_err[i]

        if !isnan(xi) && !isnan(yi) && !isnan(ei) && ei > 0
            push!(Xf, xi)
            push!(Yf, yi)
            push!(Ef, ei)
        end
    end

    return Xf, Yf, Ef
end


function perform_collapse_quality(K_vals, X, Y, Y_err, shifts, s_rel, d, a_s, data_type; Kc_offset=0, plotshow=false)

    X_shifted = X .+ shifts
    Kc_index = argmax(shifts)


    Kc_1 = K_vals[Kc_index + Kc_offset]
    Kc_2 = K_vals[Kc_index - Kc_offset]
    K_diff_mask = K_vals .> Kc_1 
    K_loc_mask  = K_vals .< Kc_2

    # Flatten once
    X_diff = vec(X_shifted[K_diff_mask, :])
    Y_diff = vec(Y[K_diff_mask, :])
    E_diff = vec(Y_err[K_diff_mask, :])

    X_loc = vec(X_shifted[K_loc_mask, :])
    Y_loc = vec(Y[K_loc_mask, :])
    E_loc = vec(Y_err[K_loc_mask, :])

    # Filter
    X_diff, Y_diff, E_diff = filter_data(X_diff, Y_diff, E_diff)
    X_loc,  Y_loc,  E_loc  = filter_data(X_loc,  Y_loc,  E_loc)

    isempty(X_diff) && error("No valid diff data")
    isempty(X_loc)  && error("No valid loc data")

    model(x, p) = @. p[1] * x + p[2]

    guess_diff = [- (d - 2), 0.0]
    guess_loc  = [2.0, 0.0]

    fit_diff = curve_fit(model, X_diff, Y_diff, guess_diff)
    fit_loc  = curve_fit(model, X_loc,  Y_loc,  guess_loc)

    coef_diff = coef(fit_diff)
    coef_loc  = coef(fit_loc)

    # Residuals (no temporary model array)
    resid_diff = @. Y_diff - (coef_diff[1] * X_diff + coef_diff[2])
    resid_loc  = @. Y_loc  - (coef_loc[1]  * X_loc  + coef_loc[2])

    μ_diff = mean(Y_diff)
    μ_loc  = mean(Y_loc)

    ss_tot_diff = sum(abs2(y - μ_diff) for y in Y_diff)
    ss_res_diff = sum(abs2, resid_diff)
    R2_diff = 1 - ss_res_diff / ss_tot_diff

    ss_tot_loc = sum(abs2(y - μ_loc) for y in Y_loc)
    ss_res_loc = sum(abs2, resid_loc)
    R2_loc = 1 - ss_res_loc / ss_tot_loc

    χ2_diff = sum(abs2(r / e) for (r,e) in zip(resid_diff, E_diff))
    χ2_loc  = sum(abs2(r / e) for (r,e) in zip(resid_loc,  E_loc))

    χ2_red_diff = χ2_diff / max(length(Y_diff) - length(coef_diff), 1)
    χ2_red_loc  = χ2_loc  / max(length(Y_loc)  - length(coef_loc),  1)

    # slope errors
    slope_diff = coef(fit_diff)[1]
    slope_loc = coef(fit_loc)[1]
    slope_diff_theory = -(d - 2)
    slope_loc_theory  = 2.0

    err_diff = abs((slope_diff - slope_diff_theory)/slope_diff_theory)
    err_loc  = abs((slope_loc  - slope_loc_theory)/slope_loc_theory)

    if plotshow == true
        plt = plot(title=latexstring("Collapse quality analysis $(data_type), \$d=$(d)\$, \$a_s=$(a_s)a_0\$"),
            xlabel=latexstring("\$\\ln(\\xi/N^{1/d})\$"), ylabel=latexstring("\$\\ln(\\Lambda)\$"))
        scatter!(plt, X_diff, Y_diff, seriestype=:scatter, label="diff side")
        plot!(plt, X_diff, model(X_diff, coef(fit_diff)), label="diff fit (slope ≈ $(round(coef(fit_diff)[1], digits=3)))")
        scatter!(plt, X_loc, Y_loc, seriestype=:scatter, label="loc side")
        plot!(plt, X_loc, model(X_loc, coef(fit_loc)), label="loc fit (slope ≈ $(round(coef(fit_loc)[1], digits=3)))")
        display(plt)
    end

    println("\n===== Collapse fit quality $(data_type), d=$(d), a_s=$(a_s) =====")
    println("χ2 red avg = $(round((χ2_red_diff + χ2_red_loc)/2, digits=4))")
    #println("χ2 red diff = $(round(χ2_red_diff, digits=4))")
    #println("χ2 red loc = $(round(χ2_red_loc, digits=4))")
    println("R² avg = $(round((R2_diff + R2_loc)/2, digits=4))")
    #println("R² diff = $(round(R2_diff, digits=4))")
    #println("R² loc = $(round(R2_loc, digits=4))")
    println("Slope err avg = $(round((err_diff + err_loc)/2, digits=4))")  
    #println("Slope error diff = $(round(err_diff, digits=4))")
    #println("Slope error loc = $(round(err_loc, digits=4))")
    println("Collapse tightness quality σ_rel = $(round(s_rel, digits=4))")
end

#critical Kc analysis from collapse shifts
"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=0, K_val_g=0, exclude_tol_frac=0.02)
    #filter points for fit    
    if n_k_filter != 0
        n_k_filter += 1  # account for Julia 1-based indexing    
        xi = xi[n_k_filter:end]
        xierr = xierr[n_k_filter:end]
        K_vals = K_vals[n_k_filter:end]
    end

    u = (1)./xi
    uerr = u.^2 .*xierr
    # Model function
    model(K, p) = abs(p[1]) .+ p[2] .* abs.(K .- p[4]).^(abs(p[3]))   #1/ξ(K) = β₀ + A|K−Kc|^{ν}
    names = ["β₀", "A", "ν", "Kc"]

    # Initial guess
    β₀₀ = minimum(u)*0.5
    A₀  = maximum(u)
    ν₀  = 1
    Kc₀ = K_vals[argmin(u) + K_val_g]  # where ξ is largest
    p0 = [β₀₀, A₀, ν₀, Kc₀]

    # Mask out values too close to trial Kc₀
    ΔK = maximum(K_vals) - minimum(K_vals)
    mask = abs.(K_vals .- Kc₀) .> exclude_tol_frac*ΔK
    Kfit, ufit, uerrfit = K_vals[mask], u[mask], uerr[mask]

    # Perform nonlinear least squares fit
    fit = curve_fit(model, Kfit, ufit, p0)
    pbest = coef(fit)
    covar = estimate_covar(fit)
    perr  = sqrt.(diag(covar))

    # goodness of fit
    residuals = ufit .- model(Kfit, pbest)
    if uerr[1] == 0.0
        χ2 = sum((residuals[2:end] ./ uerrfit[2:end]).^2)
    else
        χ2 = sum((residuals ./ uerrfit).^2)
    end
    dof = length(ufit) - length(pbest)
    χ2_red = χ2 / dof

    # Unpack results
    β0, A, ν, Kc = pbest
    err_β0, err_A, err_ν, err_Kc = perr


    return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
            err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


function perform_Kc_anal(shifts, shifts_err, K_vals; data_type="", d=3, n_k_filter=0, K_guess_index=0, show_xierr=true, save=false)
    xi = exp.(shifts)
    xierr = xi.*shifts_err
    results = fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=n_k_filter, K_val_g=K_guess_index)
    # Plot fit
    Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
    if show_xierr == true 
        plt_fit1 = plot(K_vals, xi, seriestype=:scatter, ms=6,
                        xlabel=L"κ", ylabel=L"ξ(κ)", yerror=xierr,
                        title=latexstring("\$$data_type\$, \$κ_c≈ $(round(results.Kc,digits=3))\$, \$d=$(d)\$"),
                        label="data")
    else
        plt_fit1 = plot(K_vals, xi, seriestype=:scatter, ms=6,
                        xlabel=L"κ", ylabel=L"ξ(κ)",
                        title=latexstring("\$$data_type\$, \$κ_c≈ $(round(results.Kc,digits=3))\$, \$d=$(d)\$"),
                        label="data")
    end
    ξfit1 = (1)./(results.β0 .+ results.A .* abs.(Kgrid .- results.Kc).^(abs(results.ν)))
    plot!(plt_fit1, Kgrid, ξfit1, lw=2,
            label="fit (ν ≈ $(round(abs(results.ν),digits=3)))")
    vline!(plt_fit1, [results.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
    if save == true
        savefig(plt_fit1, "plots\\Ex\\Kc_fit_$(data_type)_d=$(d)_a=$(a_s).pdf")
    end
    display(plt_fit1)

    println("\n===== Critical fit $(data_type) d=$(d), a_s=$(a_s), κ_c_offset=$(K_guess_index) =====")
    println("κc1  ≈ $(round(results.Kc, digits=3))  ± $(round(results.err_Kc, digits=3))")
    println("ν1   ≈ $(round(results.ν, digits=3))   ± $(round(results.err_ν, digits=3))")
    println("A1   ≈ $(round(results.A, digits=3))   ± $(round(results.err_A, digits=3))")
    println("β_01  ≈ $(round(results.β0, digits=3))  ± $(round(results.err_β0, digits=3))")
    println("χ²1  = $(round(results.χ2, digits=3)),  χ²_red1 = $(round(results.χ2_red, digits=3))")
    return results
end

#slopes analysis
function finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim; data_type ="", Kc, ΔKfit=0.1, plotshow=true, n_kicks_i=2, n_kicks_f=0, save=false)

    
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

        plt2 = plot(title=latexstring("Linear fits near $data_type \$a_s=$(a_s)a_0\$ \$d=$(dim)\$"),
                    xlabel=L"κ", ylabel=L"\ln{Λ(κ_c)}", legend=:topleft)
        for j in 1:N
            a, b = fit_lines[j]
            Kloc = Kfit
            yloc = a .+ b .* (Kloc .- Kc)
            scatter!(plt2, Kloc, lnΛ[Kfit_interval,j], yerror=ln_Λ_err[Kfit_interval,j], label="", lw=1.8)#t=$(round(t_vals[j],digits=3))
            plot!(plt2, Kloc, yloc, lw=2, ls=:dash, label="")
        end
        vline!(plt2, [Kc], color=:red, linestyle=:dash, label=L"\kappa_c="*"$(round(Kc,digits=3))")
        
        #
         
        # Fig. 14: ln|s| vs ln t
        plt3 = plot(xlabel=L"\ln{t}", ylabel=L"(\ln{Λ})'(κ_c)",
                    title=latexstring("Scaling of slopes $data_type \$a_s=$(a_s)a_0\$ \$d=$(dim)\$"), label="data")
        scatter!(plt3, logt, logs; yerr=logs_err, label="data", ms=6)
        plot!(plt3, logt, logs_fit, lw=2, label="fit ν≈$(round(ν,digits=3))"*" ± "*"$(round(ν_err,digits=3))")
        

        if save == true
            savefig(plt2, "plots\\Ex\\fss_linear_fits_d$(dim)_Kc$(round(Kc,digits=3)).pdf")
            savefig(plt3, "plots\\Ex\\fss_slope_scaling_d$(dim)_κc$(round(Kc,digits=3)).pdf")
        end
        display(plt2)
        display(plt3)
    end
    println("\n===== Linear finite-time-scaling results $data_type d=$(dim), a_s=$(a_s), Δκfit=$(ΔKfit) =====")
    println("ν  = $(round(ν,digits=3)) ± $(round(ν_err,digits=3))")
    println("Goodness of the fit")
    #println("χ2 loglog = $(round(chi2 ,digits=3))")
    println("χ2 red = $(round(reduced_chi2, digits=3))")
end

#Time evolution at criticality for dimension
function perform_tevol_ln(d, shifts, t_vals, data_mat, data_mat_err; Kc_offset_index=0, data_type="")

    K_c1_i = argmax(shifts) + Kc_offset_index

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
        xlabel="Time (kicks)", ylabel=latexstring("\$E(t, κ_c)\$"))
    plot!(plt1, t_vals, data_mat[K_c1_i, :], xscale=:log10, yscale=:log10, marker=:o, label="κ_c=$(round(K_vals[K_c1_i], digits=3))", ms=3)
    plot!(plt1, t_vals, exp.(model(log.(t_vals),fit_params1)), xscale=:log10, yscale=:log10, label=latexstring("\$α≈$(round(fit_params1[1], digits=3)) ± $(round(α_err, digits=3))\$"))
    display(plt1)
    println("Approximate dimension from fit: d ≈ $(round(2/fit_params1[1], digits=3)) ± $(round(2*α_err/fit_params1[1]^2, digits=3))")
end

function perform_tevol_powerlaw(dime, shifts, t_vals, data_mat, data_mat_err; transient, Kc_offset_index=0, data_type="")
    
    #filter Nkicks range
    t_vals, data_mat, data_mat_err = filter_Nkicks(t_vals, data_mat, data_mat_err; n_kicks_i=transient, n_kicks_f=0)

    K_c1_i = argmax(shifts) + Kc_offset_index
    #fit powerlaw with transient: model(t) = A * (t - t0)^α
    model(t, p) = p[1] .* abs.(t .- p[3]) .^ p[2]  # y = A * (t - t0)^α
    guess1 = [1.0, 2/dime, transient]  # Initial guess for [A, α, t0]
    fit1 = curve_fit(model, t_vals, data_mat[K_c1_i, :], guess1)
    fit_params1 = coef(fit1)
    fit_errs1 = standard_errors(fit1)
    α_err = fit_errs1[2]

    #χ^2 goodness of fit
    residuals1 = data_mat[K_c1_i, :] .- model(t_vals, fit_params1)
    χ2_1 = sum((residuals1 ./ data_mat_err[K_c1_i, :]).^2)
    dof1 = length(t_vals) - length(fit_params1) 
    χ2_red1 = χ2_1/dof1

    #plot
    plt1 = plot(title=latexstring("Time evolution of $(data_type) near \$κ_c\$"),
        xlabel="Time (kicks)", ylabel=latexstring("\$E(t, κ_c)\$"))
    plot!(plt1, t_vals, data_mat[K_c1_i, :].*1e8, marker=:o, label="κ_c=$(round(K_vals[K_c1_i], digits=3))", ms=3)
    plot!(plt1, t_vals, model(t_vals,fit_params1).*1e8, label=latexstring("\$d≈$(round(2/fit_params1[2], digits=3)) ± $(round(2*α_err/fit_params1[2]^2, digits=3))\$"))
    display(plt1)
    println("Approximate dimension from fit: d ≈ $(round(2/fit_params1[2], digits=3)) ± $(round(2*α_err/fit_params1[2]^2, digits=3))")
    println("Estimated transient time t0 ≈ $(round(fit_params1[3], digits=3)) ± $(round(fit_errs1[3], digits=3))")
    println("χ2 red =", χ2_red1)
end

function filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=1, n_kicks_f=0)
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]
    return t_vals, p2_mat, p2_err_mat
end

function perform_tevol_transient(d, shifts, t_vals, data_mat, data_mat_err; Kc_offset_index=0, data_type="", n_kicks_i=1, n_kicks_f=0)

    #filter Nkicks range
    t_vals, data_mat, data_mat_err = filter_Nkicks(t_vals, data_mat, data_mat_err; n_kicks_i=n_kicks_i, n_kicks_f=n_kicks_f)

    K_c1_i = argmax(shifts) + Kc_offset_index

    #fit powerlaw with transient: model(t) = A * (t - t0)^α
    model(t, p) = p[1] .* (t .- p[3]) .^ p[2]  # y = A * (t - t0)^α
    guess1 = [1.0, 2/d, n_kicks_i]  # Initial guess for [A, α, t0]
    fit1 = curve_fit(model, t_vals, data_mat[K_c1_i, :], guess1)
    fit_params1 = coef(fit1)    
    fit_errs1 = standard_errors(fit1)
    α_err = fit_errs1[2]

    #χ^2 goodness of fit
    residuals1 = data_mat[K_c1_i, :] .- model(t_vals,fit_params1)
    χ2_1 = sum((residuals1 ./ data_mat_err[K_c1_i, :]).^2)
    dof1 = length(t_vals) - length(fit_params1) 
    println("χ2 =", χ2_1)
    println("χ2 rel =", χ2_1/dof1)  

    #Plot
    plt1 = plot(title=latexstring("Time evolution of $(data_type) near \$κ_c\$ with transient"),
        xlabel="Time (kicks)", ylabel=latexstring("$(data_type)"))  
    plot!(plt1, t_vals, data_mat[K_c1_i, :], marker=:o, label="κ_c=$(round(K_vals[K_c1_i], digits=3))", ms=3)
    plot!(plt1, t_vals, model(t_vals,fit_params1), label=latexstring("\$α≈$(round(fit_params1[2], digits=3)) ± $(round(α_err, digits=3))\$"))
    display(plt1)
    println("Approximate dimension from fit: d ≈ $(round(2/fit_params1[2], digits=3)) ± $(round(2*α_err/fit_params1[2]^2, digits=3))")
    println("Estimated transient time t0 ≈ $(round(fit_params1[3], digits=3)) ± $(round(fit_errs1[3], digits=3))")
end

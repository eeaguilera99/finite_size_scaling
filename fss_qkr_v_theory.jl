include("fss_qkr3_timescaling_theory.jl")  # for finite_time_scaling

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
function finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim, a_s; type=true, Kc, ΔKfit=0.1, plotshow=true, n_kicks_i=1, n_kicks_f=0)

    
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

    # 3️⃣ log–log fit of |s(t)| vs ln t
    logt = log.(t_vals)
    logs = log.(abs.(s_vals))
    logs_err = s_errs ./ abs.(s_vals)   # σ(ln s) = σ_s / |s|
    A = hcat(ones(length(logt)), logt)

# Weighted linear regression using weights = 1/σ²
    W = Diagonal(1.0 ./ (logs_err.^2))
    covmat = inv(A' * W * A)
    coeff = covmat * (A' * W * logs) # weighted least squares solution
    logs_fit = A * coeff
    slope = coeff[2]
    slope_err = sqrt(covmat[2,2])
    ν = 1 / (dim*slope)
    ν_err = abs(ν^2) * slope_err / dim  # propagate: dν = -(ν²/3)*dslope

    
    intercept = coeff[1]

    # 4️⃣ Plots
    if plotshow
        #= 
        # Fig. 13: lnΛ vs K for several t
        plt1 = plot(title=latexstring("\$ln{Λ(K)}\$ vs \$κ\$, \$a_s=$(a_s)a_0\$"),
                    xlabel=L"κ", ylabel=L"\ln{Λ(K)}", legend=:topleft)
        for j in 1:N
            scatter!(plt1, K_vals, lnΛ[:,j], label="", lw=1.8)
        end
        
        vline!(plt1, [Kc], color=:red, linestyle=:dash, label=latexstring("\$κ_c=$(round(Kc,digits=3))\$"))
        display(plt1)=#
        if type==true
            title = latexstring("Slopes using \$1/n_c^2\$")
        else
            title = latexstring("Slopes using \$E_k\$")
        end

        plt2 = plot(title=latexstring("Linear fits near \$Δκ_{fit}=$(ΔKfit)\$, \$a_s=$(a_s)a_0\$"),
                    xlabel=L"κ", ylabel=L"\ln{Λ(κ)}", legend=:topleft)
        for j in 1:N
            a, b = fit_lines[j]
            Kloc = Kfit
            yloc = a .+ b .* (Kloc .- Kc)
            scatter!(plt2, Kloc, lnΛ[Kfit_interval,j], yerror=ln_Λ_err[Kfit_interval,j], label="", lw=1.8)#t=$(round(t_vals[j],digits=3))
            plot!(plt2, Kloc, yloc, lw=2, ls=:dash, label="")
        end
        vline!(plt2, [Kc], color=:red, linestyle=:dash, label=latexstring("\$κ_c=$(round(Kc,digits=3))\$"))

         
        # Fig. 14: ln|s| vs ln t
        plt3 = plot(xlabel=L"\ln{t}", ylabel=L"(\ln{Λ})'(κ_c)",
                    title=latexstring("Scaling of slopes \$a_s=$(a_s)a_0\$"), label="data")
        scatter!(plt3, logt, logs; yerr=logs_err, label="data", ms=6)
        plot!(plt3, logt, logs_fit, lw=2, label=latexstring("fit \$ν≈$(round(ν,digits=3))±$(round(ν_err,digits=3))\$"))

        display(plot(plt2, plt3, suptitle= title, layout=(1,2), size=(1000,400), 
        bottom_margin=5Plots.mm, left_margin=5Plots.mm))

    end

    return (ν=ν, err_ν=ν_err, slope=slope, slope_err=slope_err, intercept=intercept,
            s_vals=s_vals, s_errs=s_errs,
            lnΛc_vals=lnΛc_vals, fit_lines=fit_lines)
end

Kc_nc = 1.25
Kc_ek = 1.354
dim=3


results1 = finite_time_linear_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=true, loc_amp=6), nc_err_mat, dim, a_s;
                                 Kc = Kc_nc, ΔKfit = 1, n_kicks_i=18, n_kicks_f=0)
results2 = finite_time_linear_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=true, loc_amp=6), p2_err_mat, dim, a_s; type=false,
                                 Kc = Kc_ek, ΔKfit = 1, n_kicks_i=18, n_kicks_f=0)

println("\n===== Linear finite-time-scaling results =====")
println("ν  = $(round(results1.ν,digits=4)) ± $(round(results1.err_ν,digits=4))")

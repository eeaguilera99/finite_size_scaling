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
function finite_time_linear_scaling(K_vals, t_vals, p2_mat, dim; Kc, ΔKfit=0.1, plotshow=true, n_kicks_i=1, n_kicks_f=0)

    
    #filter Nkicks range
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    #p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]

    K_vals = collect(K_vals)
    t_vals = collect(t_vals)
    M, N = size(p2_mat)

    # 1️⃣ Build Λ(K,t)
    Λ = p2_mat ./ (t_vals' .^ (2/dim))
    lnΛ = log.(Λ)

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

    for (j,t) in enumerate(t_vals)
        y = lnΛ[mask_global, j]
        # linear regression y = a + b*X
        A = hcat(ones(length(X)), X)
        coeffs = A \ y
        yfit = A * coeffs
        resid = y - yfit
        σ2 = sum(resid.^2) / (length(y)-2)
        cov = σ2 * inv(A'A)
        lnΛc_vals[j] = coeffs[1]
        s_vals[j] = coeffs[2]
        s_errs[j] = sqrt(cov[2,2])
        fit_lines[j] = (coeffs[1], coeffs[2])
    end

    # 3️⃣ log–log fit of |s(t)| vs t
    logt = log.(t_vals)
    logs = log.(abs.(s_vals))
    A = hcat(ones(length(logt)), logt)
    coeff = A \ logs
    logs_fit = A * coeff
    slope = coeff[2]
    ν = 1 / (3*slope)
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

        plt2 = plot(title="Linear fits near Kc=$(round(Kc,digits=3))",
                    xlabel=L"K", ylabel=L"\ln{Λ(K)}", legend=:topleft)
        for j in 1:N
            a, b = fit_lines[j]
            Kloc = Kfit
            yloc = a .+ b .* (Kloc .- Kc)
            scatter!(plt2, Kloc, lnΛ[:,j], label="", lw=1.8)#t=$(round(t_vals[j],digits=3))
            plot!(plt2, Kloc, yloc, lw=2, ls=:dash, label="")
        end
        vline!(plt2, [Kc], color=:red, linestyle=:dash, label="Kc")
        display(plt2)
        #=
        # Fig. 14: ln|s| vs ln t
        plt2 = plot(logt, logs, seriestype=:scatter, ms=6,
                    xlabel=L"\ln{t}", ylabel=L"(\ln{Λ})'(K_c)",
                    title="Scaling of slopes", label="data")
        plot!(plt2, logt, logs_fit, lw=2, label="fit ν≈$(round(ν,digits=3))")
        display(plt2)=#
    end

    return (ν=ν, slope=slope, intercept=intercept,
            s_vals=s_vals, s_errs=s_errs,
            lnΛc_vals=lnΛc_vals)
end

Kc = 1.165
dim=3

res = finite_time_linear_scaling(K_vals, t_vals, p2_mat, dim;
                                 Kc = Kc, ΔKfit = 0.25, n_kicks_i=5)

println("\n===== Linear finite-time-scaling results =====")
println("ν = $(round(res.ν,digits=3))  (from slope = $(round(res.slope,digits=4)))")

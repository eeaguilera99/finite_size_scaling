using Plots
using LaTeXStrings

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
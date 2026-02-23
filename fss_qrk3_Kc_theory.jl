include("imp_data_theory.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling


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
                dof = length(ufit[2:end]) - length(pbest)
        else
                χ2 = sum((residuals ./ uerrfit).^2)
                dof = length(ufit) - length(pbest)
        end
        χ2_red = χ2 / dof

        # Unpack results
        β0, A, ν, Kc = pbest
        err_β0, err_A, err_ν, err_Kc = perr


        return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
                err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


dim1 = 3
dim2 = 3   # spatial dimension


# Perform collapse
t_transient = 28
avg = true
amp = 2
K_guess_index = 0# index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit
type_data = 2 # 1 for kinetic energy, 2 for nc^2, 3 for both

function perform_Kc_analysis(K_vals, t_vals, p2_mat, p2_err_mat, nc_mat, nc_err_mat; q=2, d1=dim1, d2=dim2, transient=t_transient, avg=avg, amp=amp,
         K_guess_index=K_guess_index, fit_k_filter=fit_k_filter) 
        if q == 1
                _, shifts1, _, _, _, _ = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=avg, loc_amp=amp), 
                p2_err_mat; d=d1, n_kicks_i=transient, n_kicks_f=0)
                shifts1err = shifts_parametric_mc(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=avg, loc_amp=amp), 
                p2_err_mat; d=d1, nbins=30, nmc=1000)[2]
                xi1 = exp.(shifts1)
                xi1err = xi1.*shifts1err
                results1 = fit_xi_offset_LsqFit(K_vals, xi1, xi1err; n_k_filter=fit_k_filter)
                # Plot fit
                Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
                plt_fit1 = plot(K_vals, xi1, seriestype=:scatter, ms=6,
                                xlabel=L"κ", ylabel=L"ξ(κ)",
                                title=latexstring("\$E_k\$, \$κ_c≈ $(round(results1.Kc,digits=3))\$, \$d=$(d1)\$"),
                                label="data")
                ξfit1 = (1)./(results1.β0 .+ results1.A .* abs.(Kgrid .- results1.Kc).^(abs(results1.ν)))
                plot!(plt_fit1, Kgrid, ξfit1, lw=2,
                        label="fit (ν ≈ $(round(abs(results1.ν),digits=3)))")
                vline!(plt_fit1, [results1.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
                display(plt_fit1)
                println("\n===== Critical fit results with offset and error bars =====")
                println("κc1  ≈ $(results1.Kc)  ± $(results1.err_Kc)")
                println("ν1   ≈ $(results1.ν)   ± $(results1.err_ν)")
                println("A1   ≈ $(results1.A)   ± $(results1.err_A)")
                println("β_01  ≈ $(results1.β0)  ± $(results1.err_β0)")
                println("χ²1  = $(results1.χ2),  χ²_red1 = $(results1.χ2_red)")
        end
        if q == 2
                _, shifts2, _, _, _, _ = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=avg, loc_amp=amp), 
                nc_err_mat; d=d2, n_kicks_i=transient, n_kicks_f=0)
                shifts2err = shifts_parametric_mc(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=avg, loc_amp=amp), 
                nc_err_mat; d=d2, nbins=30, nmc=1000)[2]
                xi2 = exp.(shifts2)
                xi2err = xi2.*shifts2err
                results2 = fit_xi_offset_LsqFit(K_vals, xi2, xi2err; n_k_filter=fit_k_filter, K_val_g=K_guess_index)
                # Plot fit
                Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
                plt_fit2 = plot(K_vals, xi2, yerror=xi2err, seriestype=:scatter, ms=6,
                                xlabel=L"κ", ylabel=L"ξ(κ)",
                                title=latexstring("\$a_s=$(a_s)a_0\$, \$κ_c≈ $(round(results2.Kc,digits=3))\$, \$d=$(d2)\$"),#\$1/n_c^2\$ \$a_s=$(a_s)a_0\$
                                label="data")
                ξfit2 = (1)./(results2.β0 .+ results2.A .* abs.(Kgrid .- results2.Kc).^(abs(results2.ν)))
                plot!(plt_fit2, Kgrid, ξfit2, lw=2,
                        label="fit (ν ≈ $(round(abs(results2.ν),digits=3)))")
                vline!(plt_fit2, [results2.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
                display(plt_fit2)
                println("\n===== Critical fit results with offset and error bars =====")
                println("κc2  ≈ $(results2.Kc)  ± $(results2.err_Kc)")
                println("ν2   ≈ $(results2.ν)   ± $(results2.err_ν)")
                println("A2   ≈ $(results2.A)   ± $(results2.err_A)")
                println("β_02  ≈ $(results2.β0)  ± $(results2.err_β0)")
                println("χ²2  = $(results2.χ2),  χ²_red2 = $(results2.χ2_red)") 
        end
        if q == 3
                _, shifts1, _, _, _, _ = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=avg, loc_amp=amp), 
                p2_err_mat; d=d1, n_kicks_i=transient, n_kicks_f=0)
                shifts1err = shifts_parametric_mc(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=avg, loc_amp=amp), 
                p2_err_mat; d=d1, nbins=30, nmc=1000)[2]
                xi1 = exp.(shifts1)
                xi1err = xi1.*shifts1err
                results1 = fit_xi_offset_LsqFit(K_vals, xi1, xi1err; n_k_filter=fit_k_filter)
                _, shifts2, _, _, _, _ = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=avg, loc_amp=amp),
                nc_err_mat; d=d2, n_kicks_i=transient, n_kicks_f=0)
                shifts2err = shifts_parametric_mc(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=avg, loc_amp=amp),
                nc_err_mat; d=d2, nbins=30, nmc=1000)[2]
                xi2 = exp.(shifts2)
                xi2err = xi2.*shifts2err
                results2 = fit_xi_offset_LsqFit(K_vals, xi2, xi2err; n_k_filter=fit_k_filter, K_val_g=K_guess_index)
                # Plot fit
                Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
                plt_fit1 = plot(K_vals, xi1, seriestype=:scatter, ms=6,
                                xlabel=L"κ", ylabel=L"ξ(κ)",
                                title=latexstring("\$E_k\$, \$κ_c≈ $(round(results1.Kc,digits=3))\$, \$d=$(d1)\$"),
                                label="data")
                ξfit1 = (1)./(results1.β0 .+ results1.A .* abs.(Kgrid .- results1.Kc).^(abs(results1.ν)))
                plot!(plt_fit1, Kgrid, ξfit1, lw=2,
                        label="fit (ν ≈ $(round(abs(results1.ν),digits=3)))")
                vline!(plt_fit1, [results1.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
                plt_fit2 = plot(K_vals, xi2, yerror=xi2err, seriestype=:scatter, ms=6,
                                xlabel=L"κ", ylabel=L"ξ(κ)",
                                title=latexstring("\$a_s=$(a_s)a_0\$, \$κ_c≈ $(round(results2.Kc,digits=3))\$, \$d=$(d2)\$"),#\$1/n_c^2\$ \$a_s=$(a_s)a_0\$
                                label="data")
                ξfit2 = (1)./(results2.β0 .+ results2.A .* abs.(Kgrid .- results2.Kc).^(abs(results2.ν)))
                plot!(plt_fit2, Kgrid, ξfit2, lw=2,
                        label="fit (ν ≈ $(round(abs(results2.ν),digits=3)))")
                vline!(plt_fit2, [results2.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
                display(plot(plt_fit1, plt_fit2, layout=(1,2), size=(1050,550), 
                suptitle=latexstring("\$ ξ(κ)\$ with offset, \$a_s=$(a_s)a_0\$"), 
                bottom_margin=5Plots.mm, left_margin=5Plots.mm))
                #print results
                println("\n===== Critical fit1 results with offset and error bars =====")
                println("κc1  ≈ $(results1.Kc)  ± $(results1.err_Kc)")
                println("ν1   ≈ $(results1.ν)   ± $(results1.err_ν)")
                println("A1   ≈ $(results1.A)   ± $(results1.err_A)")
                println("β_01  ≈ $(results1.β0)  ± $(results1.err_β0)")
                println("χ²1  = $(results1.χ2),  χ²_red1 = $(results1.χ2_red)")

                println("\n===== Critical fit2 results with offset and error bars =====")
                println("κc2  ≈ $(results2.Kc)  ± $(results2.err_Kc)")
                println("ν2   ≈ $(results2.ν)   ± $(results2.err_ν)")
                println("A2   ≈ $(results2.A)   ± $(results2.err_A)")
                println("β_02  ≈ $(results2.β0)  ± $(results2.err_β0)")
                println("χ²2  = $(results2.χ2),  χ²_red2 = $(results2.χ2_red)") 
        end
end

perform_Kc_analysis(K_vals, t_vals, p2_mat, p2_err_mat, nc_mat, nc_err_mat; q=type_data, d1=dim1, d2=dim2, transient=t_transient, avg=avg, amp=amp,
         K_guess_index=K_guess_index, fit_k_filter=fit_k_filter)

        


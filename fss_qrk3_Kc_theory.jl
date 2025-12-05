include("fss_qkr3_timescaling_theory.jl")  # for finite_time_scaling


"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi; exclude_tol_frac=0.05)
    xi=(1)./xi
    # Model function
    model(K, p) = abs(p[1]) .+ p[2] .* abs.(K .- p[4]).^(abs(p[3]))   #1/ξ(K) = β₀ + A|K−Kc|^{−ν}
    names = ["β₀", "A", "ν", "Kc"]

    # Initial guess
    β₀₀ = minimum(xi)*0.5
    A₀  = maximum(xi)
    ν₀  = 0.5
    Kc₀ = K_vals[argmin(xi)]  # where ξ is largest
    p0 = [β₀₀, A₀, ν₀, Kc₀]

    # Mask out values too close to trial Kc₀
    ΔK = maximum(K_vals) - minimum(K_vals)
    mask = abs.(K_vals .- Kc₀) .> exclude_tol_frac*ΔK
    Kfit, ξfit = K_vals[mask], xi[mask]

    # Perform nonlinear least squares fit
    fit = curve_fit(model, Kfit, ξfit, p0)
    pbest = coef(fit)
    covar = estimate_covar(fit)
    perr  = sqrt.(diag(covar))

    # goodness of fit
    residuals = ξfit .- model(Kfit, pbest)
    χ2 = sum((residuals ./ model(Kfit, pbest)).^2)
    dof = length(ξfit) - length(pbest)
    χ2_red = χ2 / dof

    # Unpack results
    β0, A, ν, Kc = pbest
    err_β0, err_A, err_ν, err_Kc = perr


    return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
            err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


dim1 = 3
dim2 = 3   # spatial dimension

#filter Kick strength
_, p2_mat, p2_err_mat = filter_K(K_vals, p2_mat, p2_err_mat; n_kkicks_i=2, n_kkicks_f=0)
K_vals, nc_mat, nc_err_mat = filter_K(K_vals, nc_mat, nc_err_mat; n_kkicks_i=2, n_kkicks_f=0)

# Perform collapse
res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=true, loc_amp=2), 
p2_err_mat; d=dim1, n_kicks_i=40, n_kicks_f=0)
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=true, loc_amp=2), 
nc_err_mat; d=dim2, n_kicks_i=40, n_kicks_f=0)

# ===== Example usage =====
xi1 = exp.(shifts1)
xi2 = exp.(shifts2)
results1 = fit_xi_offset_LsqFit(K_vals, xi1)
results2 = fit_xi_offset_LsqFit(K_vals, xi2)
#=
# Plot raw data
plt_raw1 = plot(K_vals, xi1, seriestype=:scatter, ms=6,
                xlabel="K", ylabel="ξ(K)", title=L" E_k", label="data")
plt_raw2 = plot(K_vals, xi2, seriestype=:scatter, ms=6,
                xlabel="K", ylabel="ξ(K)", title=L"1/nc^2", label="data")
display(plot(plt_raw1, plt_raw2, layout=(1,2), size=(1000,400), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))=#
        
# Plot fit
Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
plt_fit1 = plot(K_vals, xi1, seriestype=:scatter, ms=6,
                xlabel=L"κ", ylabel=L"ξ(κ)",
                title=latexstring("\$E_k\$, \$κ_c≈ $(round(results1.Kc,digits=3))\$, \$d=$(dim1)\$"),
                label="data")
ξfit1 = (1)./(results1.β0 .+ results1.A .* abs.(Kgrid .- results1.Kc).^(abs(results1.ν)))
plot!(plt_fit1, Kgrid, ξfit1, lw=2,
        label="fit (ν ≈ $(round(abs(results1.ν),digits=3)))")
vline!(plt_fit1, [results1.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")

plt_fit2 = plot(K_vals, xi2, seriestype=:scatter, ms=6,
                xlabel=L"κ", ylabel=L"ξ(κ)",
                title=latexstring("\$1/n_c^2\$, \$κ_c≈ $(round(results2.Kc,digits=3))\$, \$d=$(dim2)\$"),
                label="data")
ξfit2 = (1)./(results2.β0 .+ results2.A .* abs.(Kgrid .- results2.Kc).^(abs(results2.ν)))
plot!(plt_fit2, Kgrid, ξfit2, lw=2,
        label="fit (ν ≈ $(round(abs(results2.ν),digits=3)))")
vline!(plt_fit2, [results2.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")

display(plot(plt_fit1, plt_fit2, layout=(1,2), size=(1050,550), 
suptitle=latexstring("\$ ξ(κ)\$ with offset, \$a_s=$(a_s)a_0\$"), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))


#println(results)

println("\n===== Critical fit results with offset and error bars =====")
println("κc1  ≈ $(results1.Kc)  ± $(results1.err_Kc)")
println("ν1   ≈ $(results1.ν)   ± $(results1.err_ν)")
println("A1   ≈ $(results1.A)   ± $(results1.err_A)")
println("β_01  ≈ $(results1.β0)  ± $(results1.err_β0)")
println("χ²1  = $(results1.χ2),  χ²_red1 = $(results1.χ2_red)")

println("\n===== Critical fit results with offset and error bars =====")
println("κc2  ≈ $(results2.Kc)  ± $(results2.err_Kc)")
println("ν2   ≈ $(results2.ν)   ± $(results2.err_ν)")
println("A2   ≈ $(results2.A)   ± $(results2.err_A)")
println("β_02  ≈ $(results2.β0)  ± $(results2.err_β0)")
println("χ²2  = $(results2.χ2),  χ²_red2 = $(results2.χ2_red)")     

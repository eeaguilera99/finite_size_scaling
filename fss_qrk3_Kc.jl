include("fss_qkr3_timescaling.jl")  # for finite_time_scaling


"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi, xierr; exclude_tol_frac=0.02)
    u = (1)./xi
    # Model function
    model(K, p) = abs(p[1]) .+ p[2] .* abs.(K .- p[4]).^(abs(p[3]))   #1/ξ(K) = β₀ + A|K−Kc|^{ν}
    names = ["β₀", "A", "ν", "Kc"]

    # Initial guess
    β₀₀ = maximum(u)*0.5
    A₀  = minimum(u)
    ν₀  = 1
    Kc₀ = K_vals[argmin(u)]  # where ξ is largest
    p0 = [β₀₀, A₀, ν₀, Kc₀]

    # Mask out values too close to trial Kc₀
    ΔK = maximum(K_vals) - minimum(K_vals)
    mask = abs.(K_vals .- Kc₀) .> exclude_tol_frac*ΔK
    Kfit, ufit = K_vals[mask], u[mask]

    # Perform nonlinear least squares fit
    fit = curve_fit(model, Kfit, ufit, p0)
    pbest = coef(fit)
    covar = estimate_covar(fit)
    perr  = sqrt.(diag(covar))

    # goodness of fit
    residuals = ufit .- model(Kfit, pbest)
    uerr = u.^2 .*xierr
    χ2 = sum((residuals ./ uerr[mask]).^2)
    dof = length(ufit) - length(pbest)
    χ2_red = χ2 / dof

    # Unpack results
    β0, A, ν, Kc = pbest
    err_β0, err_A, err_ν, err_Kc = perr


    return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
            err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


dim = 3 # spatial dimension


# Perform collapse
res, shifts, shiftserr, X, Y, Yerr, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=2, n_kicks_f=0)
# ===== Example usage =====
xi_fit = exp.(shifts)[2:end]
xierr_fit = (xi.*shiftserr)[2:end]  #error propagation
K_vals_fit = K_vals[2:end]
results = fit_xi_offset_LsqFit(K_vals_fit, xi_fit, xierr_fit)#exclude frist point from fit

#=
# Plot raw data
plt_raw = plot(K_vals, (1)./xi, seriestype=:scatter, ms=6,
 xlabel="K", ylabel="ξ(K)", title="Raw ξ(K) data", label="data")
display(plt_raw)=#

# Plot fit
Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
plt_fit = plot(K_vals_fit, xi_fit, seriestype=:scatter, ms=6,
    xlabel=L"κ", ylabel=L"ξ(κ)",
    title=latexstring("\$ξ(κ)\$ with offset \$a_s=$(a_s)a_0\$, \$κ_c\$≈$(round(results.Kc,digits=5))"),
    label="data")
Kgrid = range(minimum(K_vals_fit), maximum(K_vals_fit), length=400)
ξfit = (1)./(results.β0 .+ results.A .* abs.(Kgrid .- results.Kc).^(abs(results.ν)))
plot!(plt_fit, Kgrid, ξfit, lw=2,
    label="fit (ν ≈ $(round(abs(results.ν),digits=3)))")
vline!(plt_fit, [results.Kc], linestyle=:dash, color=:red, label=L"\kappa_c")
display(plt_fit)
       
#println(results)

println("\n===== Critical fit results with offset and error bars =====")
println("κc  ≈ $(results.Kc)  ± $(results.err_Kc)")
println("ν   ≈ $(results.ν)   ± $(results.err_ν)")
println("A   ≈ $(results.A)   ± $(results.err_A)")
println("β_0  ≈ $(results.β0)  ± $(results.err_β0)")
println("χ²  = $(results.χ2),  χ²_red = $(results.χ2_red)")


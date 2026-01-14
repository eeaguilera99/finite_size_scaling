include("fss_qkr3_timescaling.jl")  # for finite_time_scaling


"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=2, exclude_tol_frac=0.02)

    #filter points for fit    
    if n_k_filter != 1    
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
    Kc₀ = K_vals[argmin(u)+1]  # where ξ is largest
    println(Kc₀)
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


dim = 3 # spatial dimension
#K_vals, p2_mat, p2_err_mat = filter_K(K_vals, p2_mat, p2_err_mat; n_kkicks_i=2, n_kkicks_f=0)

# Perform collapse
_, shifts, _, _, _, _, _ = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5, n_kicks_f=0)
shiftserr = shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, nbins=30, nmc=1000)[2]
# ===== Example usage =====

xi = exp.(shifts)
xierr = (xi.*shiftserr)  #error propagation
results = fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=1)#exclude frist point from fit

#=
# Plot raw data
plt_raw = plot(K_vals, (1)./xi, seriestype=:scatter, ms=6,
 xlabel="K", ylabel="ξ(K)", title="Raw ξ(K) data", label="data")
display(plt_raw)=#

# Plot fit
Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
plt_fit = plot(K_vals, xi, yerror=xierr, seriestype=:scatter, ms=6,
    xlabel=L"κ", ylabel=L"ξ(κ)",
    title=latexstring("\$ξ(κ)\$ \$a_s=$(a_s)a_0\$, \$κ_c\$≈$(round(results.Kc,digits=5))"),
    label="data")

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
#=
#save data
d1 = DataFrame(xi', :auto)
CSV.write("ξ(k)_data_220.csv", d1)=#
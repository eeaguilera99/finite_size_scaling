include("fss_qkr3_timescaling.jl")  # for finite_time_scaling

"""
Plot kinetic energy as function of time at criticality for exponent of difussion.
"""



#filter Nkicks range
t_vals, p2_mat, p2_err_mat = filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=2, n_kicks_f=0)

dim = 3
res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)


xi1 = exp.(shifts1)
K_c1_i = 10#argmax(xi1)+1

#Plot straight lines 
model(t, p) = p[1] .* t .+ p[2]  # y = m*x + b
guess1 = [2/3, 0.0]  # Initial guess for [m, b]
fit1 = curve_fit(model, log.(t_vals), log.(p2_mat[K_c1_i, :]), guess1)
fit_params1 = coef(fit1)
fit_errs1 = standard_errors(fit1)
α_err = fit_errs1[1]

#χ^2 goodness of fit
residuals1 = log.(p2_mat[K_c1_i, :]) .- model(log.(t_vals), fit_params1)
χ2_1 = sum((residuals1 ./ (p2_err_mat[K_c1_i, :]./p2_mat[K_c1_i, :])).^2)
dof1 = length(t_vals) - length(fit_params1)
println("χ2 =", χ2_1)
println("χ2 rel =", χ2_1/dof1)

#Plot
plt1 = plot(title="Kinetic energy evoltion at near Kc",
    xlabel="Time (kicks)", ylabel=L"1/n_c^2")
plot!(plt1, t_vals, p2_mat[K_c1_i, :], xscale=:log10, yscale=:log10, marker=:o, label="K=$(round(K_vals[K_c1_i], digits=3))", ms=3)
plot!(plt1, t_vals, exp.(model(log.(t_vals),fit_params1)), xscale=:log10, yscale=:log10, label=latexstring("\$α≈$(round(fit_params1[1], digits=3)) ± $(round(α_err, digits=3))\$"))
display(plt1)
println("Approximate dimension from fit: d ≈ $(round(2/fit_params1[1], digits=3)) ± $(round(2*α_err/fit_params1[1]^2, digits=3))")

#=
#Compare with localized and diffusive cases
K_diff_i = 14
K_loc_i = 3
fit2 = curve_fit(model, log.(t_vals), log.(p2_mat[K_loc_i, :]), guess1)
fit_params2 = coef(fit2)
fit3 = curve_fit(model, log.(t_vals), log.(p2_mat[K_diff_i, :]), guess1)
fit_params3 = coef(fit3)
plt2 = plot(title="Kinetic energy evoltion at different regimes",
    xlabel="t )", ylabel=L"1/n_c^2")
plot!(plt2, log.(t_vals), log.(p2_mat[K_c1_i, :]), marker=:o, label="Criticality K=$(round(K_vals[K_c1_i], digits=3))", ms=3)
plot!(plt2, log.(t_vals), model(log.(t_vals),fit_params1), lc=:black, label=latexstring("\$α≈$(round(fit_params1[1], digits=3))\$"))
plot!(plt2, log.(t_vals), log.(p2_mat[K_loc_i, :]), marker=:o, label="Localized K=$(round(K_vals[K_loc_i], digits=3))", ms=3)
plot!(plt2, log.(t_vals), model(log.(t_vals),fit_params2), lc=:black, label=latexstring("\$α≈$(round(fit_params2[1], digits=3))\$"))
plot!(plt2, log.(t_vals), log.(p2_mat[K_diff_i, :]), marker=:o, label="Diffusive K=$(round(K_vals[K_diff_i], digits=3))", ms=3)
plot!(plt2, log.(t_vals), model(log.(t_vals),fit_params3), lc=:black, label=latexstring("\$α≈$(round(fit_params3[1], digits=3))\$"))
display(plt2)=#

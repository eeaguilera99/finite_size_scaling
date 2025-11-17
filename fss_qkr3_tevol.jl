include("fss_qkr3_timescaling.jl")

dim = 3

res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=true, loc_amp=6), p2_err_mat; d=dim, n_kicks_i=7, n_kicks_f=0)
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=true, loc_amp=6), nc_err_mat; d=dim, n_kicks_i=7, n_kicks_f=0)

xi1 = exp.(shifts1)
xi2 = exp.(shifts2)
K_c1_i = argmax(xi1)+1
K_c2_i = argmax(xi2)+1

#Plot straight lines 
model(t, p) = p[1] .* t .^ (p[3]) .+ p[2]  # y = m*x^a + b
guess1 = [1.0, 0.0, 2/3]  # Initial guess for [m, b, a]
fit1 = curve_fit(model, log.(t_vals), log.(p2_mat[K_c1_i, :]), guess1)
fit_params1 = coef(fit1)
fit2 = curve_fit(model, log.(t_vals), log.(nc_mat[K_c2_i, :]), guess1)
fit_params2 = coef(fit2)

#plot data and fits
#plotlyjs()
plt1 = plot(title=L"E_k",
    xlabel="Time (kicks)", ylabel=L"E_k")
plot!(plt1, t_vals, log.(p2_mat[K_c1_i, :]), xscale=:log10, marker=:o, label="K=$(round(K_vals[K_c1_i], digits=3))", ms=3)
plot!(plt1, t_vals, model(log.(t_vals),fit_params1), xscale=:log10, label=latexstring("\$d≈$(round(2/fit_params1[3], digits=3))\$"))

plt2 = plot(title=L"1/n_c^2",
    xlabel="Time (kicks)", ylabel=L"1/n_c^2")
plot!(plt2, t_vals, log.(nc_mat[K_c2_i, :]) , xscale=:log10, marker=:o, label="K=$(round(K_vals[K_c2_i], digits=3))", ms=3)
plot!(plt2, t_vals, model(log.(t_vals),fit_params2), xscale=:log10, label=latexstring("\$d≈$(round(2/fit_params2[3], digits=3))\$"))

display(plot(plt1, plt2, suptitle=latexstring("Critical data at \$κ_c\$ with, \$a_s=$(a_s)a_0\$"), layout=(1,2), size=(1000,400), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))
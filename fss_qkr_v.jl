include("fss_qkr3_timescaling.jl")  # for fit_xi_offset_vsK

dim=3  # spatial dimension
# Perform collapse
res, shifts, X, Y, Yerr, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)
xi=exp.(shifts)

#filter Nkicks range
n_kicks_i=5
n_kicks_f=0
n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
t_vals = t_vals[n_kicks_i:n_Nkicks_f]
p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]

app_Kc_index = argmax(xi)
E_vals = p2_mat[app_Kc_index, :]

model(t,p) = p[1] .+ p[2].*t.^(2/(3*p[3]))  # model: ln(Λ) = ln(ξ0) + A * t^(1/(3*ν))

fit = curve_fit(model, t_vals, E_vals, [0, 1.0, 1.5])
ν_fit = fit.param[3]

plot(t_vals, E_vals, seriestype=:scatter, xscale=:log10, yscale=:log10, label="ν=$(ν_fit)", xlabel="t", ylabel="⟨p²⟩ at Kc",
        title="Fit of Λ vs t at Kc")
plot!(t_vals, model(t_vals, fit.param), lw=2, label="fit")

#display(plt_fit)
#fit_nu_vs_t(t_vals, Y, xi)
include("fss_qkr3_timescaling.jl")  # for fit_xi_offset_vsK

#calculate ν for Ln(Λ) vs t at Kc
function fit_nu_vs_t(t_vals, Y, xi; ngrid=400, n_kicks_i=5, n_kicks_f=0)
    #filter Nkicks range
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]

    app_Kc_index = argmax(xi)
    Ln_Λ_vals = Y[app_Kc_index, :]

    model(t,p) = t.^(1/(3*p[1]))  # model: ln(Λ) = ln(ξ0) + A * t^(1/(3*ν))

    fit = curve_fit(model, t_vals, Ln_Λ_vals, [1.0])
    ν_fit = fit.param[1]

    plt_fit = plot(t_vals, Ln_Λ_vals, seriestype=:scatter, xscale=:ln, yscale=:ln, label="ν=$(ν_fit)", xlabel="t", ylabel="ln(Λ(Kc,t))",
        title="Fit of ln(Λ) vs t at Kc")
    #plot!(plt_fit, t_vals, model(t_vals, fit.param), lw=2, label="fit")
    display(plt_fit)

    return ν_fit
end

dim=3  # spatial dimension
# Perform collapse
res, shifts, X, Y, Yerr, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)
xi=exp.(shifts)
fit_nu_vs_t(t_vals, Y, xi)
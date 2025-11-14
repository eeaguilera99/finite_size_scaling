include("fss_qkr3_timescaling.jl")

a_s = 1038  # scattering length label for plots

res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=0, loc_amp=6), p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=0, loc_amp=6), nc_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)

xi1 = exp.(shifts1)
xi2 = exp.(shifts2)
K_c1_i = argmax(xi1)+1
K_c2_i = argmax(xi2)+1


plt1 = plot(title=L"E_k",
    xlabel="Time (kicks)", ylabel=L"E_k")
plot!(plt1, t_vals, p2_mat[K_c1_i, :], xscale=:log10, yscale=:log10, marker=:o, label="K=$(round(K_vals[K_c1_i], digits=3))", ms=3)

plt2 = plot(title=L"1/n_c^2",
    xlabel="Time (kicks)", ylabel=L"1/n_c^2")
plot!(plt2, t_vals, nc_mat[K_c2_i, :], xscale=:log10, yscale=:log10, marker=:o, label="K=$(round(K_vals[K_c2_i], digits=3))", ms=3)
display(plot(plt1, plt2, suptitle=latexstring("Critical data at \$κ_c\$ with \$d=$(dim)\$, \$a_s=$(a_s)a_0\$"), layout=(1,2), size=(1000,400), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))
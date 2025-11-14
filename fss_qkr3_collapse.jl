include("fss_qkr3_timescaling.jl")

"""
    finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=30)

Perform finite-time scaling collapse of Anderson transition data.

# Arguments
- `K_vals::Vector`: Kick strengths (length M).
- `t_vals::Vector`: Times (length N).
- `p2_mat::Matrix`: M×N matrix of ⟨p²⟩ values.
- `p2_err_mat::Matrix`: M×N matrix of errors (same shape).
- `nbins::Int`: Number of bins in Y direction (default = 30).

# Returns
- `shifts::Vector`: Optimal horizontal shifts aᵢ = ln ξ(Kᵢ).
- `(X, Y)`: Arrays of logarithmic coordinates.
"""

dim = 3  # spatial dimension
a_s = 920  # scattering length label for plots

# Perform collapse
res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=0, loc_amp=6), p2_err_mat; d=dim, n_kicks_i=13, n_kicks_f=0)
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=0, loc_amp=6), nc_err_mat; d=dim, n_kicks_i=13, n_kicks_f=0)

#=
# plots for ⟨p²⟩
plt1 = plot(title=L"E_k",
    xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X1[:], Y1[i,:], marker=:o, label="")#, yerror=Yerr[i,:]κ="*string(round(K, digits=3))
end

# plots for nc^2
plt2 = plot(title=L"1/n_c^2",
    xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt2, X2[:], Y2[i,:], marker=:o, label="")#, yerror=Yerr[i,:]κ="*string(round(K, digits=3))
end

display(plot(plt1, plt2, suptitle=latexstring("Raw data  \$d=$(dim)\$, \$a_s=$(a_s)a_0\$"), layout=(1,2), size=(1000,400), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))
=#
#collapse plot
plt3 = plot(title=L"E_k",
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt3, X1[:] .+ shifts1[i], Y1[i,:], marker=:o, label="")
end

plt4 = plot(title=L"$1/n_c^2$",
    xlabel="ln(ξ/N^(1/d))")
for (i,K) in enumerate(K_vals)
    plot!(plt4, X2[:] .+ shifts2[i], Y2[i,:], marker=:o, label="")
end

display(plot(plt3, plt4, suptitle=latexstring("Collapse \$d=$(dim)\$, \$a_s=$(a_s)a_0\$"), layout=(1,2), size=(1000,400), 
bottom_margin=5Plots.mm, left_margin=5Plots.mm))
#savefig(plt4, "fss_collapsed_data_nc2_d$(dim).png")=#

#println("Fit quality 1: ", s_rel1)
#println("Fit quality 2: ", s_rel2)

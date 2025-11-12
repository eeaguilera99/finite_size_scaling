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
a_s = L"1038a_0"  # scattering length label for plots

# Perform collapse
res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(p2_mat, p=0), p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals, t_vals, apply_mov_av_matrix(nc_mat, p=0), nc_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)

# plots for ⟨p²⟩
plt1 = plot(title="Raw data "*L"E_k, a_s="*a_s,
    xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X1[:], Y1[i,:], marker=:o, label="")#, yerror=Yerr[i,:]κ="*string(round(K, digits=3))
end
display(plt1)
#savefig(plt1, "fss_raw_data_d$(dim).png")

#collapse plot
plt2 = plot(title="Data collapse d=$(dim), "*L"E_k, a_s="*a_s,
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt2, X1[:] .+ shifts1[i], Y1[i,:], marker=:o, label="")
end
display(plt2)
#savefig(plt2, "fss_collapsed_data_d$(dim).png")


# plots for nc^2
plt3 = plot(title="Raw data "*L"1/n_c^2, a_s="*a_s,
    xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt3, X2[:], Y2[i,:], marker=:o, label="")#, yerror=Yerr[i,:]κ="*string(round(K, digits=3))
end
display(plt3)
#savefig(plt3, "fss_raw_data_nc2_d$(dim).png")

#collapse plot
plt4 = plot(title="Data collapse d=$(dim), "*L"1/n_c^2, a_s="*a_s,
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt4, X2[:] .+ shifts2[i], Y2[i,:], marker=:o, label="")
end
display(plt4)
#savefig(plt4, "fss_collapsed_data_nc2_d$(dim).png")

println("Fit quality 1: ", s_rel1)
println("Fit quality 2: ", s_rel2)

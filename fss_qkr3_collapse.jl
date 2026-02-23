include("imp_data_ex.jl")
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

#function FTT_collapse()

dim = 3 # spatial dimension

# Perform collapse
res, shifts, X, Y, Yerr, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5, n_kicks_f=0)

# Plot before collapse
plt1 = plot(title=latexstring("Raw data \$d=$(dim)\$, \$a_s=$(a_s)a_0\$"),
    xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label="K="*string(round(K, digits=3)))
end
display(plt1)
#savefig(plt1, "fss_raw_data_d$(dim).png")

# Plot after collapse
plt2 = plot(title=latexstring("Data collapse \$d=$(dim)\$, \$a_s=$(a_s)a_0\$"),
    xlabel=latexstring("\$\\ln(\\xi/N^{1/d})\$"), ylabel=latexstring("\$\\ln(\\Lambda)\$"))
for (i,K) in enumerate(K_vals)
    plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label="")
end
display(plt2)
#savefig(plt2, "fss_collapsed_data_d$(dim).png")

println("Fit quality: ", s_rel)
#=
#write csv file with scaling data
df1 = DataFrame(Y, :auto)
CSV.write("logΛ_220_d=3.csv", df1)

df2 = DataFrame(X, :auto)
CSV.write("logN^1d_220_d=3.csv", df2)

df3 = DataFrame(shifts', :auto)
CSV.write("logξ_220_d=3.csv", df3)

df4 = DataFrame(Yerr', :auto)
CSV.write("logΛ_err_220_d=3.csv", df4)=#
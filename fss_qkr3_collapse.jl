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

# Perform collapse
res, shifts, X, Y, Yerr, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)
#=
# Plot before collapse
plt1 = plot(title="Raw data (before shifts)",
    xlabel="ln(t^(-1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label="")
end
display(plt1)=#

# Plot after collapse
plt2 = plot(title="Data collapse (after optimal shifts)",
    xlabel="ln(ξ/t^(1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label="")
end
display(plt2)

#quality factors
var_0 = tot_variance(shifts.*0, X, Y)[1]
var_f = tot_variance(shifts, X, Y)[1]
q_1 = var_f / var_0
q_2 = 1/s_rel

println("Fit quality: ", s_rel)
println("Quality factors: q_1 = ", q_1, ", q_2 = ", q_2)

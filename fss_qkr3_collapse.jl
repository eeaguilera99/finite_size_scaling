include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")
include("fss_qrk3_timescaling_analysis.jl")

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

# Perform collapse
#perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data; data_type="Ex nc^-2", raw=true, d=dim)
q_diff, q_err_diff, q_loc, q_err_loc = perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, dim; plotshow=true)
println("Quality of collapse: q_diff = $q_diff, q_loc = $q_loc")
println("Errors in quality (thoery): q_err_diff = $q_err_diff, q_err_loc = $q_err_loc")

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
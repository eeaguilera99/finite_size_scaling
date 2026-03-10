include("fss_qkr3_timescaling.jl")
using CSV, DataFrames, Plots, LaTeXStrings

a_s = 220
K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_matrix.csv", DataFrame; header=false))   
p2_err_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_err_matrix.csv", DataFrame; header=false))     # Errors

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

dim = 3 # spatial dimension
transients = 4 .*ones(Int, 15)#[2, 2, 3, 3, 4, 1, 4, 2, 2, 3, 4, 2, 2, 3, 3]
data_type = "Ex"
shifts_data, X_data, Y_data, Yerr_data, s_rel_data = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat, transients; d=dim, nbins=30)
#_, shiftserr_data, _ = shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, nbins=30, nmc=1000)


include("fss_qrk3_timescaling_analysis.jl")

perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data; data_type=data_type, plot_label_b=false, ploterr=false, raw=false, d=dim)
#perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data, dim, a_s, "Ex"; Kc_offset=2, plotshow=false)
include("fss_qkr3_timescaling.jl")
using CSV, DataFrames, Plots, LaTeXStrings

#a_s = "κ=0.87"
a_s = 220
K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
#K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/a_s.csv", DataFrame; header=false)))             # interactions
t_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/Number_of_kicks.csv", DataFrame; header=false)))  # Times
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

D = 3 # spatial dimension
transient = 5
data_type = "Ex"
#shifts_data, X_data, Y_data, Yerr_data, s_rel_data = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=50, d=D, n_kicks_i=transient, n_kicks_f=0)
#_, shiftserr_data, _ = shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=D, nbins=30, nmc=1000)


#Collapse snc method 
V_guess = [1, 1, 1.1, 0.5, -20, 0.1]
X_data, Y_data, Yerr_data, pbest, perr, shifts_data, χ2, χ2_red = finite_time_scaling2(K_vals, t_vals, p2_mat, p2_err_mat, D, V_guess; transient=transient)


include("fss_qrk3_timescaling_analysis.jl")

#perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data; data_type=data_type, raw=false, d=D, save=false)
#perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data, D, a_s, data_type; Kc_offset=2, plotshow=false)

perform_collapse2(K_vals, X_data, Y_data, Yerr_data, shifts_data, χ2, χ2_red, pbest, perr, D, a_s; showxi=true)
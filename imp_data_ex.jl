include("fss_qkr3_timescaling.jl")
include("fss_qrk3_timescaling_analysis.jl")
using CSV, DataFrames

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

D = 3  # spatial dimension
transient = 5
data_type = "Ex"

#Scaling data
X_data, Y_data, Yerr_data = finite_time_scaling_data(t_vals, p2_mat, p2_err_mat; d=D, n_kicks_i=transient, n_kicks_f=0)

#Scaling variance optimizaztion
#shifts_data, s_rel_data = finite_time_scaling(X_data, Y_data; nbins=50)
#_, shiftserr_data, _ = shifts_parametric_mc(t_vals, p2_mat, p2_err_mat; d=D, nbins=30, nmc=1000)

#perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data; data_type=data_type, raw=false, d=D, save=false)
#perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data, D, a_s, data_type; Kc_offset=2, plotshow=false)

#Scaliong collapse Taylor fitting
#V_guess = [1, 0.1, 0.8, 1, -20, 10] #(b1, b2, Kc, ν, F00, ξsat)
#pbest, perr, shifts_data, χ2, χ2_red = finite_time_scaling2(K_vals, t_vals, p2_mat, p2_err_mat, X_data, Y_data, Yerr_data, D, V_guess; transient=transient)
#perform_collapse2(K_vals, X_data, Y_data, Yerr_data, shifts_data, χ2, χ2_red, pbest, perr, D, a_s; showxi=true)

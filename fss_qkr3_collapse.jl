include("imp_data_ex.jl")
include("data_sampling.jl")


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

#sampling data
K_vals_sampled, p2_sampled, p2_err_sampled = finite_time_scaling_sampling(K_vals, t_vals, p2_mat, p2_err_mat; critic_estimate=false, Kc_input=1.1, N_new=5, ΔK=0.001)

#Scaling data
X_data, Y_data, Yerr_data = finite_time_scaling_data(t_vals, p2_mat, p2_err_mat; d=D, n_kicks_i=transient, n_kicks_f=0)
X_data_sampled, Y_data_sampled, Yerr_data_sampled = finite_time_scaling_data(t_vals, p2_sampled, p2_err_sampled; d=D, n_kicks_i=transient, n_kicks_f=0)


#Scaling variance optimizaztion
shifts_data, s_rel_data = finite_time_scaling(X_data, Y_data; nbins=50)
_, shiftserr_data, _ = shifts_parametric_mc(t_vals, p2_mat, p2_err_mat; d=D, nbins=30, nmc=1000)

#collapse
perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data; data_type=data_type, raw=false, d=D, save=false)
perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data, D, a_s, data_type; Kc_offset=2, plotshow=false)

#sampled collapse
shifts_data_sampled, s_rel_data_sampled = finite_time_scaling(X_data_sampled, Y_data_sampled; nbins=50)
_, shiftserr_data_sampled, _ = shifts_parametric_mc(t_vals, p2_sampled, p2_err_sampled; d=D, nbins=30, nmc=1000)
perform_collapse(K_vals_sampled, X_data_sampled, Y_data_sampled, Yerr_data_sampled, shifts_data_sampled; data_type="Ex_sampled", raw=false, d=D, save=false)
perform_collapse_quality(K_vals_sampled, X_data_sampled, Y_data_sampled, Yerr_data_sampled, shifts_data_sampled, s_rel_data_sampled, D, a_s, "Ex_sampled"; Kc_offset=2, plotshow=false)
#monte_carlo_sampling(K_vals, t_vals, p2_mat, p2_err_mat; d=D, n_kicks_i=transient, nbins=30, nmc=1000)

#Scalling collapse Taylor fitting
#V_guess = [1, 0.1, 0.8, 1, -20, 10] #(b1, b2, Kc, ν, F00, ξsat)
#pbest, perr, shifts_data, χ2, χ2_red = finite_time_scaling2(K_vals, t_vals, p2_mat, p2_err_mat, X_data, Y_data, Yerr_data, D, V_guess; transient=transient)
#perform_collapse2(K_vals, X_data, Y_data, Yerr_data, shifts_data, χ2, χ2_red, pbest, perr, D, a_s; showxi=true)


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
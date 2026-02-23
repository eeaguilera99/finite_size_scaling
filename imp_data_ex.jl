using CSV, DataFrames, Plots, LaTeXStrings

a_s = 220
K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_matrix.csv", DataFrame; header=false))   
p2_err_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_err_matrix.csv", DataFrame; header=false))     # Errors

dim = 3 # spatial dimension
transient = 5
res_data, shifts_data, X_data, Y_data, Yerr_data, s_rel_data = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=transient, n_kicks_f=0)
shiftserr = shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, nbins=30, nmc=1000)[2]
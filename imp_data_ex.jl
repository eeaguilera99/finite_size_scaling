using CSV, DataFrames, Plots, LaTeXStrings

a_s = 220
K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_matrix.csv", DataFrame; header=false))   
p2_err_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_err_matrix.csv", DataFrame; header=false))     # Errors
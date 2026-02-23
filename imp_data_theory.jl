using CSV, DataFrames, Plots, LaTeXStrings


a_s = 220
K_vals = vec(Matrix(CSV.read("dataMF/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals_0 = vec(Matrix(CSV.read("dataMF/d=5_horizontal_axis.csv", DataFrame; header=false)))  # Times
p2_mat_0 = Matrix(CSV.read("dataMF/d=3_scaled_kinetic_energy_$(a_s).csv", DataFrame; header=false))
nc_mat = Matrix(CSV.read("dataMF/d=3_scaled_nc_2_$(a_s).csv", DataFrame; header=false))              # ⟨p²⟩ values

"Theory values of time are scaled, we revert them for dimension d1
For p2 values, the matrix is scaled but also rows are t values and 
columns are k values, we revert and transpose for dimension d2"
function revert_scale(time_vals, p2_vals, nc2_vals, d1, d2)
    t_vals = exp.(time_vals .* -d1)
    p2_mat = exp.(p2_vals) .* (t_vals .^ (2/d2))
    nc_mat = exp.(nc2_vals) .* (t_vals .^ (2/d2))
    return t_vals, Matrix(p2_mat'), Matrix(nc_mat')
end 
dim1 = 5
dim2 = 3
t_vals, p2_mat, nc_mat = revert_scale(t_vals_0, p2_mat_0, nc_mat_0, dim1, dim2)
p2_err_mat = 0.01 .* p2_mat  # assume 1% error if no data
nc_err_mat = 0.01 .* nc_mat 
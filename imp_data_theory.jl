include("fss_qkr3_timescaling.jl")
using CSV, DataFrames, Plots, LaTeXStrings


a_s = 920
K_vals_data = vec(Matrix(CSV.read("dataMF/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals_0 = vec(Matrix(CSV.read("dataMF/d=5_horizontal_axis.csv", DataFrame; header=false)))  # Times
p2_mat_0 = Matrix(CSV.read("dataMF/d=3_scaled_kinetic_energy_$(a_s).csv", DataFrame; header=false))
nc_mat_0 = Matrix(CSV.read("dataMF/d=3_scaled_nc_2_$(a_s).csv", DataFrame; header=false))              # ⟨p²⟩ values

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
t_vals_data, p2_mat_data, nc_mat_data = revert_scale(t_vals_0, p2_mat_0, nc_mat_0, dim1, dim2)
p2_err_mat = 0.01 .* p2_mat_data  # assume 1% error if no data
nc_err_mat = 0.01 .* nc_mat_data 

avg = 0
t_transient = 1
#Smothening 

#Adaptive moving avg
amp = 6
#p2_mat_avg = apply_mov_av_matrix(p2_mat_data, p=avg, loc_amp=amp)
#corr_p2 = trend_rep_avg(p2_mat_data, p2_mat_avg)
nc_mat_avg = apply_mov_av_matrix(nc_mat_data, p=avg, loc_amp=amp)
#corr_nc = trend_rep_avg(nc_mat_data, nc_mat_avg)

#Perform FFT smoothing
"""freq cutoff defines a ratio of lowfreq to retain, smaller ratio means more smoothing"""
#smooth_factor = 0.1
#p2_mat_avg = apply_lowpass_fft_matrix(p2_mat, smooth_factor, p=avg)
#nc_mat_avg = apply_lowpass_fft_matrix(nc_mat, smooth_factor, p=avg)

#Scaling analysis
d1 = 3
d2 = 3 
#res1, shifts1, X1, Y1, Yerr1, s_rel1 = finite_time_scaling(K_vals_data, t_vals_data, p2_mat_avg, p2_err_mat; d=d1, n_kicks_i=t_transient, n_kicks_f=0)
#shifts1_err = shifts_parametric_mc(K_vals_data, t_vals_data, p2_mat_avg, p2_err_mat; d=d1, nbins=30, nmc=1000)[2]
res2, shifts2, X2, Y2, Yerr2, s_rel2 = finite_time_scaling(K_vals_data, t_vals_data, nc_mat_avg, nc_err_mat; d=d2, n_kicks_i=t_transient, n_kicks_f=0)
#shifts2_err = shifts_parametric_mc(K_vals_data, t_vals_data, nc_mat_avg, nc_err_mat; d=d2, nbins=30, nmc=1000)[2]
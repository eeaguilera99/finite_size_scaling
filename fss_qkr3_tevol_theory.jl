include("fss_qkr3_timescaling.jl")
include("imp_data_theory.jl")
include("fss_qrk3_timescaling_analysis.jl")

dim = 3
Kc_offset_index = 0 #index to shift critical K by, if needed

#perform_tevol(dim, shifts1, t_vals, p2_mat_avg, p2_err_mat; Kc_offset_index=Kc_offset_index, data_type="Theory \$E_k\$")
perform_tevol(dim, shifts2, t_vals, nc_mat_avg, nc_err_mat; Kc_offset_index=Kc_offset_index, data_type="Theory \$1/n_c^2\$")
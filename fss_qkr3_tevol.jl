include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("imp_data_ex.jl")  # for finite_time_scaling on experimental data
include("fss_qrk3_timescaling_analysis.jl")

"""
Plot kinetic energy as function of time at criticality for exponent of difussion.
"""
d1 = 3
Kc_offset_index = 0 #index to shift critical K by, if needed


perform_tevol(d1, shifts_data, t_vals, p2_mat, p2_err_mat; Kc_offset_index=Kc_offset_index, data_type="Ex \$1/n_c^2\$")


include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
#include("imp_data_ex.jl")  # for finite_time_scaling on experimental data
#include("fss_qrk3_Kc.jl")  # for Kc analysis
include("fss_qrk3_timescaling_analysis.jl")

"""
Plot kinetic energy as function of time at criticality for exponent of difussion.
Note: Uses Kc from fss_qrk3_Kc.jl
"""

#perform_tevol(dim, shifts_data, t_vals, p2_mat, p2_err_mat; Kc_offset_index=Kc_offset_index, data_type="Ex \$1/n_c^2\$")

#perform_tevol_transient(d1, shifts_data, t_vals, p2_mat, p2_err_mat; Kc_offset_index=Kc_offset_index, data_type="Ex \$1/n_c^2\$",
#                            n_kicks_i=transient, n_kicks_f=0)
K_1 = argmax(shifts_data) + Kc_offset_index
K_2 = 16

plt1 = scatter(t_vals, p2_mat[K_2, :], yerror=p2_err_mat[argmax(shifts_data), :], 
label=latexstring("Data \$κ_c=$(round(K_vals[argmax(shifts_data)+Kc_offset_index], digits=3))\$"), xlabel="Time (N kicks)", 
ylabel="⟨p²⟩", title=latexstring("Time evolution at criticality \$d=$dim\$, $data_type, \$a_s=$a_s\$"), legend=:topleft, xscale=:log10, yscale=:log10)
display(plt1)
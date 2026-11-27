#include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("fss_qkr3_timescaling_analysis.jl") 




Kc_offset_index = -1
 # index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit


Kc_1 = perform_Kc_anal(shifts_data, shiftserr_data, K_vals; data_type="Ex", d=D, n_k_filter=fit_k_filter, K_guess_index=Kc_offset_index, show_xierr=true, save_data=true)
#Kc_1_sampled = perform_Kc_anal(shifts_data_sampled, shiftserr_data_sampled, K_vals_sampled; data_type="Ex_sampled", d=D, n_k_filter=fit_k_filter, K_guess_index=Kc_offset_index, show_xierr=false, save=false)[4]

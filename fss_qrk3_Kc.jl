#include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("fss_qrk3_timescaling_analysis.jl") 



Kc_offset_index = 1
 # index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit


Kc_1 = perform_Kc_anal(shifts_data, shiftserr_data, K_vals; data_type="Ex", d=D, n_k_filter=fit_k_filter, K_guess_index=Kc_offset_index, show_xierr=false, save=true)[4]


#=#save data
d1 = DataFrame(xi', :auto)
CSV.write("ξ(k)_data_220.csv", d1)

d2 = DataFrame(ξfit', :auto)
CSV.write("ξ(k)_data_fit_220.csv", d2)

d3 = DataFrame(Kgrid', :auto)
CSV.write("Kvals_ξfit_220.csv", d3)

d4 = DataFrame(xierr', :auto)
CSV.write("ξerr(k)_data_220.csv", d4)=#
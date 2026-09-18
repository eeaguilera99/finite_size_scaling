##
include("imp_data_ex.jl")
include("fss_qkr3_data_sampling.jl")


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

#filter transient data
t_vals, p2_mat, p2_err_mat = filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=transient)
""" Experimental data has a transients from where kicking takes effect. We firlter out this data to analyse new simulated data for
earlier times"""

#sampling data
t_vals_new = generate_sampling_t(t_vals; t_early=200, Δt=25)
K_vals_sampled, t_vals_sampled, p2_sampled, p2_err_sampled = finite_time_scaling_sampling(K_vals, t_vals, p2_mat, p2_err_mat; Kc_input=1.191, N_new=10, ΔK=0.06, t_new=t_vals_new)

#Scaling data
X_data, Y_data, Yerr_data = finite_time_scaling_data(t_vals, p2_mat, p2_err_mat; d=D, n_kicks_i=1)
X_data_sampled, Y_data_sampled, Yerr_data_sampled = finite_time_scaling_data(t_vals_sampled, p2_sampled, p2_err_sampled; d=D, n_kicks_i=1)


#Scaling variance optimizaztion
shifts_data, s_rel_data = finite_time_scaling(X_data, Y_data; nbins=50)
#_, shiftserr_data, _ = shifts_parametric_mc(t_vals, p2_mat, p2_err_mat; d=D, nbins=30, nmc=1000)

#collapse
#perform_collapse(K_vals, X_data, Y_data, Yerr_data, shifts_data; data_type=data_type, raw=false, d=D, save=false)
#perform_collapse_quality(K_vals, X_data, Y_data, Yerr_data, shifts_data, s_rel_data, D, a_s, data_type; Kc_offset=2, plotshow=false)

#sampled collapse

shifts_data_sampled, s_rel_data_sampled = finite_time_scaling(X_data_sampled, Y_data_sampled; nbins=50)
_, shiftserr_data_sampled, _ = shifts_parametric_mc(t_vals_sampled, p2_sampled, p2_err_sampled; d=D, nbins=30, nmc=1000)

perform_collapse(K_vals_sampled, X_data_sampled, Y_data_sampled, Yerr_data_sampled, shifts_data_sampled; data_type="Ex sampled", plot_label_b=false, raw=false, d=D, save=false)
perform_collapse_quality(K_vals_sampled, X_data_sampled, Y_data_sampled, Yerr_data_sampled, shifts_data_sampled, s_rel_data_sampled, D, a_s, "Ex sampled"; Kc_offset=2, plotshow=false)
Kc_offset_index = 0
 # index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit
Kc_1_sampled = perform_Kc_anal(shifts_data_sampled, shiftserr_data_sampled, K_vals_sampled; data_type="Ex sampled", d=D, n_k_filter=fit_k_filter, K_guess_index=Kc_offset_index, show_xierr=false, save=false)[4]
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

##
#compare original and sampled data
plt1 =  plot(title=latexstring("Og data, \$d=$(D)\$"), xlabel=L"\ln N^{1/d}", ylabel=L"\ln Λ")
plot!(plt1, repeat(X_data, length(K_vals),1)', Y_data', label="")
scatter!(plt1, repeat(X_data, length(K_vals),1)', Y_data', label="") #Plot orginal data not scaled
display(plt1)
in_new = indexin(setdiff(K_vals_sampled, K_vals), K_vals_sampled) #get index of new K values in sampled data
Y_data_sampled_new = Y_data_sampled[in_new, :] #get Y values for new K values
plt2 =  plot(title=latexstring("New simulated data, \$d=$(D)\$"), xlabel=L"\ln N^{1/d}", ylabel=L"\ln Λ")
plot!(plt2, repeat(X_data_sampled, length(in_new),1)', Y_data_sampled_new', label="") #plot sampled data not scaled
scatter!(plt2, repeat(X_data_sampled, length(in_new),1)', Y_data_sampled_new', label="") #plot sampled data not scaled
xlims!(plt2, xlims(plt1))
ylims!(plt2, ylims(plt1))
display(plt2)
plt3 =  plot(title=latexstring("Localization length, \$d=$(D)\$"), xlabel=L"κ", ylabel=L"ξ")
scatter!(plt3, K_vals, exp.(shifts_data), ms=6, label="")
scatter!(plt3, K_vals_sampled[in_new], exp.(shifts_data_sampled[in_new]), ms=6, label="Sampled data", color=:red)
display(plt3)
##
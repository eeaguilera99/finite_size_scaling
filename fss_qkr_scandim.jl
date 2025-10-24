include("fss_qkr3_collapse.jl")  # for finite_time_scaling

#code to loop over dimension values for best collapse

d_vals = 2:0.1:3

for d in d_values
    shifts, X, Y, min_val = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=d)
    push!(collapse_quality, min_val)
    shift_dict[d] = shifts
    println("  d = $(round(d, digits=2)) → collapse variance = $(round(min_val, sigdigits=4))")
end
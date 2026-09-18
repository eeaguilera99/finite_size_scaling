include("fss_qkr3_timescaling_theory.jl")
using DSP

# Observable: Λ = <p^2>/t^(2/3)
function scaling_funct(mat, err_mat)
    Λ = mat ./ (t_vals' .^ (2/dim))
    Λ_err = err_mat ./ (t_vals' .^ (2/dim))
    return log.(Λ), log.(Λ_err ./ Λ)
end

function moving_average(data, window_size)
    return [mean(data[i:i+window_size-1]) for i in 1:(length(data)-window_size+1)]
end

dim = 3
nc_mat_avg = apply_mov_av_matrix(nc_mat, p=true, loc_amp=6)
nc_mat_avg2 = hcat([moving_average(nc_mat[i, :], 3) for i in axes(nc_mat,1)]...)'
nc_mat_avg3 = hcat(sgolayfilt(nc_mat[i, :], SGFilter(11, 3)) for i in axes(nc_mat,1))'
corr_nc = trend_rep_avg(nc_mat, nc_mat_avg)
t_transient = 28

# Log variables
X = -log.(t_vals' .^ (1/dim))     # 1×N
Y, Yerr = scaling_funct(nc_mat, nc_err_mat)
Y_avg, _ = scaling_funct(nc_mat_avg, nc_err_mat)


#plots

#plot individual curve 
k_index = 2
plt1 = plot(title="Averaging for κ=$(round(K_vals[k_index], digits=3))", xlabel="ln(t^(-1/d))", ylabel="ln(Λ)")
plot!(plt1, X[:], Y[k_index,:] , marker=:o, label="og")
plot!(plt1, X[:], Y_avg[k_index, :], marker=:o, label="avg")
display(plt1)
#=
#plot correlations
plt2 = plot(bar(K_vals, corr_nc), xlabel="κ", ylabel="Corr(nc, nc_avg)", ylims=(-1,1))
display(plt2)=#


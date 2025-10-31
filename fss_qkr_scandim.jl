include("fss_qkr3_timescaling.jl")  # for finite_time_scaling

#code to loop over dimension values for best collapse

d_values = 1.7:0.05:10

collapse_quality_1 = Float64[]
#shift_dict = Dict{Float64, Vector{Float64}}()

for dim in d_values
    l_res, l_shifts, l_X, l_Y, l_Y_err, sX = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=5)
    push!(collapse_quality_1, sX)
end

function dim_scan(collapse_quality)
    best_idx = argmin(collapse_quality)
    best_d = d_values[best_idx]
    best_val = collapse_quality[best_idx]


    println("\n✅ Best collapse found for:")
    println("   d_best = ", best_d)
    println("   min_val = ", best_val)

    # Plot collapse quality vs dimension
    plt_quality = plot(d_values, collapse_quality, lw=1, marker=:o,
        xlabel="Dimension d", ylabel=L"\sigma^{rel}",
        title="Quality of scaling collapse vs dimension d", label="")
    scatter!(plt_quality, [best_d], [best_val], label="Best d = $(round(best_d, digits=2))", markersize=8)
    display(plt_quality)
end

dim_scan(collapse_quality_1)

#=
# Optional: Plot best collapse curves
plt_best = plot(title="Best Data Collapse (d = $(round(best_d, digits=2)))",
    xlabel="ln(ξ/t^(1/d))", ylabel="ln(Λ)", legend=:outertopright)
for (i, K) in enumerate(K_vals)
    plot!(plt_best, X[:] .+ best_shifts[i], Y[i, :], marker=:o, label="K=$K")
end
display(plt_best)

=#
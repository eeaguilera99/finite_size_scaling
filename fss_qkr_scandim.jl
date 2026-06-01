include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("fss_qrk3_timescaling_analysis.jl")

#code to loop over dimension values for best collapse

d_vals = 1:0.5:10


#=
function dim_scan(d_values)

    collapse_quality = Float64[]
    for dim in d_values
        _, _, _, _, sX = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=50, d=dim, n_kicks_i=2)
        push!(collapse_quality, sX)
    end
    best_idx = argmin(collapse_quality)
    best_d = d_values[best_idx]
    best_val = collapse_quality[best_idx]


    println("\n✅ Best collapse found for:")
    println("   d_best = ", best_d)
    println("   min_val = ", best_val)

    # Plot collapse quality vs dimension
    #plotly()
    plt_quality = plot(d_values, collapse_quality, lw=1, marker=:o,
        xlabel="Dimension d", ylabel=L"\sigma^{rel}",
        title="Quality of scaling collapse vs dimension d", label="")
    scatter!(plt_quality, [best_d], [best_val], label="Best d = $(round(best_d, digits=2))", markersize=8)
    display(plt_quality)
end

dim_scan(d_vals)=#

function Kc_scan(d_values)
    Kc_values = Float64[]
    for dim in d_values
        shifts, _, _, _, _ = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=50, d=dim, n_kicks_i=2)
        _, shiftserr, _ = shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, nbins=30, nmc=1000)
        Kc = perform_Kc_anal(shifts, shiftserr, K_vals; data_type="Ex", d=dim, n_k_filter=0, K_guess_index=0, show_xierr=false, save=false)[4]
        push!(Kc_values, Kc)
    end

    # Plot Kc vs dimension
    plt_Kc = plot(d_values, Kc_values, lw=1, marker=:o,
        xlabel="Dimension d", ylabel=L"K_c",
        title="Critical kick strength Kc vs dimension d", label="")
    display(plt_Kc)
end
Kc_scan(d_vals)

#=
# Optional: Plot best collapse curves
plt_best = plot(title="Best Data Collapse (d = $(round(best_d, digits=2)))",
    xlabel="ln(ξ/t^(1/d))", ylabel="ln(Λ)", legend=:outertopright)
for (i, K) in enumerate(K_vals)
    plot!(plt_best, X[:] .+ best_shifts[i], Y[i, :], marker=:o, label="K=$K")
end
display(plt_best)

=#
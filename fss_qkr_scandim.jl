include("fss_qkr3_timescaling.jl")
include("fss_qrk3_timescaling_analysis.jl")
using CSV, DataFrames, Plots, LaTeXStrings


#a_s = "κ=0.87"
a_s = 220
K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
#K_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/a_s.csv", DataFrame; header=false)))             # interactions
t_vals = vec(Matrix(CSV.read("dataEX/$(a_s)/Number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_matrix.csv", DataFrame; header=false))   
p2_err_mat = Matrix(CSV.read("dataEX/$(a_s)/nc_err_matrix.csv", DataFrame; header=false))     
#code to loop over dimension values for best collapse

d_vals = 1:0.5:10
V_guess = [1, 1, 1.2, 0.5, -20, 0.1]
transient = 5


function dim_scan(d_values)

    collapse_quality = Float64[]
    for dim in d_values
        _, sX = finite_time_scaling(X_data, Y_data; nbins=50, d=dim, n_kicks_i=transient)
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

#dim_scan(d_vals)

function dim_scan2(d_vals)
    collapse_quality = Float64[]
    for dim in d_vals
        X_data, Y_data, Yerr_data = finite_time_scaling_data(t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=transient, n_kicks_f=0)
        _, _, _, _, χ = finite_time_scaling2(K_vals, t_vals, p2_mat, p2_err_mat, X_data, Y_data, Yerr_data, dim, V_guess; transient=transient)
        push!(collapse_quality, χ)
    end

    best_idx = argmin(collapse_quality)
    best_d = d_vals[best_idx]
    best_val = collapse_quality[best_idx]

    println("\n✅ Best collapse found for:")
    println("   d_best = ", best_d)
    println("   min_val = ", best_val)
    println(" $d_vals , $(round.(collapse_quality, digits=2)) ")
    println("Initial guess = $V_guess (b1, b2, Kc, ν, F00, ξsat)")

    # Plot collapse quality vs dimension
    #plotly()
    plt_quality = plot(d_vals, collapse_quality, lw=1, marker=:o,
        xlabel="Dimension \$d\$", ylabel=L"χ2_{rel}",
        title="Quality of scaling2 collapse \$vs\$ dimension \$d\$", label="")
    scatter!(plt_quality, [best_d], [best_val], label="Best d = $(round(best_d, digits=2))", markersize=8)
    display(plt_quality)
end

dim_scan2(d_vals)

function Kc_scan(d_values)
    Kc_values = Float64[]
    for dim in d_values
        shifts, _ = finite_time_scaling(X_data, Y_data; nbins=50, d=dim, n_kicks_i=transient)
        _, shiftserr, _ = shifts_parametric_mc(t_vals, p2_mat, p2_err_mat; d=dim, nbins=30, nmc=1000)
        Kc = perform_Kc_anal(shifts, shiftserr, K_vals; data_type="Ex", d=dim, n_k_filter=0, K_guess_index=0, show_xierr=false, save=false)[4]
        push!(Kc_values, Kc)
    end

    # Plot Kc vs dimension
    plt_Kc = plot(d_values, Kc_values, lw=1, marker=:o,
        xlabel="Dimension d", ylabel=L"K_c",
        title="Critical kick strength Kc vs dimension d", label="")
    display(plt_Kc)
end
#Kc_scan(d_vals)

#=
# Optional: Plot best collapse curves
plt_best = plot(title="Best Data Collapse (d = $(round(best_d, digits=2)))",
    xlabel="ln(ξ/t^(1/d))", ylabel="ln(Λ)", legend=:outertopright)
for (i, K) in enumerate(K_vals)
    plot!(plt_best, X[:] .+ best_shifts[i], Y[i, :], marker=:o, label="K=$K")
end
display(plt_best)

=#
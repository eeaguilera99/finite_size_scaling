using CSV, DataFrames, Optim, Statistics, Plots

# 1. Load data
K_vals = vec(Matrix(CSV.read("Codes\\finite_size_scaling\\data\\kappa.csv", DataFrame)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("Codes\\finite_size_scaling\\data\\number_of_kicks.csv", DataFrame)))  # Times
p2_mat = Matrix(CSV.read("Codes\\finite_size_scaling\\data\\nc_matrix.csv", DataFrame))             # ⟨p²⟩ values
p2_err_mat = Matrix(CSV.read("Codes\\finite_size_scaling\\data\\nc_err_matrix.csv", DataFrame))     # Errors

# 2. Reshape into long DataFrame
function make_long_dataframe(K_vals, t_vals, p2_mat, p2_err_mat)
    rows = Vector{NamedTuple}()
    for (i, K) in enumerate(K_vals)
        for (j, t) in enumerate(t_vals)
            push!(rows, (K=K, t=t, p2=p2_mat[i,j], p2_err=p2_err_mat[i,j]))
        end
    end
    return DataFrame(rows)
end

data = make_long_dataframe(K_vals, t_vals, p2_mat, p2_err_mat)
 

# 3. Scaling transforms
function scaling_transform(p2, t)
    return p2 / t^(2/3)
end

function scaling_variable(K, Kc, t)
    return (K - Kc) * t^(1/3)
end

# 4. Collapse error function
function collapse_error(Kc, data)
    xvals, yvals = Float64[], Float64[]
    for row in eachrow(data)
        if row.t > 0 && row.p2 > 0  # avoid invalid points
            push!(xvals, scaling_variable(row.K, Kc, row.t))
            push!(yvals, scaling_transform(row.p2, row.t))
        end
    end
    # binning
    bins = collect(-5:0.3:5)
    errors = Float64[]
    for i in 1:length(bins)-1
        idx = findall(x -> bins[i] ≤ x < bins[i+1], xvals)
        if length(idx) > 3
            push!(errors, var(yvals[idx]))
        end
    end
    return mean(errors)
end

# 5. Optimize Kc
res = optimize(Kc -> collapse_error(Kc, data), 0.5, 2.0)  # adjust interval as needed
Kc_best = Optim.minimizer(res)
println("Best Kc ≈ ", Kc_best)

# 6. Plot collapsed data
function plot_collapse(data, Kc)
    xvals, yvals = Float64[], Float64[]
    for row in eachrow(data)
        if row.t > 0 && row.p2 > 0
            push!(xvals, scaling_variable(row.K, Kc, row.t))
            push!(yvals, scaling_transform(row.p2, row.t))
        end
    end
    scatter(xvals, yvals, alpha=0.4, label="", ms=3)
    xlabel!("(K - Kc) * t^(1/3)")
    ylabel!("Log⟨p²(t)⟩ / t^(2/3)")
    title!("Finite-time scaling collapse (Kc = $(round(Kc, digits=3)))")
end

#plot_collapse(data, Kc_best)


# === 7. Plot (like Fig. 10) ===
function plot_finite_time_scaling(data)
    plt = plot(xlabel="ln(1 / t^(1/3))", ylabel="ln(Λ)",
               title="Finite-time scaling (raw, like Fig. 10 left)",
               legend=:right)

    for Kval in unique(data.K)
        subset = filter(row -> row.K == Kval && row.t > 0 && row.p2 > 0, data)
        Λ = [scaling_transform(row.p2, row.t) for row in eachrow(subset)]
        xvals = log.(1 ./ (subset.t .^ (1/3)))
        yvals = log.(Λ)
        plot!(plt, xvals, yvals, label="K=$(round(Kval, digits=2))", lw=1.5)
    end
    return plt
end

finite_time_plot = plot_finite_time_scaling(data)
display(finite_time_plot)

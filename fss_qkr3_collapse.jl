using LinearAlgebra
using Statistics
using Optim
using Plots
using LsqFit
using CSV, DataFrames

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
K_vals = vec(Matrix(CSV.read("data/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("data/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("data/nc_matrix.csv", DataFrame; header=false))             # ⟨p²⟩ values
p2_err_mat = Matrix(CSV.read("data/nc_err_matrix.csv", DataFrame; header=false))     # Errors

d = 3  # spatial dimension

#filter Nkicks range
n_Nkicks_i = 4 #index to start from
n_Nkicks_f = size(t_vals,1)-5 #index to end at
t_vals = t_vals[n_Nkicks_i:n_Nkicks_f]
p2_mat = p2_mat[:,n_Nkicks_i:n_Nkicks_f]
p2_err_mat = p2_err_mat[:,n_Nkicks_i:n_Nkicks_f]


function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=30)
    M, N = size(p2_mat)

    # Observable: Λ = <p^2>/t^(2/3)
    Λ = p2_mat ./ (t_vals' .^ (2/d))
    Λ_err = p2_err_mat ./ (t_vals' .^ (2/d))

    # Log variables
    X = -log.(t_vals') ./ 3    # 1×N
    Y = log.(Λ)                # M×N
    # Propagate errors: Δ(ln Λ) ≈ ΔΛ / Λ
    Yerr = Λ_err ./ Λ

    # Flatten for binning
    allY = vec(Y)
    y_min, y_max = minimum(allY), maximum(allY)
    bins = range(y_min, y_max; length=nbins+1)
    bin_ids = [searchsortedlast(bins, y) for y in allY]

    # Cost function: variance of shifted X within Y-bins
    function cost(a_full::Vector)
        shiftedX = vec(X .+ a_full .* ones(1,N))
        total_var = 0.0
        for b in 1:nbins
            mask = (bin_ids .== b)
            if count(mask) > 1
                xb = shiftedX[mask]
                total_var += var(xb) * count(mask)
            end
        end
        return total_var
    end

    # Fix gauge: a₁ = 0
    function constrained_cost(a_free::Vector)
        a_full = vcat(0.0, a_free)
        return cost(a_full)
    end

    # Initial guess
    a0 = zeros(M-1)

    # Minimize
    res = optimize(constrained_cost, a0, NelderMead())
    shifts = vcat(0.0, Optim.minimizer(res))
    println(string(res.minimum, " = minimum cost with shifts "))
    return res, shifts, X, Y, Yerr
end


# Perform collapse
res, shifts, X, Y, Yerr = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat)

# Plot before collapse
plt1 = plot(title="Raw data (before shifts)",
    xlabel="ln(t^(-1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt1, X[:], Y[i,:], yerror=Yerr[i,:], marker=:o, label="K=$K")
end
display(plt1)

# Plot after collapse
plt2 = plot(title="Data collapse (after optimal shifts)",
    xlabel="ln(ξ/t^(1/3))", ylabel="ln(Λ)")
for (i,K) in enumerate(K_vals)
    plot!(plt2, X[:] .+ shifts[i], Y[i,:], yerror=Yerr[i,:], marker=:o, label="K=$K")
end
display(plt2)

println("Optimal shifts ln ξ(K): ", shifts)
println("Fit quality: ", res.minimum)


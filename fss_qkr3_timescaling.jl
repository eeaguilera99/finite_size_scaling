using LinearAlgebra
using Statistics
using Optim
using Plots
using LsqFit
using CSV, DataFrames

K_vals = vec(Matrix(CSV.read("data2/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("data2/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("data2/nc_matrix.csv", DataFrame; header=false))             # ⟨p²⟩ values
p2_err_mat = Matrix(CSV.read("data2/nc_err_matrix.csv", DataFrame; header=false))     # Errors

#filter Nkicks range
n_Nkicks_i = 6 #index to start from
n_Nkicks_f = size(t_vals,1) #index to end at
t_vals = t_vals[n_Nkicks_i:n_Nkicks_f]
p2_mat = p2_mat[:,n_Nkicks_i:n_Nkicks_f]
p2_err_mat = p2_err_mat[:,n_Nkicks_i:n_Nkicks_f]

function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=30, d)
    
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
    return res, shifts, X, Y, Yerr
end
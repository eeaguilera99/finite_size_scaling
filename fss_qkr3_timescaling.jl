using LinearAlgebra
using Statistics, Random, Distributions
using Optim
using Plots
using LsqFit
using CSV, DataFrames
using LaTeXStrings
using ForwardDiff 

a_s = 220

K_vals = vec(Matrix(CSV.read("data2/$(a_s)/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals = vec(Matrix(CSV.read("data2/$(a_s)/number_of_kicks.csv", DataFrame; header=false)))  # Times
p2_mat = Matrix(CSV.read("data2/$(a_s)/nc_matrix.csv", DataFrame; header=false))             # ⟨p²⟩ values
p2_err_mat = Matrix(CSV.read("data2/$(a_s)/nc_err_matrix.csv", DataFrame; header=false))     # Errors



function filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=1, n_kicks_f=0)
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]
    return t_vals, p2_mat, p2_err_mat
end

function filter_K(Kk_vals, mat, err_mat; n_kkicks_i=1, n_kkicks_f=0)
    n_kkicks_ff = length(Kk_vals) - n_kkicks_f
    Kk_vals = Kk_vals[n_kkicks_i:n_kkicks_ff]
    mat = mat[n_kkicks_i:n_kkicks_ff,:]
    err_mat = err_mat[n_kkicks_i:n_kkicks_ff,:]
    return Kk_vals, mat, err_mat
end

function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=100, d=3, n_kicks_i=1, n_kicks_f=0)
    
    #filter Nkicks range
    t_vals, p2_mat, p2_err_mat = filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=n_kicks_i, n_kicks_f=n_kicks_f)

    M, N = size(p2_mat)

    # Observable: Λ = <p^2>/t^(2/3)
    Λ = p2_mat ./ (t_vals' .^ (2/d))
    Λ_err = p2_err_mat ./ (t_vals' .^ (2/d))

    # Log variables
    X = -log.(t_vals' .^ (1/d))     # 1×N
    Y = log.(Λ)                # M×N
    # Propagate errors: Δ(ln Λ) ≈ ΔΛ / Λ
    Yerr = Λ_err ./ Λ

    # Flatten for binning
    allY = vec(Y)
    y_min, y_max = minimum(allY), maximum(allY) #range of data Yaxis
    bins = range(y_min, y_max; length=nbins+1)
    bin_ids = [searchsortedlast(bins, y) for y in allY] #assigns each y-value to a bin number from 1 to nbins

    # Cost function: variance of shifted X within Y-bins
    function cost(a_full::Vector)
        # Shift X values by the amount specified in a_full
        shiftedX = vec(X .+ a_full .* ones(1,N))
        
        total_var = 0.0
        for b in 1:nbins
            # Find all points that fall in bin b
            mask = (bin_ids .== b)
            
            # Only consider bins with more than 1 point
            if count(mask) > 1
                # Get X values for points in this bin
                xb = shiftedX[mask]
                # Add weighted variance of X values in this bin
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
    a_free_opt = Optim.minimizer(res)

    # === Compute normalized scatter directly from res.minimum ===
    total_points = M * N
    sX = sqrt(res.minimum / total_points)
    Xp = X .+ shifts                     # shifted X matrix
    sX_rel = sX / (maximum(Xp) - minimum(Xp) + eps())

    # === Return everything
    return res, shifts, X, Y, Yerr, sX_rel
end

function shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=3, nbins=100, nmc=500, rng=MersenneTwister(0))
    M, N = size(p2_mat)
    all_shifts = zeros(nmc, M)

    for m in 1:nmc
        # Sample synthetic dataset
        noise = rand!(rng, Normal(), similar(p2_mat)) .* p2_err_mat
        p2_syn = p2_mat .+ noise

        # Ensure positivity (log will be used downstream)
        p2_syn = max.(p2_syn, eps())
        _, shifts_syn, _, _, _, _, _ = finite_time_scaling(K_vals, t_vals, p2_syn, p2_err_mat; d=d, nbins=nbins)
        all_shifts[m, :] .= shifts_syn
    end

    mean_shifts = vec(mean(all_shifts, dims=1))
    std_shifts  = vec(std(all_shifts, dims=1))
    return mean_shifts, std_shifts, all_shifts
end

function tot_variance(a_full::Vector, X::Matrix, Y::Matrix; nbins=100)# calculates rel var for arbitrary shifts
    Xp = vec(X .+ a_full .* ones(1,size(X,2)))
    
    # Flatten for binning
    allY = vec(Y)
    y_min, y_max = minimum(allY), maximum(allY) #range of data Yaxis
    bins = range(y_min, y_max; length=nbins+1)
    bin_ids = [searchsortedlast(bins, y) for y in allY]

    totw, totvar = 0.0, 0.0
    for b in 1:nbins
        mask = (bin_ids .== b)
        nb = count(mask)
        if nb > 1
            xb = vec(Xp)[mask]
            totvar += var(xb) * nb
            totw   += nb
        end
    end
    sX = sqrt(totvar / max(totw, 1.0))
    sX_rel = sX / (maximum(vec(Xp)) - minimum(vec(Xp)) + eps())

    return totvar, sX_rel
end
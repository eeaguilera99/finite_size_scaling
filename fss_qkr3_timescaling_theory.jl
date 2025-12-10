using LinearAlgebra
using Statistics
using Optim
using Plots
using LsqFit
using CSV, DataFrames
using LaTeXStrings

K_vals = vec(Matrix(CSV.read("dataMF/kappa.csv", DataFrame; header=false)))             # Kick strengths
t_vals_0 = vec(Matrix(CSV.read("dataMF/d=5_horizontal_axis.csv", DataFrame; header=false)))  # Times
p2_mat_0 = Matrix(CSV.read("dataMF/d=3_scaled_kinetic_energy_1038.csv", DataFrame; header=false))
nc_mat_0 = Matrix(CSV.read("dataMF/d=3_scaled_nc_2_1038.csv", DataFrame; header=false))              # ⟨p²⟩ values
#p2_err_mat = Matrix(CSV.read("data2/nc_err_matrix.csv", DataFrame; header=false))     # Errors
a_s = 1038

"Theory values of time are scaled, we revert them for dimension d1
For p2 values, the matrix is scaled but also rows are t values and columns are k values, we revert and transpose for dimension d2"

function revert_scale(time_vals, p2_vals, nc2_vals, d1, d2)
    t_vals = exp.(time_vals .* -d1)
    p2_mat = exp.(p2_vals) .* (t_vals .^ (2/d2))
    nc_mat = exp.(nc2_vals) .* (t_vals .^ (2/d2))
    return t_vals, p2_mat', nc_mat'
end
dim1 = 5
dim2 = 3
t_vals, p2_mat, nc_mat = revert_scale(t_vals_0, p2_mat_0, nc_mat_0, dim1, dim2)

p2_err_mat = 0.01 .* p2_mat  # assume 1% error if no data
nc_err_mat = 0.01 .* nc_mat 

#functon to apply moving average smoothing
function adaptive_moving_average(y; p=true, loc_amp=2, min_win=3, max_win=15, ε=1e-12)
    #p parameter input to enable/disable smoothing
    if p ==true
        N = length(y)
        smooth = similar(y)

        # robust scale: avoid global outlier domination
        global_scale = max(maximum(abs.(y)), ε)

        for i in 1:N
            # small probe window to estimate local amplitude (safe clamp)
            probe = max(1, i-loc_amp) : min(N, i+loc_amp)
            local_amp = maximum(y[probe]) - minimum(y[probe])

            # map local amplitude to window size (inverted: larger amp -> smaller window)
            frac = clamp(local_amp / global_scale, 0.0, 1.0)
            scaled_win = round(Int, max_win - frac * (max_win - min_win))

            # enforce bounds and oddness
            scaled_win = clamp(scaled_win, min_win, max_win)
            actual_win = isodd(scaled_win) ? scaled_win : scaled_win + 1
            actual_win = min(actual_win, max_win)              # ensure not exceed max

            hw = actual_win ÷ 2
            win_start = max(1, i - hw)
            win_stop  = min(N, i + hw)
            win = win_start:win_stop

            smooth[i] = mean(view(y, win))
        end

        return smooth
    else
        return y
    end
end

# Apply moving average to each row of a matrix
function apply_mov_av_matrix(M; p=true, loc_amp=2)
    for i in 1:size(M,1)
        M[i,:] = adaptive_moving_average(M[i,:]; p=p, loc_amp=loc_amp, min_win=3, max_win=15)
    end
    return M
end

function filter_Nkicks(tt_vals, mat, err_mat; n_kicks_i=1, n_kicks_f=0)
    n_Nkicks_f = size(tt_vals,1) - n_kicks_f #index to end at
    tt_vals = tt_vals[n_kicks_i:n_Nkicks_f]
    mat = mat[:,n_kicks_i:n_Nkicks_f]
    err_mat = err_mat[:,n_kicks_i:n_Nkicks_f]
    return tt_vals, mat, err_mat
end

function filter_K(Kk_vals, mat, err_mat; n_kkicks_i=1, n_kkicks_f=0)
    n_kkicks_ff = length(Kk_vals) - n_kkicks_f
    Kk_vals = Kk_vals[n_kkicks_i:n_kkicks_ff]
    mat = mat[n_kkicks_i:n_kkicks_ff,:]
    err_mat = err_mat[n_kkicks_i:n_kkicks_ff,:]
    return Kk_vals, mat, err_mat
end

function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=100, d, n_kicks_i=1, n_kicks_f=0)
    
    #filter Nkicks range
    n_Nkicks_f = size(t_vals,1) - n_kicks_f #index to end at
    t_vals = t_vals[n_kicks_i:n_Nkicks_f]
    p2_mat = p2_mat[:,n_kicks_i:n_Nkicks_f]
    p2_err_mat = p2_err_mat[:,n_kicks_i:n_Nkicks_f]

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

    # === Compute normalized scatter directly from res.minimum ===
    total_points = M * N
    sX = sqrt(res.minimum / total_points)
    Xp = X .+ shifts                     # shifted X matrix
    sX_rel = sX / (maximum(Xp) - minimum(Xp) + eps())

    # === Return everything
    return res, shifts, X, Y, Yerr, sX_rel
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
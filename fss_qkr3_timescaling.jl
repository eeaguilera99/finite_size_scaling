using LinearAlgebra
using Statistics, Random, Distributions
using Optim
using Plots
using LsqFit
using ForwardDiff
using FFTW 

#functon to apply moving average smoothing
function adaptive_moving_average(y; p=true, loc_amp=2, min_win=3, max_win=15, ε=1e-12)
    #p parameter input to enable/disable smoothing
    if p ==true
        N = length(y)
        smooth = similar(y)

        # robust scale: avoid global outlier domination
        global_scale = max(maximum(abs.(y)), ε)#prevent one big spike from dominating the amplitude scaling

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

            smooth[i] = mean(view(y, win))#computes the avg
        end

        return smooth
    else
        return y
    end
end

# Apply moving average to each row of a matrix
function apply_mov_av_matrix(M; p=true, loc_amp=2)
    M_avg = similar(M)
    for i in axes(M,1)
        M_avg[i,:] = adaptive_moving_average(M[i,:]; p=p, loc_amp=loc_amp, min_win=3, max_win=15)
    end
    return M_avg
end

#FFT for filtering
function lowpass_fft(y::Vector, cutoff_ratio::Float64)
    N = length(y)
    Y = fft(y)

    # Cutoff index in frequency domain
    cutoff = floor(Int, cutoff_ratio * N ÷ 2)

    # Zero high frequencies (keep 2*cutoff for symmetry)
    Y[cutoff+2:end-cutoff] .= 0

    # Inverse transform to get smoothed signal
    y_smooth = real(ifft(Y))
    return abs.(y_smooth)
end

# Apply lowpass FFT to each row of a matrix
function apply_lowpass_fft_matrix(M::Matrix, cutoff_ratio::Float64; p=true)
    if p != true
        return M
    else
        M_smooth = similar(M)
        for i in axes(M,1)
            M_smooth[i,:] = lowpass_fft(M[i,:], cutoff_ratio)
        end
        return M_smooth
    end
end

# Functions to exclude early/late time points (Nkicks) or Kkicks from analysis
function filter_Nkicks(tt_vals, mat, err_mat; n_kicks_i=1, n_kicks_f=0)
    n_Nkicks_f = size(tt_vals,1) - n_kicks_f #index to end at
    t_vals_filtered = tt_vals[n_kicks_i:n_Nkicks_f]
    mat_filtered = mat[:,n_kicks_i:n_Nkicks_f]
    err_mat_filtered = err_mat[:,n_kicks_i:n_Nkicks_f]
    return t_vals_filtered, mat_filtered, err_mat_filtered
end

function filter_K(Kk_vals, mat, err_mat; n_kkicks_i=1, n_kkicks_f=0)
    n_kkicks_ff = length(Kk_vals) - n_kkicks_f
    Kk_vals = Kk_vals[n_kkicks_i:n_kkicks_ff]
    mat = mat[n_kkicks_i:n_kkicks_ff,:]
    err_mat = err_mat[n_kkicks_i:n_kkicks_ff,:]
    return Kk_vals, mat, err_mat
end

function filter_transients_per_curve(t_vals, p2_mat, err_mat, n_kicks_i_vec; n_kicks_f_vec=zeros(Int,length(n_kicks_i_vec)))

    M = size(p2_mat,1)

    curves = []

    for i in 1:M

        start_i = n_kicks_i_vec[i]
        end_i   = length(t_vals) - n_kicks_f_vec[i]

        t_i   = t_vals[start_i:end_i]
        p2_i  = p2_mat[i,start_i:end_i]
        err_i = err_mat[i,start_i:end_i]

        push!(curves, (t=t_i, p2=p2_i, err=err_i))

    end

    return curves

end

function compute_scaling_data(curves; d=3)

    M = length(curves)

    X = Vector{Vector{Float64}}(undef,M)
    Y = Vector{Vector{Float64}}(undef,M)
    Yerr = Vector{Vector{Float64}}(undef,M)

    for i in 1:M

        t = curves[i].t
        p2 = curves[i].p2
        p2err = curves[i].err

        Λ = p2 ./ (t .^(2/d))
        Λerr = p2err ./ (t .^(2/d))

        X[i] = -log.(t) ./ d
        Y[i] = log.(Λ)
        Yerr[i] = Λerr ./ Λ

    end

    return X, Y, Yerr

end

function flatten_data(X,Y)

    X_all = Float64[]
    Y_all = Float64[]
    curve_id = Int[]

    for i in 1:length(X)
        append!(X_all, X[i])
        append!(Y_all, Y[i])
        append!(curve_id, fill(i,length(X[i])))
    end

    return X_all, Y_all, curve_id

end

#Main function for finite time scaling analysis
function finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat, n_kicks_i_vec; d=3, nbins=100)

    curves = filter_transients_per_curve(t_vals, p2_mat, p2_err_mat, n_kicks_i_vec)

    X, Y, Yerr = compute_scaling_data(curves, d=d)

    M = length(X)

    # Flatten data
    X_all, Y_all, curve_id = flatten_data(X, Y)

    # Bin in Y
    y_min, y_max = minimum(Y_all), maximum(Y_all)
    bins = range(y_min, y_max; length=nbins+1)

    bin_ids = [searchsortedlast(bins, y) for y in Y_all]

    # === Cost function ===
    function cost(a_full)

        shiftedX = similar(X_all)

        for j in eachindex(X_all)
            shiftedX[j] = X_all[j] + a_full[curve_id[j]]
        end

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

    # Gauge fixing
    function constrained_cost(a_free)

        a_full = vcat(0.0, a_free)

        return cost(a_full)

    end

    # Initial guess
    a0 = zeros(M-1)

    # Optimization
    res = optimize(constrained_cost, a0, NelderMead())

    shifts = vcat(0.0, Optim.minimizer(res))

    # === Compute collapse quality ===
    total_points = length(X_all)

    sX = sqrt(res.minimum / total_points)

    shiftedX = [X_all[i] + shifts[curve_id[i]] for i in eachindex(X_all)]

    sX_rel = sX / (maximum(shiftedX) - minimum(shiftedX) + eps())

    return shifts, X, Y, Yerr, sX_rel

end

# Function to perform parametric bootstrap for shift uncertainties
function shifts_parametric_mc(K_vals, t_vals, p2_mat, p2_err_mat; d=3, nbins=100, nmc=500, rng=MersenneTwister(0))
    M, N = size(p2_mat)
    all_shifts = zeros(nmc, M)

    for m in 1:nmc
        # Sample synthetic dataset
        noise = rand!(rng, Normal(), similar(p2_mat)) .* p2_err_mat
        p2_syn = p2_mat .+ noise

        # Ensure positivity (log will be used downstream)
        p2_syn = max.(p2_syn, eps())
        shifts_syn, _, _, _, _ = finite_time_scaling(K_vals, t_vals, p2_syn, p2_err_mat; d=d, nbins=nbins)
        all_shifts[m, :] .= shifts_syn
    end

    mean_shifts = vec(mean(all_shifts, dims=1))
    std_shifts  = vec(std(all_shifts, dims=1))
    return mean_shifts, std_shifts, all_shifts
end

# calculates rel var for arbitrary shifts
function tot_variance(a_full::Vector, X::Matrix, Y::Matrix; nbins=100)
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


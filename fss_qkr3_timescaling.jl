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

#Main function for finite time scaling analysis
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


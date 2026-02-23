include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("fss_qrk3_timescaling_analysis.jl") 

"""
Fit ξ(K) = ξ0 + A * |K - Kc|^{-ν} using LsqFit.jl

Returns best-fit parameters + standard errors from covariance matrix.
"""
function fit_xi_offset_LsqFit(K_vals, xi, xierr; n_k_filter=0, K_val_g=0, exclude_tol_frac=0.02)
        #filter points for fit    
        if n_k_filter != 0
                n_k_filter += 1  # account for Julia 1-based indexing    
                xi = xi[n_k_filter:end]
                xierr = xierr[n_k_filter:end]
                K_vals = K_vals[n_k_filter:end]
        end

    u = (1)./xi
    uerr = u.^2 .*xierr
    # Model function
    model(K, p) = abs(p[1]) .+ p[2] .* abs.(K .- p[4]).^(abs(p[3]))   #1/ξ(K) = β₀ + A|K−Kc|^{ν}
    names = ["β₀", "A", "ν", "Kc"]

    # Initial guess
    β₀₀ = minimum(u)*0.5
    A₀  = maximum(u)
    ν₀  = 1
    Kc₀ = K_vals[argmin(u) + K_val_g]  # where ξ is largest
    println(Kc₀)
    p0 = [β₀₀, A₀, ν₀, Kc₀]

    # Mask out values too close to trial Kc₀
    ΔK = maximum(K_vals) - minimum(K_vals)
    mask = abs.(K_vals .- Kc₀) .> exclude_tol_frac*ΔK
    Kfit, ufit, uerrfit = K_vals[mask], u[mask], uerr[mask]

    # Perform nonlinear least squares fit
    fit = curve_fit(model, Kfit, ufit, p0)
    pbest = coef(fit)
    covar = estimate_covar(fit)
    perr  = sqrt.(diag(covar))

    # goodness of fit
    residuals = ufit .- model(Kfit, pbest)
    if uerr[1] == 0.0
        χ2 = sum((residuals[2:end] ./ uerrfit[2:end]).^2)
    else
        χ2 = sum((residuals ./ uerrfit).^2)
    end
    dof = length(ufit) - length(pbest)
    χ2_red = χ2 / dof

    # Unpack results
    β0, A, ν, Kc = pbest
    err_β0, err_A, err_ν, err_Kc = perr


    return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
            err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


dim = 3 # spatial dimension
K_guess_index = 1 # index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit

perform_Kc_anal(shifts, shiftserr, K_vals; data_type="Ex", d=dim, n_k_filter=fit_k_filter, K_guess_index=K_guess_index)

#=#save data
d1 = DataFrame(xi', :auto)
CSV.write("ξ(k)_data_220.csv", d1)

d2 = DataFrame(ξfit', :auto)
CSV.write("ξ(k)_data_fit_220.csv", d2)

d3 = DataFrame(Kgrid', :auto)
CSV.write("Kvals_ξfit_220.csv", d3)

d4 = DataFrame(xierr', :auto)
CSV.write("ξerr(k)_data_220.csv", d4)=#
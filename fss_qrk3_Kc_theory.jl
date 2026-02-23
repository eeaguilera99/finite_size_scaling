include("imp_data_theory.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling


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
                dof = length(ufit[2:end]) - length(pbest)
        else
                χ2 = sum((residuals ./ uerrfit).^2)
                dof = length(ufit) - length(pbest)
        end
        χ2_red = χ2 / dof

        # Unpack results
        β0, A, ν, Kc = pbest
        err_β0, err_A, err_ν, err_Kc = perr


        return (β0=abs(β0), A=A, ν=abs(ν), Kc=Kc,
                err_β0=err_β0, err_A=err_A, err_ν=err_ν, err_Kc=err_Kc, χ2=χ2, χ2_red=χ2_red)
end


dim1 = 3
dim2 = 3   # spatial dimension


# Perform collapse
t_transient = 28
avg = true
amp = 2
K_guess_index = 0# index offset for initial Kc guess
fit_k_filter = 0 # number of low-K points to exclude from fit



perform_Kc_anal(shifts1, shifts1_err, K_vals; data_type="E_k", d=dim1, n_k_filter=fit_k_filter)
perform_Kc_anal(shifts2, shifts2_err, K_vals; data_type="1/n_c^2", d=dim2, n_k_filter=fit_k_filter)
        


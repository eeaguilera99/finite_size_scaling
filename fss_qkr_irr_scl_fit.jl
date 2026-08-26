include("imp_data_ex.jl")
include("fss_qrk3_timescaling_analysis.jl")

F01 = 1
F10 = 1
#F11 = 1
dim = 3
K_c_guess = 1.2
ν_guess = 1

tt_vals, _, _ = filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=5)

# Model function for fitting with corrections to scaling
function model(xy, p)
    K, t = xy[1,:], xy[2,:]

    b1  = p[1]
    #b2  = p[2]
    Kc  = p[2]
    ν   = p[3]
    F00 = p[4]

    ΔK = K .- Kc
    χK = b1 .* ΔK #.+ b2 .* ΔK.^2
    #return p[4] .+ (p[1].*((K .- p[2]))).*(t.^(1/(p[3]))).*F01
    return F00 .+
           F01 .* χK .* t.^(1.0 / (dim * ν)) #.+ p[5].*t.^(p[6]).*(F10 .+ p[1].*(p[2].- K).*(t.^(1/p[3]))).*p[7]
end

#Flatten
K_fit = vec([k for k in K_vals, t in tt_vals])
t_fit = vec([t for t in tt_vals, k in K_vals]')
Y_fit = vec(Y_data)

# Initial parameter guesses: b1, Kc, α, F00, ψ, y, F11
p0 = [1.0, K_c_guess, ν_guess, -20]#, -1, 1]  initial guesses
fit = curve_fit(model, [K_fit'; t_fit'], Y_fit, p0)
pbest = coef(fit)
perr  = sqrt.(diag(estimate_covar(fit)))
#b1, Kc, α, F00 = pbest
b1, Kc, ν, F00 = pbest #, ψ, y, F11




#quality of fit
Yerr_fit = vec(Yerr_data)
valid = isfinite.(Y_fit) .&
        isfinite.(Yerr_fit) .&
        (Yerr_fit .> 0)

residuals = Y_fit[valid] .-
            model(
                [K_fit[valid]'; t_fit[valid]'],
                pbest
            )

χ2 = sum((residuals ./ Yerr_fit[valid]).^2)

dof = sum(valid) - length(pbest)

χ2_red = χ2 / dof


#=
#plot raw data
plt1 = plot(title="Scaling funct vs K", xlabel = "κ", ylabel = "lnΛ")
for t in 1:Int(length(t_vals))
    plot!(plt1, K_vals, logΛ[:,t], seriestype=:scatter, label="data t=$(t_vals[t])", ms=3)
end
display(plt1)=#


Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
#=
#plot fit
plt2 = plot(title="Irelevant scaling fit d=$(dim)", xlabel = "κ", ylabel = "lnΛ")
for t in 1:Int(length(t_vals))
    logΛfit = model([Kgrid'; fill(t_vals[t], length(Kgrid))'], pbest)
    plot!(plt2, K_vals, logΛ[:,t], seriestype=:scatter, label="", ms=3)
    plot!(plt2, Kgrid, logΛfit, lw=2, label="")
end
display(plt2)

plt3 = plot(title="Scaling function with corrections d=$(dim)", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- Kgrid) .* (t_vals[t].^(1/ν))
    logΛfit = model([Kgrid'; fill(t_vals[t], length(Kgrid))'], pbest)
    plot!(plt3, -(1/dim).*log.(abs.(Xfit)), logΛfit, lw=2, label="")
end
display(plt3)=#

#=corrected scaling function
plt4 = plot(title="Corrected scaling function d=$(dim), \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ - correction)\$"))
correction = [(F10 + F11*(Kc - k)*t^(1/ν))*ψ*t^y for k in K_vals, t in t_vals]
logΛ_c = Y1 .- correction
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- K_vals) .* (t_vals[t].^(1/ν))
    plot!(plt4, -(1/dim).*log.(abs.(Xfit)), logΛ_c[:,t], seriestype=:scatter, ms=3, label="")
end
display(plt4)=#


#plot scaling function witouth Corrections

ΔK = K_vals .- Kc

χK = b1 .* ΔK #.+ b2 .* ΔK.^2

logxi = fill(NaN, length(K_vals))

for i in eachindex(χK)
    if abs(χK[i]) > 0
        logxi[i] = -ν * log(abs(χK[i]))
    end
end


plt5 = plot(title="Scaling by Taylor fit \$d=$(dim)\$, \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
for i in eachindex(K_vals)

    #fixed K
    K = K_vals[i]

    # Horizontal shift for each K
    X_collapse = X_data' .+ logxi[i] 

    # ACTUAL DATA
    Y_actual = Y_data[i, :]
    Y_err_actual = Yerr_data[i, :]

    #fit data
    # Use the fitted model to calculate ln Λ
    #t = tt_vals[j]

    #=logΛ_fit = model(
        [K_vals'; fill(t, length(K_vals))'],
        pbest
    )

    valid_fit = isfinite.(X_collapse) .&
            isfinite.(logΛ_fit)=#

    valid_data = isfinite.(X_collapse) .&
            isfinite.(Y_actual)

    plot!(
        plt5,
        X_collapse[valid_data],
        Y_actual[valid_data],#logΛ_fit[valid_fit],
        yerror=Y_err_actual[valid_data],
        seriestype = :scatter,
        ms = 3,
        label = ""
    )
end
display(plt5)

#plot localiation length ξ(k)
plt6 = plot(title="Localization length ξ(k) \$d=$(dim)\$, \$a_s=$(a_s)a_0\$", xlabel=latexstring("κ"), ylabel=latexstring("ξ(k)"))
plot!(plt6, K_vals, exp.(logxi), seriestype=:scatter, ms=3, label="ξ(k) data")
plot!(plt6, Kgrid, exp.(-ν*log.(abs.(b1.*(Kgrid .- Kc) .+ 0 .*(Kgrid .- Kc).^2))), lw=2, label="ξ(k) fit (ν=$(round(ν, digits=3)))")
vline!(plt6, [Kc], lw=2, ls=:dash, color=:red, label="Kc=$(round(Kc, digits=3))")
display(plt6)

#=plot irrelevant scaling fit results
plt6 = plot(title="Irrelevant scaling fit d=$(dim), \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"), xlims=(0,0.25))
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- K_vals) .* (t_vals[t].^(1/ν))
    logΛfit6 = model([K_vals'; fill(t_vals[t], length(K_vals))'], pbest) .- model([K_vals'; fill(t_vals[t], length(K_vals))'], [pbest[1], pbest[2], pbest[3], pbest[4], 0, 0, 0])
    plot!(plt6, -(1/dim).*log.(abs.(Xfit)), logΛfit6, lw=2, label="")
    println("F0+F1", model([K_vals'; fill(t_vals[t], length(K_vals))'], pbest) )
    println("F0", model([K_vals'; fill(t_vals[t], length(K_vals))'], [pbest[1], pbest[2], pbest[3], pbest[4], 0, 0, 0]))
end

display(plt6)=#



println("Collapse \$d=$(D)\$, \$a_s=$(a_s)a_0\$ fitted parameters with errors:")
println("b1 = $(round(b1, digits=3)) ± $(round(perr[1], digits=3))")
println("κ_c = $(round(Kc, digits=3)) ± $(round(perr[2], digits=3))")
println("ν = $(round(ν, digits=3)) ± $(round(perr[4], digits=3))")
println("χ2 = $(round(χ2, digits=3)), χ2_red=$(round(χ2_red, digits=3))")




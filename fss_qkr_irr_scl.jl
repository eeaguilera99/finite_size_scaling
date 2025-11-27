include("fss_qkr3_timescaling.jl")

#Irrelevant scaling parameters
dim = 3
α = 1.5087 
α_err = 0.0010071
y = -8.143
y_err = 0.0
F00 = -22.784 
F00_err = 4.5488e-3
F01 = 1
F10 = 1
F11 = -1.2289e1 
F11_err = 0.0
b0 = -4.5853e-2
b0_err = 1.3437e-3
b1 = 4.5862e-2 
b1_err = 1.3888e-3
ψ = -2.1162e1
ψ_err = 0.0
Kc = 1.0801
Kc_err = 1.6544e-3


t_vals = t_vals[5:end]
p2_mat = p2_mat[:,5:end]

function XX(K,t)
    return (b0*(Kc-K)/Kc)*t^(1/(α))#+ b1*((Kc-K)/Kc)^2
end

function scaling_funct(K, t)
    return F00 + F01*XX(K,t) + (F10 + F11*XX(K,t))*ψ*t^y
end 
function corrected_scaling_funct(K, T, M)
    F1 = [(F10 + F11*XX(k,t))*ψ*t^y for k in K, t in T]
    return M .- F1
end

function inverse_ξ(K)
    return abs(b0*(Kc-K)/Kc)^(α/dim)
end

F = [scaling_funct(k, t) for k in K_vals, t in t_vals]
  #matrix of scaled data

X = [XX(k, t) for k in K_vals, t in t_vals]  #matrix of scaled x-coordinates


res, shifts, X1, Y1, Y1err, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)


#=
plt1 = plot(title="fss_qkr_irr_scaling d=$(dim)",
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)",
    size=(800,600))
scatter!(-(α/3).*log.(abs.(X)), F, mc=:blue)
display(plt1)=#

Λ = p2_mat ./ (t_vals' .^ (2/dim))
ln_Λ_c = corrected_scaling_funct(K_vals, t_vals, log.(Λ))

#=
#plot irr scaling fit results
plt2 = plot(title=latexstring("Finite-time scaling with correction \$d=$(dim)\$"),
    xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
scatter!(vec(-(α/dim).*log.(abs.(X))), vec(ln_Λ_c), mc=:red, label="Corrected")
#scatter!(-(α/3).*log.(abs.(X)), F, mc=:blue, label="")
scatter!(vec(X1 .+ shifts), vec(Y1), mc=:blue, label="Uncorrected")
#xlims!(-1, 7)
display(plt2)=#

#plot correlation length from fit
plt3 = plot(title=latexstring("Correlation length from fit \$d=$(dim)\$"),
    xlabel=latexstring("\$κ\$"), ylabel=latexstring("\$ξ(κ)\$"))
ξ_fit = [1/inverse_ξ(k) for k in K_vals]
scatter!(K_vals, ξ_fit, mc=:green, label="ξ from fit")
K_grid = range(minimum(K_vals), maximum(K_vals), length=400)
plot!(plt3, K_grid, [1/inverse_ξ(k) for k in K_grid], lw=2, label="")
vline!(plt3, [Kc], linestyle=:dash, color=:red, label=L"\kappa_c="*"$(round(Kc,digits=3))")
display(plt3)

# One-parameter scaling collapse using ξ from fit
fit_shifts = log.( [1/inverse_ξ(k) for k in K_vals] )
plt4 = plot(title=latexstring("Collapse using \$ξ\$ from fit \$d=$(dim)\$"),
    xlabel=latexstring("ln \$(ξ/N^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
for (i,K) in enumerate(K_vals)
    plot!(plt4, X1[:] .+ fit_shifts[i], Y1[i,:], yerror=Y1err[i,:], marker=:o, label="")
end
display(plt4)

#=
# Comparison plot
plt5 = plot(title="One-parameter scaling collapse vs Corrections",
    xlabel=latexstring("ln \$(ξ/N^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
scatter!(plt5, vec(X1 .+ [fit_shifts[i] for i in 1:length(K_vals), t in 1:size(X1,2)]), vec(Y1), mc=:red, marker=:o, label="corrections")
scatter!(plt5, vec(X1 .+ [shifts[i] for i in 1:length(K_vals), t in 1:size(X1,2)]).+ 2.05, vec(Y1), mc=:blue, marker=:o, label="1-param scaling")
display(plt5)=#
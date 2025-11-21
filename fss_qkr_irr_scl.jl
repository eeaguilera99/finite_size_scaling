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
    return (b0*(K-Kc)/Kc + b1*((K-Kc)/Kc)^2)*t^(1/α)
end

function scaling_funct(K, t)
    return F00 + F01*XX(K,t) + (F10 + F11*XX(K,t))*ψ*t^y
end 
function corrected_scaling_funct(K, T, M)
    F1 = [(F10 + F11*XX(k,t))*ψ*t^y for k in K, t in T]
    return M .- F1
end

F = [scaling_funct(k, t) for k in K_vals, t in t_vals]
  #matrix of scaled data

X = [(k-Kc)/Kc*t^(1/α) for k in K_vals, t in t_vals]  #matrix of scaled x-coordinates


res, shifts, X1, Y1, Y1err, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)


#=
plt1 = plot(title="fss_qkr_irr_scaling d=$(dim)",
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)",
    size=(800,600))
scatter!(-(α/3).*log.(abs.(X)), F, mc=:blue)
display(plt1)=#

Λ = p2_mat ./ (t_vals' .^ (2/dim))
ln_Λ_c = corrected_scaling_funct(K_vals, t_vals, log.(Λ))

plt2 = plot(title=latexstring("Finite-time scaling with correction \$d=$(dim)\$"),
    xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
scatter!(vec(-(α/3).*log.(abs.(X))), vec(ln_Λ_c), mc=:red, label="Corrected")
#scatter!(-(α/3).*log.(abs.(X)), F, mc=:blue, label="")
scatter!(vec(X1 .+ shifts.+0.5), vec(Y1), mc=:green, label="Uncorrected")
display(plt2)   
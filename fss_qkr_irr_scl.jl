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

function XX(K,t)
    return b0*(K-Kc)*t^(1/α)+b1*((K-Kc)*t^(1/α))^2
end

function scaling_funct(K, t)
    return F00 + F01*XX(K,t) + (F10 + F11*XX(K,t))*ψ*t^y
end 

F = -log.(abs.([scaling_funct(k, t) for k in K_vals, t in t_vals]))  #matrix of scaled data

X = -(1/(α/3))*log.(abs.([XX(k, t) for k in K_vals, t in t_vals]))  #matrix of scaled x-coordinates

plot(title="fss_qkr_irr_scaling d=$(dim)",
    xlabel="ln(ξ/N^(1/d))", ylabel="ln(Λ)",
    size=(800,600))
scatter!(F, X)

include("fss_qkr3_timescaling.jl")

#Irrelevant scaling parameters
dim = 6

# Replace with the path to your actual CSV file
filepath = "irr_scale_data\\param_mR=1.csv"

# Define aliases
rename_map = Dict("nu" => "α", "c0" => "ψ")

# Read the file manually since it's not a standard CSV
lines = readlines(filepath)

# Iterate and define variables dynamically
for line in lines
    if occursin(":", line)
        name_val = split(line, ":")
        varname = strip(name_val[1])
        valstr = strip(name_val[2])
        
        # Rename if applicable
        varname = get(rename_map, varname, varname)

        if occursin("±", valstr)
            val, err = strip.(split(valstr, "±"))
            value = parse(Float64, val)
            error = parse(Float64, err)
        else
            value = parse(Float64, valstr)
            error = 0.0
        end

        # Define variables dynamically
        @eval begin
            $(Symbol(varname)) = $value
            $(Symbol(varname * "_err")) = $error
        end
    end
end

#=
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
=#
ν = α/dim

t_vals = t_vals[5:end]
p2_mat = p2_mat[:,5:end]

function χ(K)
    return (b0*-(Kc-K) + 0*(-(Kc-K))^2 + 0*(-(Kc-K))^2)# + b1*((Kc-K)/Kc)^2
end

function scaling_funct(K, t)
    return F00 + F01*χ(K)*t^(1/(α)) + (F10 + F11*χ(K)*t^(1/(α)))*ψ*t^y
end 
function corrected_scaling_funct(K, T, M)
    F1 = [(F10 + F11*χ(k)*t^(1/(α)))*ψ*t^y for k in K, t in T]
    return M .- F1
end

function inverse_ξ(K)
    return abs(χ(K))^(ν)
end


# One-parameter scaling collapse
res, shifts, X1, Y1, Y1err, s_rel = finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; d=dim, n_kicks_i=1, n_kicks_f=0)

# Compute scaled data
F = [scaling_funct(k, t) for k in K_vals, t in t_vals]
X = [χ(k)*t^(1/α) for k in K_vals, t in t_vals]  #matrix of scaled x-coordinates
ln_Λ_c = corrected_scaling_funct(K_vals, t_vals, Y1) #matrix of corrected lnΛ values


#plot irr scaling fit results
plt1 = plot(title=latexstring("Finite-time scaling with correction \$d=$(dim)\$"),
    xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
scatter!(vec(-(ν).*log.(abs.(X))), vec(ln_Λ_c), mc=:red, label="Corrected")
#scatter!(-(α/3).*log.(abs.(X)), F, mc=:blue, label="")
scatter!(vec(X1 .+ shifts) .+ 1.21 , vec(Y1), mc=:blue, label="Uncorrected")
#xlims!(-1, 7)
display(plt1)

#=
#plot correlation length from fit
plt2 = plot(title=latexstring("Correlation length from fit \$d=$(dim)\$"),
    xlabel=latexstring("\$κ\$"), ylabel=latexstring("\$ξ(κ)\$"))
ξ_fit = [1/inverse_ξ(k) for k in K_vals]
scatter!(K_vals, ξ_fit, mc=:green, label="ξ from fit")
K_grid = range(minimum(K_vals), maximum(K_vals), length=400)
plot!(plt2, K_grid, [1/inverse_ξ(k) for k in K_grid], lw=2, label="")
vline!(plt2, [Kc], linestyle=:dash, color=:red, label=L"\kappa_c="*"$(round(Kc,digits=3))")
display(plt2)


# One-parameter scaling collapse using ξ from fit
fit_shifts = log.( [1/inverse_ξ(k) for k in K_vals] )
plt3 = plot(title=latexstring("Collapse using \$ξ\$ from fit \$d=$(dim)\$"),
    xlabel=latexstring("ln \$(ξ/N^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
for (i,K) in enumerate(K_vals)
    plot!(plt3, X1[:] .+ fit_shifts[i], Y1[i,:], yerror=Y1err[i,:], marker=:o, label="")
end
display(plt3)


# Comparison plot
plt4 = plot(title="One-parameter scaling collapse vs Corrections",
    xlabel=latexstring("ln \$(ξ/N^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
scatter!(plt4, vec(X1 .+ [fit_shifts[i] for i in 1:length(K_vals), t in 1:size(X1,2)]), vec(Y1), mc=:red, marker=:o, 
label=latexstring("\$ξ\$ from corrections"))
scatter!(plt4, vec(X1 .+ [shifts[i] for i in 1:length(K_vals), t in 1:size(X1,2)]).+ 2.3, vec(Y1), mc=:blue, marker=:o, label="1-param scaling")
display(plt4)=#
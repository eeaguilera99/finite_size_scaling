include("fss_qkr3_timescaling.jl")
include("fss_qrk3_timescaling_analysis.jl")

F01 = 1
F10 = 1
#F11 = 1
dim = 3
K_c_guess = 1.2

t_vals, _, _ = filter_Nkicks(t_vals, p2_mat, p2_err_mat; n_kicks_i=5, n_kicks_f=0)

# Model function for fitting with corrections to scaling
function model(xy, p)
    K, t = xy[1,:], xy[2,:]
    return p[4] .+ p[1].*((p[2] .- K)).*(t.^(1/p[3])).*F01 #.+ p[5].*t.^(p[6]).*(F10 .+ p[1].*(p[2].- K).*(t.^(1/p[3]))).*p[7]
end

#Flatten
X, Y = vec([k for k in K_vals, t in t_vals]), vec([t for t in t_vals, K in K_vals])

# Initial parameter guesses: b1, Kc, α, F00, ψ, y, F11
p0 = [1, K_c_guess, dim*0.5, -24.0]#, 2, -1, 1]  initial guesses
fit = curve_fit(model, [X'; Y'], vec(Y_data), p0)
pbest = coef(fit)
b1, Kc, α, F00 = pbest #, ψ, y, F11


#quality of fit
residuals = vec(Y_data) .- model([X'; Y'], pbest)
χ2 = sum((residuals ./ vec(Yerr_data)).^2)
dof = length(vec(Y_data)) - length(pbest)


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
    Xfit = (Kc .- Kgrid) .* (t_vals[t].^(1/α))
    logΛfit = model([Kgrid'; fill(t_vals[t], length(Kgrid))'], pbest)
    plot!(plt3, -(1/dim).*log.(abs.(Xfit)), logΛfit, lw=2, label="")
end
display(plt3)=#

#=corrected scaling function
plt4 = plot(title="Corrected scaling function d=$(dim), \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ - correction)\$"))
correction = [(F10 + F11*(Kc - k)*t^(1/α))*ψ*t^y for k in K_vals, t in t_vals]
logΛ_c = Y1 .- correction
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- K_vals) .* (t_vals[t].^(1/α))
    plot!(plt4, -(1/dim).*log.(abs.(Xfit)), logΛ_c[:,t], seriestype=:scatter, ms=3, label="")
end
display(plt4)=#


#plot scaling function witouth Corrections
plt5 = plot(title="Scaling function without corrections d=$(dim), \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"))
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- K_vals) .* (t_vals[t].^(1/α))
    logΛ_nc = model([K_vals'; fill(t_vals[t], length(K_vals))'], [pbest[1], pbest[2], pbest[3], pbest[4]])
    plot!(plt5, -(1/dim).*log.(abs.(Xfit)), logΛ_nc, seriestype=:scatter, ms=3, label="")
end
display(plt5)

#=plot irrelevant scaling fit results
plt6 = plot(title="Irrelevant scaling fit d=$(dim), \$a_s=$(a_s)a_0\$", xlabel=latexstring("ln \$(ξ/t^{1/d})\$"), ylabel=latexstring("ln \$(Λ)\$"), xlims=(0,0.25))
for t in 1:Int(length(t_vals))
    Xfit = (Kc .- K_vals) .* (t_vals[t].^(1/α))
    logΛfit6 = model([K_vals'; fill(t_vals[t], length(K_vals))'], pbest) .- model([K_vals'; fill(t_vals[t], length(K_vals))'], [pbest[1], pbest[2], pbest[3], pbest[4], 0, 0, 0])
    plot!(plt6, -(1/dim).*log.(abs.(Xfit)), logΛfit6, lw=2, label="")
    println("F0+F1", model([K_vals'; fill(t_vals[t], length(K_vals))'], pbest) )
    println("F0", model([K_vals'; fill(t_vals[t], length(K_vals))'], [pbest[1], pbest[2], pbest[3], pbest[4], 0, 0, 0]))
end

display(plt6)=#

println("Fitted parameters with errors:")
println("b1=$(b1)")
println("Kc=$(Kc)")
println("α=$(α)")
println("χ2=$(χ2), χ2_red=$(χ2/dof)")




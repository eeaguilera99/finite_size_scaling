include("fss_qkr3_timescaling.jl")

F01 = 1
F10 = 1
dim = 3

t_vals = t_vals[5:end]
p2_mat = p2_mat[:,5:end]

# Model function for fitting with corrections to scaling
function model(xy, p)
    K, t = xy[1,:], xy[2,:]
    return p[4] .+ p[1].*((1 .- K./p[2])).*(t.^(1/p[3])).*F01 .+ p[5].*t.^(p[6]).*(F10 .+ p[1].*(1 .- K./p[2]).*(t.^(1/p[3])).*p[7])
end

#Flatten
X , Y = vec([k for k in K_vals, t in t_vals]), vec([t for t in t_vals, K in K_vals])

Λ = p2_mat ./ (t_vals' .^ (2/dim))
logΛ = log.(Λ)

b1_0 = 
p0 = [1, 1.2, 1.5, 1, 1, -1, -1]  # initial guesses
fit = curve_fit(model, [X'; Y'], vec(logΛ), p0)
pbest = coef(fit)
b1, Kc, α, F00, ψ, y, F11 = pbest

#plot raw data
plt1 = plot(title="Scaling funct vs K", xlabel = "κ", ylabel = "lnΛ")
for t in 1:length(t_vals)
    plot!(plt1, K_vals, logΛ[:,t], seriestype=:scatter, label="data t=$(t_vals[t])", ms=3)
end
display(plt1)

#plot fit
plt2 = plot(title="Irelevant scaling fit d=$(dim)", xlabel = "κ", ylabel = "lnΛ")
for t in 1:length(t_vals)
    Kgrid = range(minimum(K_vals), maximum(K_vals), length=400)
    logΛfit = model([Kgrid'; fill(t_vals[t], length(Kgrid))'], pbest)
    plot!(plt2, K_vals, logΛ[:,t], seriestype=:scatter, label="", ms=3)
    plot!(plt2, Kgrid, logΛfit, lw=2, label="")
end
display(plt2)

#scatter!(vec(-(α/dim).*log.(abs.(X))), vec(logΛ), mc=:blue, label="Data")
#plot!(vec(-(α/dim).*log.(abs.(X))), vec(model([X'; Y'], pbest)), mc=:blue, label="Fit")


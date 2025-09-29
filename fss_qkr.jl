
using CSV
using DataFrames
using FiniteSizeScaling
using Plots

# Read number of kicks (column vector)
number_of_kicks = CSV.read("Codes\\finite_size_scaling\\data\\Number_of_kicks.csv", DataFrame; header=false)[!, 1]

# Read kick strength (column vector)
kick_strength = CSV.read("Codes\\finite_size_scaling\\data\\kappa.csv", DataFrame; header=false)[!, 1]

# Read energy matrix (rows: kick strength, columns: number of kicks)
energy = Matrix(CSV.read("Codes\\finite_size_scaling\\data\\nc_matrix.csv", DataFrame; header=false))

# Read error matrix (same shape as energy)
error = Matrix(CSV.read("Codes\\finite_size_scaling\\data\\nc_err_matrix.csv", DataFrame; header=false))

data_with_err = [] #Array of tuples (x, y, error, param), x is kick strength, y is energy, param is number of kicks
fit_weights = []
global_kicks = size(kick_strength,1)-2 #custom number of kicks strength to consider
global_nkicks = size(energy,2) #custom number of kicks to consider
for col in 5:global_nkicks
    kicks = Float64.(kick_strength[1:global_kicks])
    energies = Float64.(energy[1:global_kicks, col])
    errors = Float64.(error[1:global_kicks, col])
    n_kicks = Float64(number_of_kicks[col])
    push!(data_with_err, (kicks, energies, errors, n_kicks))
    #push!(fit_weights, 1.0 ./ (error[1:global_kicks, col] .^ 2))
    push!(fit_weights, ones(size(error[1:global_kicks, col])))
end

# Scaled functions
function x_scaled(X, t, v1)
    return ((X.-v1).*(t^(-1/3)))
end
function y_scaled(Y, t, v1)
    return (Y*(t^(-2/3)))
end

#call fss
scaled_data, residuals, min_res, best_v1 = fss_one_var(data=data_with_err, 
xs=x_scaled, ys=y_scaled, v1i=0.5, v1f=2, n1=500, p=6, weights=fit_weights)


# Custom function to plot scaled_data using Plots

function plot_scaled_data(scaled_data; xlabel="Scaled x", ylabel="Scaled y", 
    title="Finite Size Scaling Plot κ_c=$(best_v1)")
    plt = plot()
    for d in scaled_data
        x = d[1]
        y = d[2]
        err = length(d) >= 3 ? d[3] : nothing
        param = length(d) >= 4 ? d[4] : "param"
        if err !== nothing
            plot!(plt, x, y, yerror=err, seriestype=:scatter, label="Nkicks = $(param)")
        else
            plot!(plt, x, y, seriestype=:scatter, label="Nkicks = $(param)")
        end
    end
    xlabel!(plt, xlabel)
    ylabel!(plt, ylabel)
    title!(plt, title)
    display(plt)
end

# Custom function to plot data with log-log axes using Plots
function plot_data_log(data)
    plt = plot()
    for d in data
        x = d[1]
        y = d[2]
        err = length(d) >= 3 ? d[3] : nothing
        param = length(d) >= 4 ? d[4] : "param"
        if err !== nothing
            plot!(plt, x, y, yerror=err, seriestype=:scatter, label="Nkicks = $(param)")
        else
            plot!(plt, x, y, seriestype=:scatter, label="Nkicks = $(param)")
        end
    end
    xlabel!(plt, xlabel)
    ylabel!(plt, ylabel)
    if v1_opt !== nothing
        title!(plt, "$(title) (optimal v₁ = $(v1_opt))")
    else
        title!(plt, title)
    end
    display(plt)
end

# Example usage:
#plot_data(data_with_err)
#plot_scaled_data(scaled_data, xlabel="(κ - κ_c)  t^(-1/3)", ylabel="⟨p²⟩ * t^(-2/3)")
plot_data_log(data_with_err, xlabel=" t^(-1/3)", ylabel="⟨p²⟩ * t^(-2/3)")

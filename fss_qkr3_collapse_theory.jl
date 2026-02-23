include("imp_data_theory.jl")
include("fss_qkr3_timescaling.jl")

"""
    finite_time_scaling(K_vals, t_vals, p2_mat, p2_err_mat; nbins=30)

Perform finite-time scaling collapse of Anderson transition data.

# Arguments
- `K_vals::Vector`: Kick strengths (length M).
- `t_vals::Vector`: Times (length N).
- `p2_mat::Matrix`: M×N matrix of ⟨p²⟩ values.
- `p2_err_mat::Matrix`: M×N matrix of errors (same shape).
- `nbins::Int`: Number of bins in Y direction (default = 30).

# Returns
- `shifts::Vector`: Optimal horizontal shifts aᵢ = ln ξ(Kᵢ).
- `(X, Y)`: Arrays of logarithmic coordinates.
"""

perform_collapse(K_vals, X1, Y1, Yerr1, shifts1, s_rel1; data_type=" Th E_k", raw=true, d=dim1)
perform_collapse(K_vals, X2, Y2, Yerr2, shifts2, s_rel2; data_type=" Th nc^-2", raw=true, d=dim2)



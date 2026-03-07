#include("imp_data_ex.jl")
include("fss_qkr3_timescaling.jl")  # for finite_time_scaling
include("fss_qrk3_timescaling_analysis.jl")

"""
finite_time_linear_scaling(K_vals, t_vals, p2_mat; Kc, ΔKfit=0.3, plotshow=true)

Implements the Lemarié finite-time-scaling method:

1. Compute Λ = <p²>/t^(2/3).
2. For each t, fit ln Λ ≈ ln Λc + s(t)*(K−Kc)
   using only points |K−Kc| < ΔKfit.
3. Plot ln Λ(K) vs K (figure 13-like).
4. Plot ln|s(t)| vs ln t and extract ν from the slope (figure 14-like).

Returns a NamedTuple with ν, slope, intercept, and the vectors of s(t).
"""


Kc = round(Kc_1, digits=3) # use Kc from analysis
dim = 5
t_transient = 2
finite_time_linear_scaling(K_vals, t_vals, p2_mat, p2_err_mat, dim;
                                 Kc = Kc, ΔKfit = 0.5, n_kicks_i=t_transient, n_kicks_f=0)


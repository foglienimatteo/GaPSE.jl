

data = readdlm(FILE_BACKGROUND, comments = true)
data_dict = Dict([name => reverse(data[:, i]) for (i, name) in enumerate(NAMES_BACKGROUND)]...)

D = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. D"])
f = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"])
ℋ = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"] .* data_dict["gr.fac. D"])
ℋ_p(s) = derivative(spl, s)
s_b = Spline1D(data_dict["comov. dist."] .* h_0, [2.0 / 5.0 for i in 1:length(data_dict["comov. dist."])])
f_evo = 0

ps = readdlm(FILE_PS, comments = true)
ps_dict = Dict([name => ps[:, i] for (i, name) in enumerate(NAMES_PS)]...)

PK = Spline1D(ps_dict["k (h/Mpc)"], ps_dict["P (Mpc/h)^3"])

N = 1024                                # number of points to use in the Fourier transform
kmax = ps_dict["k (h/Mpc)"][end]        # maximum k-value
kmin = ps_dict["k (h/Mpc)"][begin]      # minimum k-value
s0 = 1 / kmax;                            # minimum r-value (should be ~1/kmax)

I00 = Spline1D(xicalc(PK, 0, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I20 = Spline1D(xicalc(PK, 2, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I40 = Spline1D(xicalc(PK, 4, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I02 = Spline1D(xicalc(PK, 0, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I22 = Spline1D(xicalc(PK, 2, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I31 = Spline1D(xicalc(PK, 3, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I13 = Spline1D(xicalc(PK, 1, 3; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I11 = Spline1D(xicalc(PK, 1, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)

σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), kmin, kmax)
σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), kmin, kmax)
σ_2 = quadgk(q -> PK(q) / (2 * π^2), kmin, kmax)
σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), kmin, kmax)

ℛ(s) = 5 * s_b(s) + (2 - 5 * s_b(s)) / (ℋ(s) * s) + ℋ_p(s)/(ℋ(s)^2) - f_evo


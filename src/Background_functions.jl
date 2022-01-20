

data = readdlm(FILE_BACKGROUND, comments = true)
data_dict = Dict([name => reverse(data[:, i]) for (i, name) in enumerate(NAMES_BACKGROUND)]...)

D = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. D"])
f = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"])
H_com = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"] .* data_dict["gr.fac. D"])


ps = readdlm(FILE_PS, comments = true)
ps_dict = Dict([name => ps[:, i] for (i, name) in enumerate(NAMES_PS)]...)

PK = Spline1D(ps_dict["k (h/Mpc)"], ps_dict["P (Mpc/h)^3"])

N = 1024     # number of points to use in the Fourier transform
kmax = 1e3   # maximum k-value
kmin = 1e-5  # minimum k-value
s0 = 1e-3;   # minimum r-value (should be ~1/kmax)

I00 = Spline1D(xicalc(PK, 0, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I20 = Spline1D(xicalc(PK, 2, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I40 = Spline1D(xicalc(PK, 4, 0; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I02 = Spline1D(xicalc(PK, 0, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I22 = Spline1D(xicalc(PK, 2, 2; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I31 = Spline1D(xicalc(PK, 3, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I13 = Spline1D(xicalc(PK, 1, 3; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)
I11 = Spline1D(xicalc(PK, 1, 1; N = N, kmin = kmin, kmax = kmax, r0 = s0)...)





F_map_data = readdlm(FILE_F_MAP, comments = true)
F_map_data_dict = Dict([name => F_map_data[2:end, i] for (i, name) in enumerate(NAMES_F_MAP)]...)

F_grid = RectangleGrid(unique(F_map_data_dict["mu"]), unique(F_map_data_dict["x"]))
spline_F(x, μ) = GridInterpolations.interpolate(F_grid, F_map_data_dict["F"], [μ, x])


data = readdlm(FILE_BACKGROUND, comments = true)
N_Z_MIN_10 = findfirst(z -> z < 10, data[:, 1])
println("N_Z_MIN_10 = ", N_Z_MIN_10)
data_dict = Dict([name => reverse(data[:, i][N_Z_MIN_10:end]) for (i, name) in enumerate(NAMES_BACKGROUND)]...)

D = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. D"])
f = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"])
ℋ = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["gr.fac. f"] .* data_dict["gr.fac. D"])
ℋ_p(s) = Dierckx.derivative(ℋ, s)
s_b = Spline1D(data_dict["comov. dist."] .* h_0, [2.0 / 5.0 for i in 1:length(data_dict["comov. dist."])])
s_of_z = Spline1D(data_dict["z"], data_dict["comov. dist."] .* h_0)
z_of_s = Spline1D(data_dict["comov. dist."] .* h_0, data_dict["z"])
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

σ_0 = quadgk(q -> PK(q) * q^2 / (2 * π^2), kmin, kmax)[1]
σ_1 = quadgk(q -> PK(q) * q / (2 * π^2), kmin, kmax)[1]
σ_2 = quadgk(q -> PK(q) / (2 * π^2), kmin, kmax)[1]
σ_3 = quadgk(q -> PK(q) / (2 * π^2 * q), kmin, kmax)[1]

ℛ(s) = 5 * s_b(s) + (2 - 5 * s_b(s)) / (ℋ(s) * s) + ℋ_p(s) / (ℋ(s)^2) - f_evo

s2(s1, s, μ) = √(s1^2 + s^2 + 2 * s1 * s * μ)
y(s1, s, μ) = (μ * s + s1) / s2(s, s1, μ)


s_min = s_of_z(z_MIN)
s_max = s_of_z(z_MAX)
ϕ(s; s_min = s_min, s_max = s_max) = s_min < s < s_max ? 1.0 : 0.0
W(θ; θ_max = θ_MAX) = 0 < θ < θ_max ? 1.0 : 0.0
function V(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     r1, r2 = s_min * sin(θ_max), s_max * sin(θ_max)
     d1, d2 = s_min * cos(θ_max), s_max * cos(θ_max)
     calotta_up = π / 3 * (s_max - d2)^2 * (2 * s_max + d2)
     calotta_down = π / 3 * (s_min - d1)^2 * (2 * s_min + d1)
     tronco_cono = π / 3 * (r1^2 + r1 * r2 + r2^2) * (s_max - s_min) * cos(θ_max)
     return calotta_up + tronco_cono - calotta_down
end

function z_eff(s_min = s_min, s_max = s_max, θ_max = θ_MAX)
     int_w2 = 2 * π * (1 - cos(θ_max))
     int_z_ϕ2 = quadgk(s -> s^2 * z_of_s(s), s_min, s_max)[1]
     int_ϕ2 = quadgk(s -> s^2, s_min, s_max)[1]

     return int_z_ϕ2 / int_ϕ2
end
s_eff = s_of_z(z_eff())


println("z_eff = ", z_eff())
println("s_eff = ", s_eff)

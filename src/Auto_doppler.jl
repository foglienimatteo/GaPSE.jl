

function J00(s1, s2, y)
     1 / 45 * f(s1) * f(s2) * ℋ(s1) * ℋ(s2) * ℛ(s1) * ℛ(s2) *
     (y^2 * s1 * s2 - 2 * y * (s1^2 + s2^2) + 3 * s1 * s2)
end

function J02(s1, s2, y)
     2 / 63 * f(s1) * f(s2) * ℋ(s1) * ℋ(s2) * ℛ(s1) * ℛ(s2) *
     (y^2 * s1 * s2 - 2 * y * (s1^2 + s2^2) + 3 * s1 * s2)
end

function J04(s1, s2, y)
     1 / 105 * f(s1) * f(s2) * ℋ(s1) * ℋ(s2) * ℛ(s1) * ℛ(s2) *
     (y^2 * s1 * s2 - 2 * y * (s1^2 + s2^2) + 3 * s1 * s2)
end

function J20(s1, s2, y)
     s = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
     1 / 3 * y * s^2 * f(s1) * f(s2) * ℋ(s1) * ℋ(s2) * ℛ(s1) * ℛ(s2)
end

function J31(s1, s2, y)
     - y * f(0) * ℋ(0) * s1^2 * f(s1) * ℛ(s1) * (ℛ(s2) - 5*s_b(s2) + 2) 
end

function J11(s1, s2, y)
    1 / 5 * y * f(0) * ℋ(0) * s1^2 * f(s1) * ℋ(s1) * ℛ(s1) * (ℛ(s2) - 5 * s_b(s2) + 2)
end

function J13(s1, s2, y)
     1 / 5 * y * f(0) * ℋ(0) * s1^2 * f(s1) * ℋ(s1) * ℛ(s1) * (ℛ(s2) - 5 * s_b(s2) + 2)
end

function Jσ2(s1, s2, y)
     1 / 3 * y * f(0)^2 * ℋ(0)^2 * (ℛ(s1) - 5 * s_b(s1) + 2) * (ℛ(s2) - 5 * s_b(s2) + 2)
end


##########################################################################################92


function ξ_doppler(s1, s2, y)

     #=
     D(s1) * D(s2) * ( J00(s1, s2, y) * I00(s) + J02(s1, s2, y) * I20(s) +
          J04(s1, s2, y) * I40(s) + J20(s1, s2, y) * I02(s) ) + 
     D(s1) * ( J31(s1, s2, y) * I13(s1) + J11(s1, s2, y) * I11(s1) + 
          J13(s1, s2, y) * I31(s1) ) +
     D(s2) * (J31(s2, s1, y) * I13(s2) + J11(s2, s1, y) * I11(s2) + 
          J13(s2, s1, y) * I31(s2)) +
     Jσ2(s1, s2, y) * σ_2
     =#

     #=
     D1, D2 = D(s1), D(s2)
     f0, f1, f2 = f(0), f(s1), f(s2)
     ℋ0, ℋ1, ℋ2 = ℋ(0), ℋ(s1), ℋ(s2)
     ℛ0, ℛ1, ℛ2 = ℛ(0), ℛ(s1), ℛ(s2)
     s_b1, s_b2 = s_b(s1), s_b(s2)

     s = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
     prefac = f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     J00 = 1 / 45 * prefac * c1
     J02 = 2 / 63 * prefac * c1
     J04 = 1 / 105 * prefac * c1
     J20 = 1 / 3 * y * s^2 * c1
     =#

     #=
     prefac_12 = f0 * ℋ0 * s1^2 * f1 * ℛ1
     prefac_21 = f0 * ℋ0 * s2^2 * f2 * ℛ2
     parenth_1 = (ℛ1 - 5 * s_b1 + 2)
     parenth_2 = (ℛ2 - 5 * s_b2 + 2)

     J31_12 = -y * prefac_12 * parenth_2
     J11_12 = 1 / 5 * y * prefac_12 * ℋ1 * parenth_2
     J13_12 = 1 / 5 * y * prefac_12 * ℋ1 *parenth_2

     J31_21 = -y * prefac_21 * parenth_1
     J11_21 = 1 / 5 * y * prefac_21 * ℋ2 * parenth_1
     J13_21 = 1 / 5 * y * prefac_21 * ℋ2 *parenth_1

     Jσ2 = 1 / 3 * y * f0 * ℋ0^2 * parenth_1 * parenth_2
     =#

     #return D1 * D2 * (J00 * I00(s) + J02 * I20(s) + J04 * I40(s) + J20 * I02(s))


     s = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2
     c2 = (1.0 / 3.0) * y * s^2

     D1 = D(s1)
     D2 = D(s2)
     f1 = f(s1)
     f2 = f(s2)
     H1 = ℋ(s1)
     H2 = ℋ(s2)
     R1 = 1 - 1.0 / (H1 * s1)
     R2 = 1 - 1.0 / (H2 * s2)
     prefac = D1 * D2 * f1 * f2 * R1 * R2 * H1 * H2

     parenth = (1.0 / 45.0) * I00(s) + (2.0 / 63.0) * I20(s) + (1.0 / 105.0) * I40(s)

     return prefac * (c1 * parenth + c2 * I02(s))
end


function int_on_mu_doppler(s1, s, μ)
     if ϕ(s2(s1, s, μ)) > 0
          return ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ)) * spline_F(s / s1, μ)
     else
          return 0.0
     end
end


function integral_on_mu_doppler(s1, s; kwargs...)
     return quadgk(μ -> int_on_mu_doppler(s1, s, μ), -1, 1; kwargs...)[1]
end

function map_integral_on_mu_doppler(s1 = s_eff; kwargs...)
     ss = 10 .^ range(-1, 3, length = 100)
     xis = [integral_on_mu_doppler(s1, s; kwargs...) for s in ss]
     return (ss, xis)
end


function PS_doppler(L::Integer = 0; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128, kwargs...)
     if ϕ(s_eff) > 0
          println("im in")
          ks, pks = xicalc(s -> 2*π^2*integral_on_mu_doppler(s_eff, s; kwargs...), L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          println("im out")
          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L ) .* pks
          end
     else
          return 0
     end
end

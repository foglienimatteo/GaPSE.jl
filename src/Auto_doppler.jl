
##########################################################################################92


function ξ_doppler(s1, s2, y)

     D1, D2 = D(s1), D(s2)

     f1, ℋ1 = f(s1), ℋ(s1)#, ℋ0
     f2, ℋ2 = f(s2), ℋ(s2)#, ℋ0
     #f1, ℋ1, ℋ1_p, s_b1 = f(s1), ℋ(s1), ℋ_p(s1), s_b(s1)
     #f2, ℋ2, ℋ2_p, s_b2 = f(s2), ℋ(s2), ℋ_p(s2), s_b(s2)
     #ℛ1, ℛ2 = ℛ(s1, ℋ1, ℋ1_p, s_b1), ℛ(s2, ℋ2, ℋ2_p, s_b2)
     ℛ1, ℛ2 = 1 - 1 / (ℋ1 * s1), 1 - 1 / (ℋ2 * s2)

     s = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     parenth = I00(s) / 45.0 + I20(s) / 31.5 + I40(s) / 105.0

     first = prefac * (c1 * parenth + I02(s) * y * s^2 / 3.0)

     return first

     #=
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
     =#
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


# mean time of evaluation: 141 seconds
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



function integrand_ξ_lensing(χ1, χ2, s1, s2, y)
     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)
     a(χ) = 1 / (1 + z_of_s(χ))

     D1, D2 = D(χ1), D(χ2)
     ℋ0 = ℋ(0)
     s_b1, s_b2 = s_b(s1), s_b(s2)
     psb1, psb2 = 5 * s_b1 - 2, 5 * s_b2 - 2

     denomin = s1 * s2 * Δχ^4 * a(χ1) * a(χ2)
     factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2

     new_J00 = -(3 / 4) * χ1^2 * χ2^2 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - 9 * χ1 * χ2 * y^2 - 7 * χ1 * χ2)
     new_J02 = -(3 / 2) * χ1^2 * χ2^2 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - 3 * χ1 * χ2 * y^2 - 5 * χ1 * χ2)
     new_J31 = 9 * y * Δχ^6
     new_J22 = (9 / 4) * χ1 * χ2 * (
                    2 * (χ1^4 + χ2^4) * (7y^2 - 3)
                    -
                    16y * χ1 * χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
                    +
                    χ1^2 * χ2^2 * (11y^4 + 14y^2 + 23)
               )


     return factor / denomin * (
          new_J00 * I00(Δχ) + new_J02 * I20(Δχ) +
          new_J31 * I13(Δχ) + new_J22 * I22(Δχ)
     )
end


function ξ_lensing(s1, s2, y; rtol = 1e-2, atol = 1e-3, kwargs...)
     my_int(var) = integrand_ξ_lensing(var[1], var[2], s1, s2, y)
     a = [0.0, 0.0]
     b = [s1, s2]
     return hcubature(my_int, a, b; rtol = rtol, atol = atol, kwargs...)[1]
end

function int_on_mu_lensing(s1, s, μ)
     if ϕ(s2(s1, s, μ)) > 0
          return ξ_lensing(s1, s2(s1, s, μ), y(s1, s, μ)) * spline_F(s / s1, μ)
     else
          return 0.0
     end
end

function integral_on_mu_lensing(s1, s; kwargs...)
     return quadgk(μ -> int_on_mu_lensing(s1, s, μ), -1, 1; atol=1e-2, rtol=1e-3, kwargs...)[1]
end


function map_integral_on_mu_lensing(s1 = s_eff; kwargs...)
     ss = 10 .^ range(-1, 3, length = 100)
     xis = [integral_on_mu_lensing(s1, s; kwargs...) for s in ss]
     return (ss, xis)
end



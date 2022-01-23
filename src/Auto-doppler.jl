

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
     s = √(s1^2 + s2^2 - 2 * s1 * s2 * y)
     
     D(s1) * D(s2) * ( J00(s1, s2, y) * I00(s) + J02(s1, s2, y) * I20(s) +
          J04(s1, s2, y) * I40(s) + J20(s1, s2, y) * I02(s) ) + 
     D(s1) * ( J31(s1, s2, y) * I13(s1) + J11(s1, s2, y) * I11(s1) + 
          J13(s1, s2, y) * I31(s1) ) +
     D(s2) * (J31(s2, s1, y) * I13(s2) + J11(s2, s1, y) * I11(s2) + 
          J13(s2, s1, y) * I31(s2)) +
     Jσ2(s1, s2, y) * σ_2
end





function PS_doppler(L::Integer = 0; int_s_min = 1e-2, int_s_max = 2*s_max, N = 1024, kwargs...)
     if ϕ(s_eff) > 0
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_μ(ξ_doppler, s_eff, s; L=L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1/int_s_max)
     
          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2) * s_eff^2) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L * s_eff^2) .* pks
          end
     else
          return 0
     end
end

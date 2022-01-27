

function integrand_ξ_lensing(χ1, χ2, s1, s2, y; tol=1e-2)

     (s(s1,s2,y)>tol*s1) || (return 0.0)

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     D1, D2 = D(χ1), D(χ2)
     a_χ1, a_χ2 = 1 / (1 + z_of_s(χ1)), 1 / (1 + z_of_s(χ2))
     #psb1, psb2 = -1.0, -1.0
     #s_b1, s_b2 = s_b(s1), s_b(s2)
     #psb1, psb2 = 5 * s_b1 - 2, 5 * s_b2 - 2

     denomin = s1 * s2 * Δχ^4 * a_χ1 * a_χ2
     factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2)
     #factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2
     χ1χ2 = χ1 * χ2

     new_J00 = -0.75 * χ1χ2^2 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
     new_J02 = -1.5 * χ1χ2^2 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
     new_J31 = 9 * y * Δχ^6
     new_J22 = 2.25 * χ1χ2 * (
                    2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
                    -
                    16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
                    +
                    χ1χ2^2 * (11y^4 + 14y^2 + 23)
               )


     return s1^2 * factor / denomin * (
          new_J00 * I00(Δχ) + new_J02 * I20(Δχ) +
          new_J31 * I13(Δχ) + new_J22 * I22(Δχ)
     )
end


function ξ_lensing(s1, s2, y; kwargs...)
     #=
     μs = -1:μ_step:1
     xs = x1:x_step:x2
     xs_grid = [x for x = xs for μ = μs]
     μs_grid = [μ for x = xs for μ = μs]

     time_1 = time()
     new_F(x, μ) = F(x, μ; kwargs...)
     Fs_grid = @showprogress map(new_F, xs_grid, μs_grid)
     time_2 = time()
     =#
     
     my_int(var) = integrand_ξ_lensing(var[1], var[2], s1, s2, y)
     a = [0.0, 0.0]
     b = [s1, s2]
     return hcubature(my_int, a, b; kwargs...)
end

function int_on_mu_lensing(s1, s, μ; L::Integer=0)
     if ϕ(s2(s1, s, μ)) > 0
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensing(s1, s2(s1, s, μ), y(s1, s, μ); rtol=1e-2, atol=1e-4)
          #println("int = $int")
          return int .* (spline_F(s / s1, μ) * Pl(μ, L))
     else
          return 0.0
     end
end

function integral_on_mu_lensing(s1, s; L::Integer = 0, kwargs...)
     int = quadgk(μ -> int_on_mu_lensing(s1, s, μ; L=L)[1], -1, 1; kwargs...)
     println("int = $int")
     return int
end


function map_integral_on_mu_lensing(s1 = s_eff; 
               L::Integer = 0, pr::Bool=true, kwargs...)

     t1 = time()
     ss = 10 .^ range(-1, 3, length = 100)
     vec = [integral_on_mu_lensing(s1, s; L=L, kwargs...) for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu_lensing [in s] = $(t2-t1)\n")
     return (ss, xis, xis_err)
end


function PS_lensing(; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true, kwargs...)

     if ϕ(s_eff) > 0
          t1 = time()
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_mu_lensing(s_eff, s; L = L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          t2 = time()
          pr && println("\ntime needed for PS_lensing [in s] = $(t2-t1)\n")

          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L) .* pks
          end
     else
          throw(ErrorException(
               "ϕ(s_eff) should be >0, but in s_eff=$s_eff " *
               "is ϕ(s_eff) = $(ϕ(s_eff))! \n"
          ))
     end
end



##########################################################################################92



function print_map_int_on_mu_lensing(out::String; L::Integer = 0,
     s1 = s_eff, pr::Bool = true, kwargs...)

     t1 = time()
     vec = map_integral_on_mu_lensing(s1; L = L, pr = pr, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of xi_lensing.")
          parameters_used(io)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for (i, key) in enumerate(keys(kwargs))
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\ns \t xi \t xi_error")
          for (s, xi, xi_err) in zip(vec[1], vec[2], vec[3])
               println(io, "$s \t $xi \t $(xi_err)")
          end
     end
end

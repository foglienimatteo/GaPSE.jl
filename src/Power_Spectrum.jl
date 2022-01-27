

IMPLEMENTED_GR_EFFECTS = ["auto_doppler", "auto_lensing"]
IMPLEMENTED_XI_FUNCS = [integral_on_mu_doppler, integral_on_mu_lensing]
dict_gr_xi = Dict([a => b for (a, b) in zip(IMPLEMENTED_GR_EFFECTS, IMPLEMENTED_XI_FUNCS)]...)


# mean time of evaluation: 2 seconds!!!
function PS(f_in::Union{Function,Dierckx.Spline1D};
     int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true)

     t1 = time()
     ks, pks = xicalc(x -> 2 * π^2 * f_in(x), L, 0; N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
     t2 = time()
     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(t2-t1)\n")

     if iseven(L)
          return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
     else
          return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L) .* pks
     end
end


function PS(in::String; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true)

     xi_table = readdlm(in, comments = true)
     ss = convert(Vector{Float64}, xi_table[2:end, 1])
     fs = convert(Vector{Float64}, xi_table[2:end, 2])
     f_in = Spline1D(ss, fs)

     return PS(f_in; int_s_min = int_s_min, int_s_max = int_s_max, N = N, L = L, pr = pr)
end


function print_PS(out::String, in::String;
     int_s_min = 1e-2, int_s_max = 2 * s_max, s1=s_eff,
     N = 128, L::Integer = 0, pr::Bool = true, kwargs...)

     time_1 = time()

     vec = if isfile(in)
               pr && println("\nI'm computiong the PS from the file $in")
               PS(in, int_s_min = int_s_min, int_s_max = int_s_max,
               N = N, L = L, pr = pr, kwargs...)
     else
          if in ∈ IMPLEMENTED_GR_EFFECTS
               pr && println("\nI'm computiong the PS for the $in GR effect.")
               t1 = time()
               ss = 10 .^ range(-1, 3, length = 100)
               v = [dict_gr_xi[in](s1, s; L = L, kwargs...) for s in ss]
               xis, xis_err = [x[1] for x in v], [x[2] for x in v]
               t2 = time()
               pr && println("\ntime needed to create the xi map [in s] = $(t2-t1)\n")
               f_in = Spline1D(ss, xis)
               PS(f_in, int_s_min = int_s_min, int_s_max = int_s_max,
                    N = N, L = L, pr = pr, kwargs...)
          else
               throw(ErrorException(
                    "$in is neither a GR impemented effect or a file.\n" *
                    "\t The implemented GR effects are currently: \n"*
                    "\t $(IMPLEMENTED_GR_EFFECTS)"
               ))
          end
     end

     time_2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is the Power Spectrum computation of ")
          parameters_used(io)
          println(io, "#\n# For this PS computation we set: ")
          println(io, "# \t int_s_min = $int_s_min \t int_s_max = $int_s_max ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for (i, key) in enumerate(keys(kwargs))
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\nk \t P")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t $pk")
          end
     end
end

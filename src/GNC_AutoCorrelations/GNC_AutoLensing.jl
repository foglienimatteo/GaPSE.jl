
function integrand_ξ_GNC_Lensing(
    IP1::Union{Point,DP}, IP2::Union{Point,DP},
    P1::Union{Point,DP}, P2::Union{Point,DP},
    y, cosmo::Union{Cosmology, DC}; Δχ_min::AbstractFloat=DevFloat(1e-1), 
    b1=nothing, b2=nothing, s_b1=nothing, s_b2=nothing, 𝑓_evo1=nothing, 𝑓_evo2=nothing,
    s_lim=nothing, obs::Union{Bool,Symbol}=:noobsvel) where {DP<:DevPoint, DC<:DevCosmology}

    
    s1 = P1.comdist
    s2 = P2.comdist
    
    χ1, D1, a1 = IP1.comdist, IP1.D, IP1.a
    χ2, D2, a2 = IP2.comdist, IP2.D, IP2.a

    Ω_M0 = cosmo.params.Ω_M0
    
    s_b_s1 = isnothing(s_b1) ? cosmo.params.s_b1 : s_b1
    s_b_s2 = isnothing(s_b2) ? cosmo.params.s_b2 : s_b2

    Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
    Δχ = Δχ_square > 0 ? √(Δχ_square) : 0

    denomin = s1 * s2 * a1 * a2
    factor = ℋ0^4 * Ω_M0^2 * D1 * (s1 - χ1) * D2 * (s2 - χ2) * (5 * s_b_s1 - 2) * (5 * s_b_s2 - 2)
    
    #if Δχ > Δχ_min
        
        #χ1χ2 = χ1 * χ2
        
        #new_J00 = -3 / 4 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
        #new_J02 = -3 / 2 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
        new_J31 = 9 * y * Δχ^2
        #new_J22 = 9 / 4 * χ1χ2 / Δχ^4 * (
        #    2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
        #    - 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
        #    + χ1χ2^2 * (11y^4 + 14y^2 + 23)
        #)
        
        #new_J00, new_J02, new_J22 = 1.0f0, 1.0f0, 1.0f0
        
        new_J00, new_J02, new_J22 = new_J31, new_J31, new_J31
        #=
        new_J00 = begin
            new_J00_a = 8 * y * (χ1^2 + χ2^2)
            new_J00_b = - χ1χ2 * (9 * y^2 + 7)
            new_J00_sum = new_J00_a + new_J00_b
            #eps(new_J00_a) ≈ abs.(new_J00_sum) ? 0.0 : -3 / 4 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * new_J00_sum
            eps(new_J00_a) ≈ abs.(new_J00_sum) ? zero(s1) : -3  * χ1χ2^2 / Δχ^4 * (y^2 - 1) * new_J00_sum /4
        end 

        new_J02 = begin
            new_J02_a = 4 * y * (χ1^2 + χ2^2)
            new_J02_b = -χ1χ2 * (3 * y^2 + 5)
            new_J02_sum = new_J02_a + new_J02_b
            #eps(new_J02_a) ≈ abs.(new_J02_sum) ? 0.0 : -3 / 2 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * new_J02_sum
            eps(new_J02_a) ≈ abs.(new_J02_sum) ? zero(s1) : -3  * χ1χ2^2 / Δχ^4 * (y^2 - 1) * new_J02_sum / 2
        end
      
        new_J22 = begin
            new_J22_a = 2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
            new_J22_b = - 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
            new_J22_c = χ1χ2^2 * (11y^4 + 14y^2 + 23)
            new_J22_sum = new_J22_a + new_J22_b + new_J22_c
            #eps(new_J22_b) ≈ abs.(new_J22_sum) ? 0.0 : 9 / 4 * χ1χ2 / Δχ^4 * new_J22_sum
            eps(new_J22_b) ≈ abs.(new_J22_sum) ? zero(s1) : 9 * χ1χ2 / Δχ^4 * new_J22_sum / 4
        end
        =#
        
        #new_J22 = log10(abs(new_J22_sum)) < log10(abs(new_J22_b)) - 15 ? 0.0 : new_J22_coeff * new_J22_sum
        
        #if(log10(abs(new_J22_sum)) < log10(abs(new_J22_b)) - 14.5)
        #    0.0
        #else
        #    new_J22_coeff * new_J22_sum
        #end
        
        #if(log(abs(new_J22_a)) > 13 && log(abs(new_J22_b)) > 13 && log(abs(new_J22_c)) > 13 && log(abs(new_J22_sum)) < )
        

        I00 = cosmo.tools.I00(Δχ)
        I20 = cosmo.tools.I20(Δχ)
        I13 = cosmo.tools.I13(Δχ)
        I22 = cosmo.tools.I22(Δχ)
 
        return factor / denomin * (
            new_J00 * I00 + new_J02 * I20 +
            new_J31 * I13 + new_J22 * I22
        )

        # The problem is new_J22: its parenthesis does not cancel out completely
        #  sometimes for y->1.0 and Δχ smalls (below 0.2 tipically)
        #if( abs(y-1.0)<1e-3 && abs(χ1 - 2356.267102533347) < 1  && abs(2356.1023012765418 - χ2) < 1)
            #println(
            #    """
            #    s1 = $s1 \t s2 = $s2 \t χ1 = $χ1 \t χ2 = $χ2 \t y = $y \t Δχ = $Δχ
            #    I00 = $I00 \t I20 = $I20 \t I13 = $I13 \t I22 = $I22
            #    new_J00 = $new_J00 \t new_J02 = $new_J02 \t new_J31 = $new_J31 \t new_J22 = $new_J22
            #
            #    χ1χ2 = $(χ1χ2) \t Δχ^4 = $(Δχ^4)
            #    new_J22_coeff = $(9 / 4 * χ1χ2 / Δχ^4)
            #    new_J22_first = $(2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)) \t log10(abs(new_J22_a)) = $(log10(abs(new_J22_a)))
            #    new_J22_second = $(- 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)) \t log10(abs(new_J22_b)) = $(log10(abs(new_J22_b)))
            #    new_J22_third = $(χ1χ2^2 * (11y^4 + 14y^2 + 23)) \t log10(abs(new_J22_c)) = $(log10(abs(new_J22_c)))
            #    new_J22_sum = $(new_J22_sum) \t log10(abs(new_J22_sum)) = $(log10(abs(new_J22_sum)))
            #    res = $res
            #    \n---
            #    """
            #)
        #end

        #=
    else
        #println("s1 = $s1 \t s2 = $s2")
        #println("χ1 = $χ1 \t χ2 = $χ2")
        #3 / 5 * (5 * cosmo.tools.σ_2 + 6 * cosmo.tools.σ_0 * χ2^2)
        #res = 3 * cosmo.tools.σ_2 + 6 / 5 * χ1^2 * cosmo.tools.σ_0
        #println("res = $res")
        #res

        #3 * cosmo.tools.σ_2 + 6 / 5 * χ1^2 * cosmo.tools.σ_0
        return factor / denomin * 3 * cosmo.tools.σ_2 + 6 * χ1^2 * cosmo.tools.σ_0 / 5
    end
    =#
end

function integrand_ξ_GNC_Lensing(
    χ1::AbstractFloat, χ2::AbstractFloat,
    s1::AbstractFloat, s2::AbstractFloat,
    y, cosmo::Cosmology;
    kwargs...)

    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
    return integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end




##########################################################################################92



function ξ_GNC_Lensing(P1::Union{Point,DP}, P2::Union{Point,DP}, y, cosmo::Cosmology;
    en::AbstractFloat=1e6, N_χs_2::Int=100, suit_sampling::Bool=true, 
    backend=CPU(),  kwargs...) where DP<:DevPoint

    χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
    #χ2s = P2.comdist .* range(1e-5, 1, length = N_χs_2 + 7)
    χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

    if backend==false

        IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
        IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

        int_ξs = [
          GaPSE.integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
          for IP1 in IP1s, IP2 in IP2s
        ]

        res = trapz((χ1s, χ2s), int_ξs)
        return res

    else
        #=
        # with this I get: "Argument 8 to your kernel function is of type GaPSE.Cosmology, which is not a bitstype"
        int_ξs = KernelAbstractions.zeros(backend, Float32, N_χs_2, N_χs_2)

        kernel! = kernel_2d!(backend)
        kernel!(int_ξs, GaPSE.integrand_ξ_GNC_Lensing, P1, P2, y, cosmo, N_χs_2, kwargs...; ndrange=size(int_ξs))
        KernelAbstractions.synchronize(backend)

        res = trapz((χ1s, χ2s), reshape(int_ξs, N_χs_2, N_χs_2))
        return res
        =#

        #= 
        # with this everything works as expected
        χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
        #χ2s = P2.comdist .* range(1e-5, 1, length = N_χs_2 + 7)
        χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

        IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
        IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

        int_ξs = [
            en * GaPSE.integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
            for IP1 in IP1s, IP2 in IP2s
        ]

        res = trapz((χ1s, χ2s), int_ξs)
        return res / en
        =#
        
        χ1s = P1.comdist .* range(1e-6, 1, length=N_χs_2)
        #χ2s = P2.comdist .* range(1e-5, 1, length = N_χs_2 + 7)
        χ2s = P2.comdist .* range(1e-6, 1, length=N_χs_2)

        #IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
        #IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

        int_ξs = KernelAbstractions.zeros(backend, Float32, N_χs_2, N_χs_2)
        #tmp = KernelAbstractions.zeros(length(IP1s))

        #i = 1
        #for j in 1:length(IP2s)
        #  int_ξ_Lensing[i,j] = @oneapi items=30 en * GaPSE.integrand_ξ_GNC_Lensing_oneapi(IP1s, IP2s[j], P1, P2, y, cosmo; kwargs...)
        #end
        #  global i+=1

        #function vadd(tmp, IP1s, IP2s;kwargs...)
        #   i = get_global_id()
        #   #@inbounds c[i] = a[i] + b[i]
        #   tmp[i] = en * GaPSE.integrand_ξ_GNC_Lensing(IP1s[i], IP2s, P1, P2, y, cosmo;kwargs...)
        #   return 
        #end

        #a = oneArray(rand(10));

        #b = oneArray(rand(10));

        #c = similar(a);
        #for j in 1:length(IP2s)
        # @oneapi items=10 vadd(tmp, IP1s, IP2s[j])
        # [int_ξs[i, j] =  tmp[i] for i in 1:legnth(IP1s)]
        #end
        #@kernel function mykernel!(int_ξs, integrand_ξ_GNC_Lensing, IP1s, IP2s, P1, P2, y, cosmo, kwargs...)
        #@kernel function mykernel!(int_ξs, integrand_ξ_GNC_Lensing, P1, P2, y, cosmo, kwargs...) 
        #  i, j = @index(Global, NTuple)
        #  IP1 = GaPSE.Point(P1.comdist * lr(1e-6, 1, N_χs_2, i), cosmo)
        #  IP2 = GaPSE.Point(P2.comdist * lr(1e-6, 1, N_χs_2, j), cosmo)
        #  #IP1 = GaPSE.Point(χ1s[i], cosmo)
        #  #IP2 = GaPSE.Point(χ2s[j], cosmo) 
        #  int_ξs[i,j] = integrand_ξ_GNC_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
        #end
        devcosmo = adapt(backend, cosmo)
        IP1s = adapt(backend, [adapt(backend, GaPSE.Point(x, cosmo)) for x in χ1s])
        IP2s = adapt(backend, [adapt(backend, GaPSE.Point(x, cosmo)) for x in χ2s])
        devIP1, devIP2 =  adapt(backend,P1), adapt(backend,P2)
        int_f(IP1, IP2, devIP1, devIP2, y, devcosmo) = integrand_ξ_GNC_Lensing(IP1, IP2, devIP1, devIP2, y, devcosmo)

        @kernel function mykernel!(int_f, int_ξs, IP1s, IP2s, devIP1, devIP2, y, devcosmo)
            i, j = @index(Global, NTuple)
            #IP1 = GaPSE.Point(P1.comdist * lr(1e-6, 1, N_χs_2, i), cosmo)
            #IP2 = GaPSE.Point(P2.comdist * lr(1e-6, 1, N_χs_2, j), cosmo)
            #IP1 = GaPSE.Point(χ1s[i], cosmo)
            #IP2 = GaPSE.Point(χ2s[j], cosmo) 
            int_ξs[i, j] = int_f(IP1s[i], IP2s[j], devIP1, devIP2, y, devcosmo)
        end
        #Array(a) .+ Array(b) == Array(c)

        #kernel!(int_ξs, GaPSE.integrand_ξ_GNC_Lensing, IP1s, IP2s, P1, P2, y, cosmo, kwargs...; ndrange=size(int_ξs))
        compiled_kernel! = mykernel!(backend, 64)
        compiled_kernel!(int_f, int_ξs, IP1s, IP2s, devIP1, devIP2, y, devcosmo; ndrange=size(int_ξs))
        KernelAbstractions.synchronize(backend)
        hostint_ξs = Array(int_ξs)

        res = trapz((χ1s, χ2s), reshape(hostint_ξs,N_χs_2, N_χs_2))
        #println("res = $res")
        return res
        
    end
end


function ξ_GNC_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
    P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
    return ξ_GNC_Lensing(P1, P2, y, cosmo; kwargs...)
end


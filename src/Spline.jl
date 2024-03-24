# This works, and it's the best one!
struct MySpline
    xs::Vector{Float64}
    coeffs::Vector{Tuple{Float64,Float64,Float64,Float64}}
    N::Int64

    function MySpline(xs::Vector{T1}, ys::Vector{T2}; bc::String="Error", ic::String="Natural") where {T1<:Real,T2<:Real}
        @assert length(xs) == length(ys) "xs and ys have different lengths!"
        @assert length(xs) > 3 "cannot interpolate with less than 4 pairs of data"

        # Boundary conditions
        @assert bc == "Error" || bc=="error" "Only bc=\"Error\"/\"error\" is allowed at the moment, bc=$bc is not valid."


        N, Δ = length(xs), diff(xs)
        γ, δ = similar(Δ), similar(xs)
        B, C, D = similar(xs), similar(xs), similar(xs)
        # B, C and D should be "similar(Δ)", but we want to use C[i] = g''(x[i]) to compute B and D, 
        # and their last terms B[N] & D[N] require g''(x[N]), that would not be included in C if
        # it would span from 1 to N-1; we will cut them at the end

        # Initial conditions
        b_1, b_N, c_1, a_N, d_1, d_N =
            if ic == "Natural"
                1, 1, 0, 0, 0, 0
            elseif ic == "Parabolic"
                1, 1, -1, -1, 0, 0
            elseif ic == "ThirdDerivative"
                η_in = [(ys[i+1] - ys[i]) / (xs[i+1] - xs[i]) for i in 1:3]
                η_in_2 = [(η_in[i+1] - η_in[i]) / (xs[i+2] - xs[i]) for i in 1:2]
                η_in_3 = (η_in_2[2] - η_in_2[1]) / (xs[4] - xs[1])

                η_end = [(ys[N-3+i] - ys[N-4+i]) / (xs[N-3+i] - xs[N-4+i]) for i in 1:3]
                η_end_2 = [(η_end[i+1] - η_end[i]) / (xs[N-2+i] - xs[N-4+i]) for i in 1:2]
                η_end_3 = (η_end_2[2] - η_end_2[1]) / (xs[N] - xs[N-3])

                -Δ[1] / 6, -Δ[N-1] / 6, Δ[1] / 6, Δ[N-1] / 6, Δ[1]^2 * η_in_3, -Δ[N-1]^2 * η_end_3
            else
                vic = ["Natural", "Parabolic", "ThirdDerivative"]
                throw(AssertionError("ic=$ic is not valid. Valid inital conditions are: $vic"))
            end

        γ[1], δ[1] = c_1 / b_1, d_1 / b_1
        for i in 2:N-1
            ζ = (ys[i+1] - ys[i]) / Δ[i] - (ys[i] - ys[i-1]) / Δ[i-1]
            den = 2 * (Δ[i] + Δ[i-1]) - Δ[i-1] * γ[i-1]
            γ[i] = Δ[i] / den
            δ[i] = (6 * ζ - Δ[i-1] * δ[i-1]) / den
        end
        δ[N] = (d_N - a_N * δ[N-1]) / (b_N - a_N * γ[N-1])

        @assert ic != "Natural" || δ[N] ≈ 0 " δ[N] = $(δ[N]) !=0  with ic=$ic"

        #β[N] = b_N * β[N-1] - a_N * γ[N-1]
        #δ[N] = (d_N * β[N-1] - a_N * δ[N-1])/β[N]

        C[N] = δ[N] / 2
        for i in N-1:-1:1
            C[i] = (δ[i] - 2 * γ[i] * C[i+1]) / 2
            B[i] = (ys[i+1] - ys[i]) / Δ[i] - Δ[i] * (C[i+1] + 2 * C[i]) / 3
            D[i] = (C[i+1] - C[i]) / (3 * Δ[i])
        end
        #C = C[begin:end-1] #not needed, the cycle inside new does that

        new(xs, [(ys[i], B[i], C[i], D[i]) for i in 1:N-1], N)
    end
end

function (S::MySpline)(x)
    i = searchsortedlast(S.xs, x)
    @assert 1 ≤ i ≤ S.N "BC Error: $(S.xs[1]) ≤ $x ≤ $(S.xs[end]) does not hold!"
    x_i = S.xs[i]
    #a, b, c, d = S.coeffs[i]
    #return a + b * (x - x_i) + c * (x - x_i)^2 + d * (x - x_i)^3
    #return @evalpoly(x-x_i, a, b, c, d)
    return @evalpoly(x - x_i, S.coeffs[i]...)
end

function derivative(S::MySpline, x::T; nu::Int=1) where {T<:Real}
    @assert nu > 0 "The derivative order nu must be >0; nu=$nu is not valid!"
    i = searchsortedlast(S.xs, x)
    @assert 1 ≤ i ≤ S.N "BC Error: $(S.xs[1]) ≤ $x ≤ $(S.xs[end]) does not hold!"
    x_i = S.xs[i]
    if nu > 3
        return 0
    else
        coeffs = S.coeffs[i]
        return @evalpoly(x - x_i, coeffs[nu+1:end] .* [factorial(i) / factorial(i - nu) for i in nu:length(coeffs)-1]...)
    end
end

# not efficient the following!
function derivative(S::MySpline, x::Vector{T}; nu::Int=1) where {T<:Real}
    return [derivative(S, p; nu=nu) for p in x]
end

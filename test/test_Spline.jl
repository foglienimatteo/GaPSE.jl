

@testset "test Spline constructor" begin
    @test_throws MethodError GaPSE.MySpline([x for x in 1:5], [1.0 * x for x in 1:5]; nu=2)
    @test_throws AssertionError GaPSE.MySpline([1.0*x for x in 1:5], [1.0*x for x in 1:6])
    @test_throws MethodError GaPSE.MySpline([1.0*x for x in 1:5], [1.0*x for x in 1:6]; nu=2)
    @test_throws AssertionError GaPSE.MySpline([1.0*x for x in 1:2], [1.0*x for x in 1:2])
    @test_throws MethodError GaPSE.MySpline([1.0*x for x in 1:3], [1.0*x for x in 1:3]; nu=2)
    @test_throws AssertionError GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; bc="idk")
    @test_throws AssertionError GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; ic="idk")

    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]) isa Any
    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; bc="error") isa Any
    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; bc="Error") isa Any
    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; ic="Natural") isa Any
    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; ic="Parabolic") isa Any
    @test GaPSE.MySpline([1.0*x for x in 1:6], [1.0*x for x in 1:6]; ic="ThirdDerivative") isa Any
end

@testset "test Spline evaluation" begin
    @testset "zeros" begin 
        xs = [1.0 * x for x in 1:10]
        best_sp_1 = GaPSE.MySpline(xs, xs)

        @test_throws AssertionError best_sp_1(11.0)
        @test_throws AssertionError best_sp_1(-1.0)
        @test best_sp_1(10.0) ≈ 10.0
        @test best_sp_1(10) ≈ 10
        @test best_sp_1(1.0) ≈ 1.0
        @test best_sp_1(5.0) ≈ 5.0
        @test best_sp_1(7) ≈ 7
        @test best_sp_1(4.32) ≈ 4.32

        @test_throws AssertionError best_sp_1.([x for x in 1:0.1:11])
        @test_throws AssertionError best_sp_1.([x for x in -1:0.1:10])
        @test all(best_sp_1.([x for x in 1:10]) .≈ xs)
        @test all(best_sp_1.(xs) .≈ xs)
    end

    @testset "test Spline evaluation - linear range" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        #xs = [1.0*x for x in 10 .^ range(0, 2, length=1000) ] 
        xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1];

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative");

        dierckx_ys = dierckx_sp.(rs[N1:N2])
        my_sp_1_ys = best_sp_1.(rs[N1:N2])
        my_sp_2_ys = best_sp_2.(rs[N1:N2])
        my_sp_3_ys = best_sp_3.(rs[N1:N2])

        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])

    end

    @testset "test Spline evaluation - log range" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        xs = [1.0*x for x in 10 .^ range(0, 2, length=N) ] 
        #xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = dierckx_sp.(rs[N1:N2])
        my_sp_1_ys = best_sp_1.(rs[N1:N2])
        my_sp_2_ys = best_sp_2.(rs[N1:N2])
        my_sp_3_ys = best_sp_3.(rs[N1:N2])

        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])
    end
end 


@testset "test Spline derivative" begin
    @testset "zeros" begin 
        xs = [1.0 * x for x in 1:10]
        best_sp_1 = GaPSE.MySpline(xs, xs)

        @test_throws AssertionError GaPSE.derivative(best_sp_1, 11.0)
        @test_throws AssertionError GaPSE.derivative(best_sp_1, -1.0)
        @test_throws AssertionError GaPSE.derivative(best_sp_1, 9.0; nu=0)
        @test GaPSE.derivative(best_sp_1, 10.0) ≈ 1.0
        @test GaPSE.derivative(best_sp_1, 10) ≈ 1
        @test GaPSE.derivative(best_sp_1, 1.0) ≈ 1.0
        @test GaPSE.derivative(best_sp_1, 7) ≈ 1
        @test GaPSE.derivative(best_sp_1, 4.32) ≈ 1.0

        @test_throws AssertionError GaPSE.derivative(best_sp_1, [x for x in 1:0.1:11])
        @test_throws AssertionError GaPSE.derivative(best_sp_1, [x for x in -1:0.1:10])
        @test all(GaPSE.derivative(best_sp_1, xs) .≈ 1.0)
        @test all(GaPSE.derivative(best_sp_1, [x for x in 1:10]) .≈ 1.0)
    end

    @testset "test Spline derivative - linear range - nu=1" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        #xs = [1.0*x for x in 10 .^ range(0, 2, length=1000) ] 
        xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2])

        # d^1y/dx^1 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=1))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=1)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2])
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=1)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])

    end

    @testset "test Spline derivative - linear range - nu=2" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        #xs = [1.0*x for x in 10 .^ range(0, 2, length=1000) ] 
        xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=2)

        # d^2y/dx^2 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=2))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=2)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=2)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=2)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])

    end

    @testset "test Spline derivative - linear range - nu=3" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        #xs = [1.0*x for x in 10 .^ range(0, 2, length=1000) ] 
        xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=3)

        # d^3y/dx^3 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=3))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=3)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=3)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=3)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])

    end

    @testset "test Spline derivative - linear range - nu=4,5,6" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        #xs = [1.0*x for x in 10 .^ range(0, 2, length=1000) ] 
        xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        #dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        #dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=3)
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=4)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=5)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=6)

        @test all(my_sp_1_ys .≈ 0.0)
        @test all(my_sp_2_ys .≈ 0.0)
        @test all(my_sp_3_ys .≈ 0.0)
    end

    @testset "test Spline derivative - log range - nu=1" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        xs = [1.0 * x for x in 10 .^ range(0, 2, length=N)]
        #xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=1)

        # d^1y/dx^1 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=1))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=1)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=1)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=1)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])
    end

    @testset "test Spline derivative - log range - nu=2" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        xs = [1.0 * x for x in 10 .^ range(0, 2, length=N)]
        #xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=2)

        # d^2y/dx^2 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=2))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=2)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=2)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=2)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])
    end

    @testset "test Spline derivative - log range - nu=3" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 15, 15
        #lp,hp=0.005,0.005

        xs = [1.0 * x for x in 10 .^ range(0, 2, length=N)]
        #xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=3)

        # d^3y/dx^3 crosses zero inside the interval, and next to a zero a purely
        # relative tolerance is meaningless: the two splines agree to ~1e-5 in absolute
        # terms, but |y1-y2|/|y1| is unbounded as y1 -> 0. The floor below is the typical
        # size of this derivative over the interval that is actually sampled.
        ATOL = RTOL * maximum(abs, Dierckx.derivative(dierckx_sp,
            collect(range(xs[N1], xs[N2+1], length=101)); nu=3))
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=3)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=3)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=3)

        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_1_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_2_ys)])
        @test all([isapprox(y1, y2; rtol=RTOL, atol=ATOL) for (y1, y2) in zip(dierckx_ys, my_sp_3_ys)])
    end

    @testset "test Spline derivative - log range - nu=4,5,6" begin
        RTOL = 1e-4
        N = 10000
        N1, N2 = 10, 10
        #lp,hp=0.005,0.005

        xs = [1.0 * x for x in 10 .^ range(0, 2, length=N)]
        #xs = [1.0 * x for x in 1:1:N]
        ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

        #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

        # fixed seed: rs samples a random point inside each knot interval, so without
        # it this testset checks a different point on every run
        Random.seed!(20250921)
        rs = [xs[i] + (xs[i+1] - xs[i]) * rand() for i in 1:length(xs)-1]

        #dierckx_sp = Dierckx.Spline1D(xs, ys; bc="error")
        best_sp_1 = GaPSE.MySpline(xs, ys, ic="Natural")
        best_sp_2 = GaPSE.MySpline(xs, ys, ic="Parabolic")
        best_sp_3 = GaPSE.MySpline(xs, ys, ic="ThirdDerivative")

        #dierckx_ys = Dierckx.derivative(dierckx_sp, rs[N1:N2]; nu=3)
        my_sp_1_ys = GaPSE.derivative(best_sp_1, rs[N1:N2]; nu=4)
        my_sp_2_ys = GaPSE.derivative(best_sp_2, rs[N1:N2]; nu=5)
        my_sp_3_ys = GaPSE.derivative(best_sp_3, rs[N1:N2]; nu=6)

        @test all(my_sp_1_ys .≈ 0.0)
        @test all(my_sp_2_ys .≈ 0.0)
        @test all(my_sp_3_ys .≈ 0.0)
    end
end




@testset "test Spline - linear range" begin
    RTOL = 1e-4
    N = length(xs)
    N1, N2 = 10, 10
    #lp,hp=0.005,0.005

    #xs = [x for x in 10 .^ range(0, 2, length=1000) ] 
    xs = [1.0 * x for x in 1:1:10000]
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

@testset "test Spline - log range" begin
    RTOL = 1e-4
    N = length(xs)
    N1, N2 = 10, 10
    #lp,hp=0.005,0.005

    xs = [x for x in 10 .^ range(0, 2, length=10000) ] 
    #xs = [1.0 * x for x in 1:1:10000]
    ys = sin.(xs / 987.452) + cos.(xs .^ 9.32) .* (xs) .^ (0.39761)

    #N1, N2 = Int64(ceil(N*lp)), N-Int64(floor(N*hp))

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


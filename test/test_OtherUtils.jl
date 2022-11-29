# -*- encoding: utf-8 -*-
#
# This file is part of GaPSE
# Copyright (C) 2022 Matteo Foglieni
#
# GaPSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GaPSE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#


@testset "test warning" begin
     @test "WARNING: ciao\n" == @capture_out GaPSE.warning("ciao")
end


@testset "test check_compatible_dicts" begin
     REF = Dict(
          :N => 1024::Integer, :fit_min => 0.05::Float64,
          :fit_max => 0.5::Float64, :con => true::Bool, :name => "file_name.txt"::String)

     var = 3.0
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:ma => 1.0, :con => true), "x")
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(var => 1.0))
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict("con" => true))

     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:con => 1.0), "x")
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:con => true, :N => 3.14))
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:ma => 1.0), "x")

     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1.0 + 3 * im), "x")

     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1, :fit_max => 17.43), "x"))
     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1e-2, :fit_max => 17)))
     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict()))
end


@testset "test my_println_vec" begin
     vec = [x for x in 1:0.1:4]

     out_1 = @capture_out begin
          GaPSE.my_println_vec(vec, "vector"; N=8)
     end

     @test out_1 == "vector = [\n" *
                    "1.0 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , \n" *
                    "1.8 , 1.9 , 2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , \n" *
                    "2.6 , 2.7 , 2.8 , 2.9 , 3.0 , 3.1 , 3.2 , 3.3 , \n" *
                    "3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 , 4.0 , \n" *
                    "];\n"

     out_2 = @capture_out begin
          GaPSE.my_println_vec(vec, "vector"; N=3)
     end
     @test out_2 == "vector = [\n" *
                    "1.0 , 1.1 , 1.2 , \n" *
                    "1.3 , 1.4 , 1.5 , \n" *
                    "1.6 , 1.7 , 1.8 , \n" *
                    "1.9 , 2.0 , 2.1 , \n" *
                    "2.2 , 2.3 , 2.4 , \n" *
                    "2.5 , 2.6 , 2.7 , \n" *
                    "2.8 , 2.9 , 3.0 , \n" *
                    "3.1 , 3.2 , 3.3 , \n" *
                    "3.4 , 3.5 , 3.6 , \n" *
                    "3.7 , 3.8 , 3.9 , \n" *
                    "4.0 , \n" *
                    "];\n"
end

@testset "test parent_directory" begin
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads/file.txt") == "/Users/matteofoglieni/Downloads/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads/") == "/Users/matteofoglieni/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads") == "/Users/matteofoglieni/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/") == "/Users/"
     @test GaPSE.parent_directory("/Users/matteofoglieni") == "/Users/"
     @test GaPSE.parent_directory("/Users/") == "/"
     @test GaPSE.parent_directory("/Users") == "/"

     @test GaPSE.parent_directory("matteofoglieni/Downloads/file.txt") == "matteofoglieni/Downloads/"
     @test GaPSE.parent_directory("matteofoglieni/Downloads/") == "matteofoglieni/"
     @test GaPSE.parent_directory("matteofoglieni/Downloads") == "matteofoglieni/"
     @test GaPSE.parent_directory("matteofoglieni/") == "./"
     @test GaPSE.parent_directory("matteofoglieni") == "./"
     @test GaPSE.parent_directory("file.txt") == "./"
end

@testset "test check_parent_directory" begin
     @test_throws AssertionError GaPSE.check_parent_directory("/Users/matteofoglieni/notadirectory/file.txt")
     @test_throws AssertionError GaPSE.check_parent_directory("notadirectory/file.txt")

     #@test isnothing(GaPSE.check_parent_directory("/Users/matteofoglieni/notadirectory/"))
     @test isnothing(GaPSE.check_parent_directory("notadirectory/"))
     #@test isnothing(GaPSE.check_parent_directory("/Users/"))
     @test isnothing(GaPSE.check_parent_directory("/"))
end

@testset "test return_namefile" begin
     @test_throws AssertionError GaPSE.return_namefile("/Users/matteofoglieni/Downloads/")
     @test_throws AssertionError GaPSE.return_namefile("Downloads/")
     @test_throws AssertionError GaPSE.return_namefile("./Downloads/")

     @test GaPSE.return_namefile("file") == "file"
     @test GaPSE.return_namefile("file.boh") == "file.boh"
     @test GaPSE.return_namefile("/Users/matteo.foglieni/ciao.file") == "ciao.file"
     @test GaPSE.return_namefile("matteo.foglieni/ciao.file.boh") == "ciao.file.boh"
     @test GaPSE.return_namefile("/Users/matteofoglieni/Downloads/file.txt") == "file.txt"
     @test GaPSE.return_namefile("./file.txt") == "file.txt"
     @test GaPSE.return_namefile("./file.dat") == "file.dat"
     @test GaPSE.return_namefile("file.txt") == "file.txt"
     @test GaPSE.return_namefile("file.dat") == "file.dat"
end


@testset "test check_namefile" begin
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file")
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file.boh")
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteo.foglieni/ciao.file")
     @test_throws AssertionError GaPSE.check_namefile("matteo.foglieni/ciao.file.boh")
     @test_throws AssertionError GaPSE.check_namefile("file")
     @test_throws AssertionError GaPSE.check_namefile("file.boh")
     @test_throws AssertionError GaPSE.check_namefile("./file")
     @test_throws AssertionError GaPSE.check_namefile("./file.boh")

     @test isnothing(GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file.txt"))
     @test isnothing(GaPSE.check_namefile("./file.txt"))
     @test isnothing(GaPSE.check_namefile("./file.dat"))
     @test isnothing(GaPSE.check_namefile("file.txt"))
     @test isnothing(GaPSE.check_namefile("file.dat"))
end

@testset "test check_group" begin
     @test_throws AssertionError GaPSE.check_group("prova")
     @test_throws AssertionError GaPSE.check_group("ld")
     @test_throws AssertionError GaPSE.check_group("generical")

     @test isnothing(GaPSE.check_group("LD"))
     @test isnothing(GaPSE.check_group("GNC"))
     @test isnothing(GaPSE.check_group("GNCxLD"))
     @test isnothing(GaPSE.check_group("LDxGNC"))
     @test isnothing(GaPSE.check_group("generic"))
end

@testset "test check_fileisingroup" begin
     name = "test_check_fileisingroup.txt"
     isfile(name) && rm(name)
     open(name, "w") do io
          println(io, "# boh")
          for i in 1:10
               println(io, "$i $(2*i) $(3*i)")
          end
     end

     @test_throws AssertionError GaPSE.check_fileisingroup(name, "prova")
     @test_throws AssertionError GaPSE.check_fileisingroup(name, "LD")
     @test isnothing(GaPSE.check_fileisingroup(name, "generic"))
     rm(name)

     open(name, "w") do io
          println(io, "# boh")
          for i in 1:10
               println(io, "$i $(2*i) $(3*i)")
          end
          println(io, "13 14 ")
     end

     @test_throws AssertionError GaPSE.check_fileisingroup(name, "generic")
     rm(name)
end


##########################################################################################92


@testset "test number_to_string" begin
     @test GaPSE.number_to_string(3) == "3"
     @test GaPSE.number_to_string(-2.15) == "-2.15"
     @test GaPSE.number_to_string(2.15 * im) == "2.15im"
     @test GaPSE.number_to_string(0.0 + 2.15 * im) == "2.15im"
     @test GaPSE.number_to_string(-0.0 - 2.15 * im) == "-2.15im"
     @test GaPSE.number_to_string(-3.1415 - 2.15 * im) == "-3.1415-2.15im"
     @test GaPSE.number_to_string(3.1415 + 2.15 * im) == "3.1415+2.15im"
end




##########################################################################################92


@testset "test readxy" begin
     open("some_data.txt", "w") do io
          println(io, "# lets create a small comment")
          println(io, "# made of two lines")
          for i in 1:15
               println(io, " \t $i \t\t $(i+2.5)")
          end
     end

     xs, ys = GaPSE.readxy("some_data.txt", comments=true)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2.5
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2)")
          end
     end

     xs, ys = GaPSE.readxy("some_data.txt", comments=false)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2)")
          end
     end

     xs, ys = GaPSE.readxy("some_data.txt"; xdt=Int, ydt=Int, comments=false)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2
     end

     rm("some_data.txt")
end




@testset "test readxall" begin
     open("some_data.txt", "w") do io
          println(io, "# lets create a small comment")
          println(io, "# made of two lines")
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys = GaPSE.readxall("some_data.txt", comments=true)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[1][i] ≈ i + 2
          @test ys[2][i] ≈ 0.1
          @test ys[3][i] ≈ 0.1 * i
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys = GaPSE.readxall("some_data.txt", comments=false)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[1][i] ≈ i + 2
          @test ys[2][i] ≈ 0.1
          @test ys[3][i] ≈ 0.1 * i
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys = GaPSE.readxall("some_data.txt"; xdt=Int, ydt=Float64, comments=false)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[1][i] ≈ i + 2
          @test ys[2][i] ≈ 0.1
          @test ys[3][i] ≈ 0.1 * i
     end

     rm("some_data.txt")
     println(ys)
end


@testset "test readxyall" begin
     open("some_data.txt", "w") do io
          println(io, "# lets create a small comment")
          println(io, "# made of two lines")
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys, zs = GaPSE.readxyall("some_data.txt", comments=true)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2
          @test zs[1][i] ≈ 0.1
          @test zs[2][i] ≈ 0.1 * i
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys, zs = GaPSE.readxyall("some_data.txt", comments=true)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2
          @test zs[1][i] ≈ 0.1
          @test zs[2][i] ≈ 0.1 * i
     end

     rm("some_data.txt")

     ################

     open("some_data.txt", "w") do io
          for i in 1:15
               println(io, " \t $i \t\t $(i+2) \t 0.1 \t $(0.1*i)")
          end
     end

     xs, ys, zs = GaPSE.readxyall("some_data.txt"; xdt=Int64, ydt=Float32, zdt=Float64, comments=true)
     for i in 1:15
          @test xs[i] ≈ i
          @test ys[i] ≈ i + 2
          @test zs[1][i] ≈ 0.1
          @test zs[2][i] ≈ 0.1 * i
     end

     rm("some_data.txt")
     println(ys)
end


##########################################################################################92


@testset "test sample_subdivision_begin" begin

     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 0.9, 2.0)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 1.1)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=0.0)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=1.0)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=-0.2)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=1.3)
     @test_throws AssertionError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; N=3)
     @test_throws ArgumentError GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=1.0, ass=false)

     vec_1 = GaPSE.sample_subdivision_begin(1.0, 1.2, 2.0; frac_begin=0.9, N=100)
     res_vec_1 = [1.0, 1.0022222222222221, 1.0044444444444445, 1.0066666666666666,
          1.008888888888889, 1.011111111111111, 1.0133333333333334, 1.0155555555555555,
          1.0177777777777777, 1.02, 1.0222222222222221, 1.0244444444444445, 1.0266666666666666,
          1.028888888888889, 1.031111111111111, 1.0333333333333334, 1.0355555555555556,
          1.0377777777777777, 1.04, 1.0422222222222222, 1.0444444444444445, 1.0466666666666666,
          1.048888888888889, 1.051111111111111, 1.0533333333333332, 1.0555555555555556,
          1.0577777777777777, 1.06, 1.0622222222222222, 1.0644444444444445, 1.0666666666666667,
          1.068888888888889, 1.0711111111111111, 1.0733333333333333, 1.0755555555555556,
          1.0777777777777777, 1.08, 1.0822222222222222, 1.0844444444444445, 1.0866666666666667,
          1.0888888888888888, 1.0911111111111111, 1.0933333333333333, 1.0955555555555556,
          1.0977777777777777, 1.1, 1.1022222222222222, 1.1044444444444443, 1.1066666666666667,
          1.1088888888888888, 1.1111111111111112, 1.1133333333333333, 1.1155555555555556,
          1.1177777777777778, 1.12, 1.1222222222222222, 1.1244444444444444, 1.1266666666666667,
          1.1288888888888888, 1.1311111111111112, 1.1333333333333333, 1.1355555555555557,
          1.1377777777777778, 1.14, 1.1422222222222222, 1.1444444444444444, 1.1466666666666667,
          1.1488888888888888, 1.1511111111111112, 1.1533333333333333, 1.1555555555555554,
          1.1577777777777778, 1.16, 1.1622222222222223, 1.1644444444444444, 1.1666666666666667,
          1.1688888888888889, 1.1711111111111112, 1.1733333333333333, 1.1755555555555555,
          1.1777777777777778, 1.18, 1.1822222222222223, 1.1844444444444444, 1.1866666666666668,
          1.1888888888888889, 1.191111111111111, 1.1933333333333334, 1.1955555555555555,
          1.1977777777777778, 1.2, 1.28, 1.36, 1.44, 1.52, 1.6, 1.68, 1.76, 1.84, 1.92, 2.0]

     @test length(vec_1) == length(res_vec_1)
     for (v1, v2) in zip(vec_1, res_vec_1)
          @test v1 ≈ v2
     end
end


@testset "test sample_subdivision_middle" begin

     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 0.9, 1.4, 2.0)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.2, 1.1)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 1.3)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=0.0)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=1.0)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=-0.2)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=1.3)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; rel_frac_begin=0.0)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; rel_frac_begin=1.0)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; rel_frac_begin=-0.2)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; rel_frac_begin=1.3)
     @test_throws AssertionError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; N=5)
     @test_throws MethodError GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_begin=1.0, ass=false)

     vec_1 = GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=0.70, N=100, rel_frac_begin=nothing)
     res_vec_1 = [1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.2028571428571428,
          1.2057142857142857, 1.2085714285714286, 1.2114285714285715, 1.2142857142857142,
          1.217142857142857, 1.22, 1.2228571428571429, 1.2257142857142858, 1.2285714285714286,
          1.2314285714285715, 1.2342857142857142, 1.237142857142857, 1.24, 1.2428571428571429,
          1.2457142857142858, 1.2485714285714287, 1.2514285714285713, 1.2542857142857142,
          1.2571428571428571, 1.26, 1.262857142857143, 1.2657142857142858, 1.2685714285714285,
          1.2714285714285714, 1.2742857142857142, 1.2771428571428571, 1.28, 1.282857142857143,
          1.2857142857142858, 1.2885714285714285, 1.2914285714285714, 1.2942857142857143,
          1.2971428571428572, 1.3, 1.302857142857143, 1.3057142857142856, 1.3085714285714285,
          1.3114285714285714, 1.3142857142857143, 1.3171428571428572, 1.32, 1.322857142857143,
          1.3257142857142856, 1.3285714285714285, 1.3314285714285714, 1.3342857142857143,
          1.3371428571428572, 1.34, 1.3428571428571427, 1.3457142857142856, 1.3485714285714285,
          1.3514285714285714, 1.3542857142857143, 1.3571428571428572, 1.36, 1.3628571428571428,
          1.3657142857142857, 1.3685714285714285, 1.3714285714285714, 1.3742857142857143,
          1.3771428571428572, 1.38, 1.3828571428571428, 1.3857142857142857, 1.3885714285714286,
          1.3914285714285715, 1.3942857142857144, 1.3971428571428572, 1.4, 1.4260869565217391,
          1.4521739130434783, 1.4782608695652173, 1.5043478260869565, 1.5304347826086957,
          1.5565217391304347, 1.5826086956521739, 1.608695652173913, 1.6347826086956523,
          1.6608695652173913, 1.6869565217391305, 1.7130434782608697, 1.7391304347826086,
          1.7652173913043478, 1.791304347826087, 1.817391304347826, 1.8434782608695652,
          1.8695652173913044, 1.8956521739130434, 1.9217391304347826, 1.9478260869565218,
          1.9739130434782608, 2.0]

     @test length(vec_1) == length(res_vec_1)
     for (v1, v2) in zip(vec_1, res_vec_1)
          @test v1 ≈ v2
     end

     vec_2 = GaPSE.sample_subdivision_middle(1.0, 1.2, 1.4, 2.0; frac_middle=0.70, N=30, rel_frac_begin=0.7)
     res_vec_2 = [1.0, 1.0285714285714285, 1.0571428571428572, 1.0857142857142856, 1.1142857142857143,
          1.1428571428571428, 1.1714285714285715, 1.2, 1.2095238095238094, 1.2190476190476192,
          1.2285714285714286, 1.2380952380952381, 1.2476190476190476, 1.2571428571428571,
          1.2666666666666666, 1.276190476190476, 1.2857142857142858, 1.2952380952380953,
          1.3047619047619048, 1.3142857142857143, 1.3238095238095238, 1.3333333333333333,
          1.3428571428571427, 1.3523809523809525, 1.361904761904762, 1.3714285714285714,
          1.380952380952381, 1.3904761904761904, 1.4, 1.6, 1.8, 2.0]
     @test length(vec_2) == length(res_vec_2)
     for (v1, v2) in zip(vec_2, res_vec_2)
          @test v1 ≈ v2
     end
end


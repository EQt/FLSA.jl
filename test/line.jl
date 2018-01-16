module TestLine
using Base.Test
if VERSION >= v"0.7-" using Random end
using FLSA

@testset "line" begin
    @testset "3 line knot" begin
        y = [1.0; 2.0; 3.5]
        x = FLSA.dp_line_naive(y, 1.0)
        @test x ≈ [2.0; 2.0; 2.5]
    end

    @testset "line: 3 knots, fast" begin
        y = [1.0, 2.0, 3.5]
        x = FLSA.dp_line(y, 1.0)
        @test x ≈ [2.0; 2.0; 2.5]
    end

    @testset "line: 5 knots same" begin
        y = [1.0; 4.0; 3.5; 3.8; 4.1]
        z = FLSA.dp_line_naive(y, 1.0)
        x = FLSA.dp_line(y, 1.0)
        @test x ≈ z
    end

    @testset "line: 3 rand knots same" begin
        y = rand(srand(13), 20)
        z = FLSA.dp_line_naive(y, 0.3)
        x = FLSA.dp_line(y, 0.3)
        d = 3
        @test round.(x,d) ≈ round.(z,d)
    end
end
end

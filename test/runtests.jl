using Test
using SabbatucciMACs

@testset "All" begin
    @testset "Edge Energies" begin
        @test edgeenergy(39, 17) == 3.35
        @test hasedge(39,17)
        @test edgeenergy(45,9) == 307.2
        @test hasedge(45, 9)
        @test edgeenergy(4,2) == 4.342
        @test hasedge(4,2)
        @test_throws KeyError edgeenergy(98,28)
        @test !hasedge(98,28)
        @test edgeenergy(12,5) == 3.986
        @test hasedge(12, 5)
        @test_throws KeyError edgeenergy(12,6)
        @test !hasedge(12, 6)
        @test edgeenergy(34, 6) == 166.5
        @test edgeenergy(23, 5) == 66.3
        @test edgeenergy(99, 29) == 3.0
        @test edgeenergy(99, 1) == 1.38436e5
    end

    @testset "MAC" begin
        @test isapprox(SabbatucciMACs.mac(29, 3000.0), 715, atol=1.0)
        @test isapprox(SabbatucciMACs.mac(92, 30000.0), 38.7, atol=0.1)
        @test isapprox(SabbatucciMACs.mac(1, 100.0), 11616, atol=1.0)
        @test isapprox(SabbatucciMACs.mac(10, 1200.0), 4719, atol=1.0)
        @test isapprox(SabbatucciMACs.mac(99, 100_000.0), 2.08, atol=0.01)
        @test isapprox(SabbatucciMACs.mac(99, 10.0), 20079, atol=1.0)
    end

    @testset "Misc" begin
        @test eachelement() == Base.OneTo(99)
        @test atomicweight(26)==55.8452
        @test atomicweight(92)==238.028913
        @test atomicweight(13)==26.98153857
        @test eachedge(26) == Set{Int}( ( 7, 4, 9, 10, 2, 3, 5, 8, 6, 1) )
        @test eachedge(92) == Set{Int}( ( 2, 11, 7, 9, 25, 29, 8, 3, 20, 14, 12, 18, 16, 21, 26, 10, 19, 17, 22, 6, 24, 4, 5, 13, 27, 15, 1 ))
    end
end

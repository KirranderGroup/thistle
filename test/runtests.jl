using Test

include(joinpath(@__DIR__, "..", "src", "molpro_input.jl"))

@testset "Geometry allocation" begin
    x_pos, y_pos, z_pos, atomname = allocate_geometry_arrays(3)

    @test eltype(x_pos) == Float64
    @test size(x_pos) == (3,)
    @test size(y_pos) == (3,)
    @test size(z_pos) == (3,)
    @test length(atomname) == 3

    atomname[2] = "He"
    x_pos[2], y_pos[2], z_pos[2] = (1.0, 2.0, 3.0)

    @test atomname[2] == "He"
    @test x_pos[2] == 1.0
    @test y_pos[2] == 2.0
    @test z_pos[2] == 3.0
end

@testset "MO matrix builder" begin
    coeffs = collect(1.0:6.0)
    mo = build_mo_matrix(coeffs, 2, 3)

    @test mo == [1.0 3.0 5.0; 2.0 4.0 6.0]
    @test eltype(mo) == Float64

    @test_throws ArgumentError build_mo_matrix(coeffs, 3, 3)
end


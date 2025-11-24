using Test

include("../src/molpro_input.jl")

@testset "Geometry builder" begin
    lines = ["H 0.0 1.0 2.0", "He 3.0 4.0 5.0"]
    geom = build_geometry(lines)
    @test geom.atomname == ["H", "He"]
    @test geom.x_pos == [0.0, 3.0]
    @test geom.y_pos == [1.0, 4.0]
    @test geom.z_pos == [2.0, 5.0]
end

@testset "MO matrix builder" begin
    coeffs = 1.0:6.0
    mo_numbers = [1, 1, 1, 2, 2, 2]
    real_indices = [2, 1, 3]
    mat = build_mo_matrix(coeffs, mo_numbers, real_indices)
    @test size(mat) == (2, 3)
    @test mat == [3.0 1.0 5.0; 4.0 2.0 6.0]
end

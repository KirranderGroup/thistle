using Test

include(joinpath(@__DIR__, "..", "src", "Integrals.jl"))
include(joinpath(@__DIR__, "..", "src", "MainCalculation.jl"))

include(joinpath(@__DIR__, "..", "src", "molpro_input.jl"))
include(joinpath(@__DIR__, "..", "src", "CIUtils.jl"))

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

@testset "tot_integral_k_ijkr control flow" begin
    mu = [0.1, 0.2]
    l = m = n = fill(0, 2)
    group_start = [1, 2]
    group_count = [1, 1]
    hx, hy, hz = 0.05, 0.02, 0.01
    h = hypot(hx, hypot(hy, hz))
    ones3 = ones(1, 1, 1)
    zcontr = fill(0.5, 1, 1, 1, 1)
    zcontr2 = fill(0.25, 1, 1, 1, 1)

    res = Integrals.tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                                        hx, hy, hz, h,
                                        ones3, ones3, ones3, ones3, ones3, ones3,
                                        1, 1, 2, 2, zcontr, zcontr2;
                                        cutoff1=1e-14, cutoff2=1e-14)

    expected_scale = zcontr[1] + zcontr2[1]
    expected = map(mu) do μ
        Integrals.spherical_besselj_series(0, μ * h)[1] * expected_scale
    end
    @test res ≈ expected

    l = [0, 0]
    m = [0, 1]
    n = [0, 0]
    group_start = fill(1, 4)
    group_count = fill(2, 4)
    wide = ones(3, 3, 3)
    zcontr = ones(2, 2, 2, 2)
    zcontr2 = -ones(2, 2, 2, 2)

    res = Integrals.tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                                        hx, hy, hz, h,
                                        wide, wide, wide, wide, wide, wide,
                                        1, 1, 1, 1, zcontr, zcontr2;
                                        cutoff1=1e-14, cutoff2=1e-14)

    @test all(res .== 0.0)
end

@testset "form_factor unit scaling" begin
    mu = [0.1, 0.3]  # Å⁻¹
    l = m = n = fill(0, 2)
    group_start = [1, 2]
    group_count = [1, 1]

    geom = Integrals.IntegralGeometry([0.0 1.0; 0.0 1.0], zeros(2, 2), zeros(2, 2))

    bohr_to_ang = 0.529177210903
    hx, hy, hz, h = Integrals.pairwise_offsets(geom, 1, 1, 2, 2)
    scaled_h = bohr_to_ang * h

    dx = dy = dz = ones(1, 1, 1)
    zc = fill(0.25, 1, 1, 1, 1)

    direct = Integrals.tot_integral_k_ijkr(mu, l, m, n, group_start, group_count,
                                           bohr_to_ang * hx, bohr_to_ang * hy, bohr_to_ang * hz, scaled_h,
                                           dx, dy, dz, dx, dy, dz,
                                           1, 1, 2, 2, zc, zc;
                                           cutoff1=1e-14, cutoff2=1e-14)

    via_wrapper = Integrals.form_factor(mu, geom, l, m, n, group_start, group_count,
                                        dx, dy, dz, dx, dy, dz,
                                        1, 1, 2, 2, zc, zc;
                                        unit_scale=bohr_to_ang,
                                        cutoff1=1e-14, cutoff2=1e-14)

    @test via_wrapper ≈ direct
    @test all(via_wrapper .> 0)
end

@testset "total_scattering_calculation aggregation" begin
    q = [0.1, 0.2]
    l = m = n = fill(0, 2)
    group = fill(1, 2)  # single symmetry block with two basis functions
    geom = Integrals.IntegralGeometry([0.0 0.2; 0.2 0.0], zeros(2, 2), zeros(2, 2))
    dd = ones(1, 1, 1)

    # Use simple density tensors so every quartet contribution is identical
    zcontr = fill(0.5, 2, 2, 2, 2)
    zcontr2 = fill(0.25, 2, 2, 2, 2)

    groups = MainCalculation.group_metadata(group)
    manual = zeros(Float64, length(q))
    for gi in 1:length(groups.group_start), gj in 1:length(groups.group_start),
        gk in 1:length(groups.group_start), gr in 1:length(groups.group_start)
        hx, hy, hz, h = Integrals.pairwise_offsets(geom, gi, gj, gk, gr)
        h <= 0 && continue
        manual .+= Integrals.tot_integral_k_ijkr(q, l, m, n, groups.group_start, groups.group_count,
                                                 hx, hy, hz, h,
                                                 dd, dd, dd, dd, dd, dd, gi, gj, gk, gr, zcontr, zcontr2;
                                                 cutoff1=1e-14, cutoff2=1e-14)
    end

    threaded = MainCalculation.total_scattering_calculation(1, q, geom, l, m, n, group, dd, dd, dd;
                                                             zcontr=zcontr, zcontr2=zcontr2,
                                                             cutoff1=1e-14, cutoff2=1e-14)

    @test threaded ≈ manual
end

@testset "CI bookkeeping" begin
    confs = [
        1 0 1 0;
        0 1 1 0;
    ]

    ep2, ndiff = CIUtils.maxcoincidence(confs)
    @test ep2 == [1 1; 1 1]
    @test ndiff == [0 2; 2 0]

    civs = [0.8 0.1; 0.6 0.2]
    onerdm = CIUtils.create_onerdm(confs, civs, ndiff, ep2)
    @test size(onerdm) == (2, 2)
    @test onerdm[1, 1] ≈ 1.96 atol=1e-12
    @test onerdm[2, 2] ≈ 1.0 atol=1e-12
    @test onerdm[1, 2] == 0.0

    twordm = CIUtils.create_twordm(confs, civs, ndiff, ep2)
    @test size(twordm) == (2, 2, 2, 2)
    @test twordm[2, 1, 1, 2] ≈ 1.0 atol=1e-12
    @test twordm[1, 1, 1, 1] == 0.0
end

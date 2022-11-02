using ASEconvert
using AtomsBase
using Test
using Unitful
using LinearAlgebra
import ExtXYZ

# TODO Test reduced dimension

function test_approx_eq(s::AbstractSystem, t::AbstractSystem;
                        atol=1e-14, ignore_missing=false)
    # TODO Introduce an == / ≈ method in the AbstractSystem and use it here
    if !ignore_missing
        @test ismissing(velocity(s)) == ismissing(velocity(t))
    end

    @test maximum(norm, position(s) - position(t)) < atol * u"Å"
    @test maximum(norm, bounding_box(s) - bounding_box(t)) < atol * u"Å"

    if !ismissing(velocity(s)) && !ismissing(velocity(t))
        @test maximum(norm, velocity(s) - velocity(t)) < atol * u"Å/s"
    end

    for method in (atomic_mass, atomic_symbol, atomic_number, boundary_conditions)
        @test method(s) == method(t)
    end

    extra_properties = (:magnetic_moment, )
    for prop in extra_properties
        for (at_s, at_t) in zip(s, t)
            if hasproperty(at_s, prop) && hasproperty(at_t, prop)
                @test getproperty(at_s, prop) == getproperty(at_t, prop)
            end
        end
    end

    if s isa FlexibleSystem && t isa FlexibleSystem
        @test s.data == t.data  # Check extra data
    end
end

function make_test_system()
    # Store some random data in an AtomsBase flexible system
    n_atoms    = 5
    positions  = [randn(3) for _ = 1:n_atoms]u"Å"
    velocities = [randn(3) for _ = 1:n_atoms]u"Å/s"
    symbols    = [:H, :H, :C, :N, :He]
    magnetic_moments = [0.0, 0.0, 1.0, -1.0, 0.0]

    atoms = map(1:n_atoms) do i
        Atom(symbols[i], positions[i], velocities[i];
             magnetic_moment=magnetic_moments[i])
    end
    box = [[1.50304, 0.850344, 0.717239],
           [0.36113, 0.008144, 0.814712],
           [0.06828, 0.381122, 0.129081]]u"Å"
    bcs = [Periodic(), Periodic(), DirichletZero()]
    atomic_system(atoms, box, bcs, extra_data=42)
end

@testset "ASEconvert.jl" begin
    # Store some random data in an AtomsBase flexible system
    n_atoms    = 5
    positions  = [randn(3) for _ = 1:n_atoms]u"Å"
    velocities = [randn(3) for _ = 1:n_atoms]u"Å/s"
    symbols    = [:H, :H, :C, :N, :He]
    magnetic_moments = [0.0, 0.0, 1.0, -1.0, 0.0]

    atoms = map(1:n_atoms) do i
        Atom(symbols[i], positions[i], velocities[i];
             magnetic_moment=magnetic_moments[i])
    end
    box = [[1.50304, 0.850344, 0.717239],
           [0.36113, 0.008144, 0.814712],
           [0.06828, 0.381122, 0.129081]]u"Å"
    bcs = [Periodic(), Periodic(), DirichletZero()]
    system = atomic_system(atoms, box, bcs, extra_data=42)

    # TODO Check handling of missing for velocities

    @testset "Conversion to ASE" begin
        ase_atoms = convert_ase(system)

        D = 3
        for i = 1:D
            @test pyconvert(Vector, ase_atoms.cell[i - 1]) ≈ ustrip.(u"Å", box[i]) atol=1e-14
        end

        for (i, atom) in enumerate(ase_atoms)
            @test pyconvert(Vector, atom.position) ≈ ustrip.(u"Å", positions[i]) atol=1e-14
            @test pyconvert(String, atom.symbol) == string(symbols[i])
            @test pyconvert(Float64, atom.magmom) == magnetic_moments[i]
            @test(pyconvert(Vector, ase_atoms.get_velocities()[i - 1])
                  ≈ ustrip.(u"Å/s", velocities[i]), atol=1e-14)
        end
        @test pyconvert(Vector, ase_atoms.pbc) == [true, true, false]
        @test pyconvert(Int, ase_atoms.info["extra_data"]) == 42
    end

    @testset "Conversion AtomsBase -> ASE -> AtomsBase" begin
        newsystem = pyconvert(FlexibleSystem, convert_ase(system))
        test_approx_eq(system, newsystem)
    end

    @testset "Construction of bulk systems in ASE" begin
        bulk_Fe = pyconvert(AbstractSystem, ase.build.bulk("Fe"; cubic=true))

        a = 2.87u"Å"
        @test bounding_box(bulk_Fe) == a .* [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
        @test atomic_symbol(bulk_Fe) == [:Fe, :Fe]
        @test position(bulk_Fe) == [[0.0, 0.0, 0.0], [1.435, 1.435, 1.435]]u"Å"
        @test velocity(bulk_Fe) == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]u"Å/s"
    end

    @testset "Writing files using ASE" begin mktempdir() do d
        file = joinpath(d, "output.xyz")
        ase.io.write(file, convert_ase(system))

        newsystem = ExtXYZ.Atoms(ExtXYZ.read_frame(file))
        test_approx_eq(system, newsystem; atol=1e-7, ignore_missing=true)
    end end
end

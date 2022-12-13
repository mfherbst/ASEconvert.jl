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

    @test maximum(norm, position(s)     - position(t))     < atol * u"Å"
    @test maximum(norm, bounding_box(s) - bounding_box(t)) < atol * u"Å"

    if !ismissing(velocity(s)) && !ismissing(velocity(t))
        @test maximum(norm, velocity(s) - velocity(t)) < atol * u"Å/s"
    end

    for method in (atomic_mass, atomic_symbol, atomic_number, boundary_conditions)
        @test method(s) == method(t)
    end

    extra_atomic_props = (:charge, :covalent_radius, :vdw_radius, :magnetic_moment)
    for prop in extra_atomic_props
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

function make_test_system(; store_velocities=true)
    # Generate some random data to store in Atoms
    n_atoms    = 5
    positions  = [randn(3) for _ = 1:n_atoms]u"Å"
    velocities = [randn(3) for _ = 1:n_atoms]u"Å/s"
    symbols    = [:H, :H, :C, :N, :He]
    numbers    = [1, 1, 6, 5, 2]
    charges    = [2, 1, 3.0, -1.0, 0.0]u"e_au"
    masses     = randn(n_atoms)u"u"
    vdw_radii  = randn(n_atoms)u"Å"
    covalent_radii   = randn(n_atoms)u"Å"
    magnetic_moments = [0.0, 0.0, 1.0, -1.0, 0.0]

    atoms = map(1:n_atoms) do i
        atargs = (atomic_number=numbers[i],
                  atomic_symbol=symbols[i],
                  atomic_mass=masses[i],
                  vdw_radius=vdw_radii[i],
                  covalent_radius=covalent_radii[i],
                  charge=charges[i],
                  magnetic_moment=magnetic_moments[i])
        if store_velocities
            Atom(symbols[i], positions[i], velocities[i]; atargs...)
        else
            Atom(symbols[i], positions[i]; atargs...)
        end
    end
    box = [[1.50304, 0.850344, 0.717239],
           [0.36113, 0.008144, 0.814712],
           [0.06828, 0.381122, 0.129081]]u"Å"
    bcs = [Periodic(), Periodic(), DirichletZero()]

    atomic_system(atoms, box, bcs;
                  extra_data=42, charge=-1u"e_au", multiplicity=2)
end

@testset "ASEconvert.jl" begin
    # TODO Check handling of missing for velocities

    @testset "Conversion to ASE" begin
        system = make_test_system()
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


        # TODO Test the other properties
    end

    @testset "Conversion to ASE ignores unitful quantities" begin
        dropsystem1 = atomic_system(atoms, box, bcs, extra_data=42, len=12u"m")
        dropsystem2 = atomic_system(atoms, box, bcs, extra_data=42, mass=[2u"kg", 1u"kg"])
        for sys in (dropsystem1, dropsystem2)
            ase_atoms = @test_logs((:warn, r"Unitful quantities are not yet supported"),
                                   convert_ase(sys))
            @test pyconvert(Int, ase_atoms.info["extra_data"]) == 42
        end
    end

    @testset "Conversion AtomsBase -> ASE -> AtomsBase" begin
        system = make_test_system()
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

    @testset "Writing files using ASE" begin
        system = make_test_system()
        mktempdir() do d
            file = joinpath(d, "output.xyz")
            ase.io.write(file, convert_ase(system))

            newsystem = ExtXYZ.Atoms(ExtXYZ.read_frame(file))
            test_approx_eq(system, newsystem; atol=1e-7, ignore_missing=true)
        end
    end
end

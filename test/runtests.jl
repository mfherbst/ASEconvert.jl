using ASEconvert
using AtomsBase
using AtomsBaseTesting
using AtomsCalculators
using AtomsCalculators.Testing
using PythonCall
using Test
using Unitful
using UnitfulAtomic


@testset "ASEconvert.jl" begin
    function make_ase_system(args...; drop_atprop=Symbol[], kwargs...)
        # ASE does not support vdw_radius and covalent_radius
        dropkeys = [:covalent_radius, :vdw_radius]
        make_test_system(args...; drop_atprop=append!(drop_atprop, dropkeys), kwargs...)
    end

    # TODO Test reduced dimension
    @testset "Conversion to ASE (3D, with velocity)" begin
        (; system, cell_vectors, periodicity, atprop, sysprop) = make_test_system()
        ase_atoms = @test_logs((:warn, r"Skipping atomic property vdw_radius"),
                               (:warn, r"Skipping atomic property covalent_radius"),
                               match_mode=:any, convert_ase(system))

        D = 3
        for i = 1:D
            @test pyconvert(Vector, ase_atoms.cell[i - 1]) ≈ ustrip.(u"Å", cell_vectors[i]) atol=1e-14
        end
        @assert periodicity == (true, true, false)
        @test pyconvert(Vector, ase_atoms.pbc) == [true, true, false]

        for (i, atom) in enumerate(ase_atoms)
            @test(pyconvert(Vector, atom.position)
                  ≈ ustrip.(u"Å", atprop.position[i]), atol=1e-14)
            @test(pyconvert(Vector, ase_atoms.get_velocities()[i - 1])
                  ≈ ustrip.(sqrt(u"eV"/u"u"), atprop.velocity[i]), atol=1e-12)

            @test pyconvert(String,  atom.symbol) == string(atomic_symbol(atprop.species[i]))
            @test pyconvert(Int,     atom.number) == atomic_number(atprop.species[i])
            @test pyconvert(Float64, atom.mass)   == ustrip(u"u", atprop.mass[i])
            @test pyconvert(Float64, atom.magmom) == atprop.magnetic_moment[i]
            @test pyconvert(Float64, atom.charge) == ustrip(u"e_au", atprop.charge[i])
        end
        @test pyconvert(Int, ase_atoms.info["extra_data"])   == sysprop.extra_data
        @test pyconvert(Int, ase_atoms.info["multiplicity"]) == sysprop[:multiplicity]
        @test pyconvert(Float64, ase_atoms.info["charge"])   == ustrip(u"e_au", sysprop[:charge])
    end

    @testset "Conversion to ASE (without velocities)" begin
        system = make_ase_system(drop_atprop=[:velocity]).system
        ase_atoms = convert_ase(system)
        for (i, atom) in enumerate(ase_atoms)
            @test iszero(pyconvert(Vector, ase_atoms.get_velocities()[i - 1]))
        end
    end

    @testset "Warning about non-default neutron count" begin
        species = [ChemicalSpecies(:H), ChemicalSpecies(:He, n_neutrons=1),
                   ChemicalSpecies(:Li, n_neutrons=6),
                   ChemicalSpecies(:Be), ChemicalSpecies(:B)]
        system = make_ase_system(; extra_atprop=(; species)).system
        ase_atoms = @test_logs((:warn, r"Atom 2 has a non-default neutron count."),
                               (:warn, r"Atom 3 has a non-default neutron count."),
                               match_mode=:any, convert_ase(system))
        expected_symbols = ["H", "He", "Li", "Be", "B"]
        for (i, atom) in enumerate(ase_atoms)
            @test pyconvert(Int, atom.number)    == i
            @test pyconvert(String, atom.symbol) == expected_symbols[i]
        end
    end

    @testset "Ignoring of unitful quantities" begin
        dropsystem1 = make_ase_system(; extra_sysprop=(; len=12u"m")).system
        dropsystem2 = make_ase_system(; extra_sysprop=(; mass=[2u"kg", 1u"kg"])).system
        for sys in (dropsystem1, dropsystem2)
            ase_atoms = @test_logs((:warn, r"Unitful quantities are not yet supported"),
                                   match_mode=:any, convert_ase(sys))
            @test pyconvert(Int, ase_atoms.info["extra_data"]) == 42
        end
    end

    @testset "Conversion AtomsBase -> ASE -> AtomsBase" begin
        system = make_ase_system().system
        newsystem = pyconvert(FlexibleSystem, convert_ase(system))
        test_approx_eq(system, newsystem)
    end

    @testset "Construction of bulk systems in ASE" begin
        bulk_Fe = pyconvert(AbstractSystem, ase.build.bulk("Fe"; cubic=true))

        a = 2.87u"Å"
        @test cell_vectors(bulk_Fe) == a .* ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])
        @test atomic_symbol(bulk_Fe, :) == [:Fe, :Fe]
        @test position(bulk_Fe, :) == [[0.0, 0.0, 0.0], [1.435, 1.435, 1.435]]u"Å"
        @test velocity(bulk_Fe, :) == [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]*sqrt(u"eV"/u"u")
    end

    @testset "Writing / reading files using ASE" begin
        system = make_ase_system().system
        mktempdir() do d
            file = joinpath(d, "output.xyz")
            ase.io.write(file, convert_ase(system))
            newsystem = pyconvert(FlexibleSystem, ase.io.read(file))
            test_approx_eq(system, newsystem; rtol=1e-6)
        end
    end

    @testset "Writing files using ASE, reading using ExtXYZ" begin
        import ExtXYZ

        # Unfortunately the conventions used in ExtXYZ for storing extra
        # properties are *not* the same in ASE and ExtXYZ, therefore writing
        # in one and reading in another is not loss-free. We thus drop some
        # special things (atomic masses, magmoms, charges)
        species = ChemicalSpecies.([:H, :H, :C, :N, :He])
        masses  = [mass(sp) for sp in species]
        extra_atprop  = (; species, mass=masses)

        system = make_ase_system(; extra_atprop, drop_atprop=[:velocity]).system
        mktempdir() do d
            file = joinpath(d, "output.xyz")
            ase.io.write(file, convert_ase(system))
            newsystem = ExtXYZ.Atoms(ExtXYZ.read_frame(file))

            # Keys ExtXYZ exposes as atomkeys, but they are not standard atomkeys any more.
            extrakeys_extxyz = [:atomic_number, :atomic_symbol, :masses]
            test_approx_eq(system, newsystem; rtol=1e-6,
                           ignore_atprop=[:initial_charges, :momenta,
                                          extrakeys_extxyz...,
                                          :charge, :initial_magmoms, :magnetic_moment])
        end
    end

    @testset "Test ASEcalculator" begin
        # Setup LJ calculator in ASE
        ase_lj = pyimport("ase.calculators.lj")
        ε = ustrip(u"eV", 125.7u"K" * u"k")
        σ = ustrip(3.345u"Å")
        calc_lj = ASEcalculator(ase_lj.LennardJones(; epsilon=ε, sigma=σ))

        # Build Argon supercell
        ase_system = ase.build.bulk("Ar") * pytuple((5, 5, 5))
        system = pyconvert(AbstractSystem, ase_system)

        # Check we have not made a stupid coding mistake
        @test -0.39120116 ≈ austrip(AtomsCalculators.potential_energy(system, calc_lj))

        # Check the full AtomsCalculator test suite passes
        test_energy_forces_virial(system, calc_lj; rtol=1e-8)
    end
end

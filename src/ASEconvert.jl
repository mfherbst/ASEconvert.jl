module ASEconvert
using PythonCall
using AtomsBase
using Unitful

export convert_ase
export ase
export pyconvert, pytuple  # Reexport from PythonCall
export AbstractSystem      # Reexport from AtomsBase

"""
Global constant representing the `ase` python module available from Julia.
"""
const ase = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(ase, pyimport("ase"))
    PythonCall.pyconvert_add_rule("ase.atoms:Atoms", AbstractSystem, ase_to_system)

    # Make a bunch of submodules available
    for sub in ("ase.io", "ase.build", "ase.lattice", "ase.visualize")
        pyimport(sub)
    end
end

function ase_to_system(S::Type{<:AbstractSystem}, ase_atoms::Py)
    box = [pyconvert(Vector, ase_atoms.cell[i])u"Å" for i = 0:2]

    atnums     = pyconvert(Vector, ase_atoms.get_atomic_numbers())
    atsyms     = pyconvert(Vector, ase_atoms.symbols)
    atmasses   = pyconvert(Vector, ase_atoms.get_masses())
    positions  = pyconvert(Matrix, ase_atoms.get_positions())
    velocities = pyconvert(Matrix, ase_atoms.get_velocities())
    magmoms    = pyconvert(Vector, ase_atoms.get_initial_magnetic_moments())
    charges    = pyconvert(Vector, ase_atoms.get_initial_charges())
    ase_info   = pyconvert(Dict{String,Any}, ase_atoms.info)

    # Extract atomic info keys
    ase_info_atomic = Dict{Symbol,Any}()
    for (k, v) in ase_info
        if startswith(k, "atom__")
            key = Symbol(k[7:end])
            if key in (:vdw_radius, :covalent_radius)
                ase_info_atomic[key] = v * u"Å"
            else
                @warn "Ignoring unknown atomic key: $k"
            end
        end
    end

    atoms = map(1:length(atnums)) do i
        extra = Dict(k => v[i] for (k, v) in ase_info_atomic)
        AtomsBase.Atom(atnums[i], positions[i, :]u"Å", velocities[i, :]u"Å/s";
                       atomic_symbol=Symbol(atsyms[i]),
                       atomic_number=atnums[i],
                       atomic_mass=atmasses[i]u"u",
                       magnetic_moment=magmoms[i],
                       charge=charges[i]u"e_au",
                       extra...)
    end

    # Parse extra data in info struct
    info = Dict{Symbol, Any}()
    for (k, v) in ase_info
        if startswith(k, "atom__")
            continue
        elseif k == "charge"
            info[Symbol(k)] = v * u"e_au"
        else
            info[Symbol(k)] = v
        end
    end

    bcs = [p ? Periodic() : DirichletZero() for p in pyconvert(Vector, ase_atoms.pbc)]
    PythonCall.pyconvert_return(atomic_system(atoms, box, bcs; info...))
end

"""
    convert_ase(system::AbstractSystem)

Convert a passed `system` (which satisfies the AtomsBase.AbstractSystem interface) to an
`ase.Atoms` datastructure. Conversions to other ASE objects from equivalent Julia objects
may be added in the future.
"""
function convert_ase(system::AbstractSystem{D}) where {D}
    D != 3 && @warn "1D and 2D systems not yet fully supported."

    n_atoms = length(system)
    pbc     = map(isequal(Periodic()), boundary_conditions(system))
    symbols = map(String, atomic_symbol(system))
    numbers = atomic_number(system)
    masses  = atomic_mass(system)

    cell = zeros(3, 3)
    for (i, v) in enumerate(bounding_box(system))
        cell[i, 1:D] = ustrip.(u"Å", v)
    end

    positions = zeros(n_atoms, 3)
    for at = 1:n_atoms
        positions[at, 1:D] = ustrip.(u"Å", position(system, at))
    end

    velocities = nothing
    if !ismissing(velocity(system))
        velocities = zeros(n_atoms, 3)
        for at = 1:n_atoms
            velocities[at, 1:D] = ustrip.(u"Å/s", velocity(system, at))
        end
    end

    charges = nothing
    if any(hasproperty(atom, :charge) for atom in system)
        charges = [ustrip(u"e_au", atom.charge) for atom in system]
    end

    magmoms = nothing
    if any(hasproperty(atom, :magnetic_moment) for atom in system)
        magmoms = [atom.magnetic_moment for atom in system]
    end


    # Map extra system properties
    # TODO Probably we need some mechanism to map keys which happen to be "official"
    #      AtomsBase keys to their ASE counterparts.
    info = Dict{String, Any}()
    if system isa FlexibleSystem
        # Here we can retrieve extra data
        # TODO not a good idea to directly access the field
        # TODO Implement and make use of a property interface on the system level
        for (k, v) in system.data
            if k == :charge
                info[string(k)] = ustrip(u"e_au", v)
            elseif v isa Quantity || (v isa AbstractArray && eltype(v) <: Quantity)
                @warn("Unitful quantities are not yet supported in convert_ase. " *
                      "Ignoring key $k")
            else
                info[string(k)] = v
            end
        end
    end

    # Map extra atomic properties not available as native ASE atomic properties
    # With the upcoming interface to access generic properties one could also map
    # all extra atomic properties here. For now we just support a few (on top
    # of the :charge and :magnetic_moments natively supported by ASE).
    for key in (:covalent_radius, :vdw_radius)
        if any(hasproperty(atom, key) for atom in system)
            info["atom__" * string(key)] = ustrip(u"Å", getproperty(atom, key))
        end
    end

    ase.Atoms(; symbols, positions, numbers, masses, magmoms, charges,
              cell, pbc, velocities, info)
end

# TODO Could have a convert_ase(Vector{AbstractSystem}) to make an ASE trajectory
# TODO Could have a convert_ase(Vector{Vector{Unitful}}) to make an ASE cell
# TODO Could have a way to make an ASE calculator from an InteratomicPotential object

end

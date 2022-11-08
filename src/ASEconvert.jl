module ASEconvert
using PythonCall
using AtomsBase
using Unitful

export convert_ase
export ase
export pyconvert

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
    positions  = pyconvert(Matrix, ase_atoms.get_positions())
    velocities = pyconvert(Matrix, ase_atoms.get_velocities())
    magmoms    = pyconvert(Vector, ase_atoms.get_initial_magnetic_moments())

    atoms = [AtomsBase.Atom(atnums[i], positions[i, :]u"Å", velocities[i, :]u"Å/s";
                            magnetic_moment=magmoms[i])
             for i = 1:length(atnums)]

    # Parse extra data in info struct
    info = pyconvert(Dict{Symbol,Any}, ase_atoms.info)

    # TODO We should have some mechanism to extract "official" ASE system keys not supported
    #      in AtomsBase and put them into the system-level data.

    bcs = [p ? Periodic() : DirichletZero() for p in pyconvert(Vector, ase_atoms.pbc)]
    PythonCall.pyconvert_return(atomic_system(atoms, box, bcs; info...))
end

function convert_ase(system::AbstractSystem{D}) where {D}
    D != 3 && @warn "1D and 2D systems not yet fully supported."

    n_atoms = length(system)
    pbc     = map(isequal(Periodic()), boundary_conditions(system))
    symbols = map(String, atomic_symbol(system))

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

    # Map extra atomic properties
    # Note: This is needed since ASE only support *some* atomic properties.
    # In theory it is possible to support others using the system-level info dictionary,
    # but this gets really messy and probably unreliable, so we don't do that.
    #
    # mapkeys maps from the atomic key in the AtomsBase convention to the ase.Atoms
    # key and the default value (including units) in case an atom does not set this property.
    mapkeys = Dict(:magnetic_moment => (:magmoms, 0.0), )
    extra = Dict{Symbol,Any}()
    for (key_atbase, (key_ase, atom_default)) in pairs(mapkeys)
        if any(at -> hasproperty(at, key_atbase), system)
            data_unit = unit(atom_default)
            atpropdata = map(system) do atom
                if hasproperty(atom, key_atbase)
                    ustrip(data_unit, getproperty(atom, key_atbase))
                else
                    atom_default
                end
            end
            extra[key_ase] = atpropdata
        end
    end

    # Map extra system properties
    # TODO Probably we need some mechanism to map keys which happen to be "official"
    #      ASE keys to their ASE counterparts.
    info = Dict{String, Any}()
    if system isa FlexibleSystem
        # Here we can retrieve extra data
        # TODO not a good idea to directly access the field
        # TODO Implement and make use of a property interface on the system level
        for (k, v) in system.data
            if v isa Quantity || (v isa AbstractArray && eltype(v) <: Quantity)
                @warn("Unitful quantities are not yet supported in convert_ase. " *
                      "Ignoring key $k")
            else
                info[string(k)] = v
            end
        end
    end

    ase.Atoms(; symbols, positions, cell, pbc, velocities, info, extra...)
end

# TODO Could have a convert_ase(Vector{AbstractSystem}) to make an ASE trajectory

end

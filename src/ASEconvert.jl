module ASEconvert
using PythonCall
using AtomsBase
using Unitful
using UnitfulAtomic
import PeriodicTable

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

# For ASE units, see https://wiki.fysik.dtu.dk/ase/ase/units.html
# In particular note that uTime = u"Å" * sqrt(u"u" / u"eV") and thus
const uVelocity = sqrt(u"eV" / u"u")


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

    atoms = map(1:length(atnums)) do i
        AtomsBase.Atom(atnums[i], positions[i, :]u"Å", velocities[i, :] * uVelocity;
                       atomic_symbol=Symbol(atsyms[i]),
                       atomic_number=atnums[i],
                       atomic_mass=atmasses[i]u"u",
                       magnetic_moment=magmoms[i],
                       charge=charges[i]u"e_au")
    end

    # Parse extra data in info struct
    info = Dict{Symbol, Any}()
    for (k, v) in ase_info
        if k == "charge"
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
may be added as additional methods in the future.
"""
function convert_ase(system::AbstractSystem{D}) where {D}
    D != 3 && @warn "1D and 2D systems not yet fully supported."

    n_atoms = length(system)
    pbc     = map(isequal(Periodic()), boundary_conditions(system))
    numbers = atomic_number(system)
    masses  = ustrip.(u"u", atomic_mass(system))

    symbols_match = [
        PeriodicTable.elements[atnum].symbol == string(atomic_symbol(system, i))
        for (i, atnum) in enumerate(numbers)
    ]
    if !all(symbols_match)
        @warn("Mismatch between atomic numbers and atomic symbols, which is not " *
              "supported in ASE. Atomic numbers take preference.")
    end

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
            velocities[at, 1:D] = ustrip.(uVelocity, velocity(system, at))
        end
    end

    # We don't map any extra atom properties, which are not available in ASE as this
    # only causes a mess: ASE could do something to the atoms, but not taking
    # care of the extra properties, thus rendering the extra properties invalid
    # without the user noticing.
    charges = nothing
    magmoms = nothing
    for key in atomkeys(system)
        if key in (:position, :velocity, :atomic_symbol, :atomic_number, :atomic_mass)
            continue  # Already dealt with
        elseif key == :charge
            charges = ustrip.(u"e_au", system[:, :charge])
        elseif key == :magnetic_moment
            magmoms = system[:, :magnetic_moment]
        else
            @warn "Skipping atomic property $key, which is not supported in ASE."
        end
    end

    # Map extra system properties
    info = Dict{String, Any}()
    for (k, v) in pairs(system)
        if k in (:bounding_box, :boundary_conditions)
            continue
        elseif k in (:charge, )
            info[string(k)] = ustrip(u"e_au", v)
        elseif v isa Quantity || (v isa AbstractArray && eltype(v) <: Quantity)
            @warn("Unitful quantities are not yet supported in convert_ase. " *
                  "Ignoring key $k")
        else
            info[string(k)] = v
        end
    end

    ase.Atoms(; positions, numbers, masses, magmoms, charges,
              cell, pbc, velocities, info)
end

# TODO Could have a convert_ase(Vector{AbstractSystem}) to make an ASE trajectory
# TODO Could have a convert_ase(Vector{Vector{Unitful}}) to make an ASE cell
# TODO Could have a way to make an ASE calculator from an InteratomicPotential object

end

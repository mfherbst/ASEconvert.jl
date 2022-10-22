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
             for i in 1:length(atnums)]

    bcs = [p ? Periodic() : DirichletZero() for p in pyconvert(Vector, ase_atoms.pbc)]

    # Parse extra data in info struct
    info = pyconvert(Dict{Symbol, Any}, ase_atoms.info)
    PythonCall.pyconvert_return(atomic_system(atoms, box, bcs; info...))
end

function convert_ase(system::AbstractSystem{D}) where {D}
    n_atoms = length(system)

    cell = zeros(3, 3)
    for (i, v) in enumerate(bounding_box(system))
        cell[i, 1:D] = ustrip.(u"Å", v)
    end

    positions = zeros(n_atoms, 3)
    for at in 1:n_atoms
        positions[at, 1:D] = ustrip.(u"Å", position(system, at))
    end

    velocities = nothing
    if !ismissing(velocity(system))
        velocities = zeros(n_atoms, 3)
        for at in 1:n_atoms
            velocities[at, 1:D] = ustrip.(u"Å/s", velocity(system, at))
        end
    end

    pbc = map(isequal(Periodic()), boundary_conditions(system))
    symbols = String.(atomic_symbol(system))
    ase.Atoms(; symbols, positions, cell, pbc, velocities)
end

# TODO Should also have a more specialised version on FlexibleSystem, which can fetch more data.

# TODO Could have a convert_ase(Vector{AbstractSystem}) to make an ASE trajectory

end

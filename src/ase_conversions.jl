import PeriodicTable

# Some general notes:
#    - ASE has a cell object, but does neither make the distinction
#      between isolated and periodic cells nor do cells store
#      the boundary conditions like we do in AtomsBase. So we will
#      *not* support conversion to and from cells for now.
#    - Isotope handling is not yet fully implemented in AtomsBase,
#      so this will be ignored in this interface.
#    - Atomic numbers and atomic symbols cannot differ in ASE,
#      so this is not supported in the interface.

# TODO One could maybe also implement this as a wrapper, where data
#      is not copied, but updated in-place in the python object.

# For ASE units, see https://wiki.fysik.dtu.dk/ase/ase/units.html
# In particular note that uTime = u"Å" * sqrt(u"u" / u"eV") and thus
const uVelocity = sqrt(u"eV" / u"u")


function ase_to_system(S::Type{<:AbstractSystem}, ase_atoms::Py)
    cell_vectors = [pyconvert(Vector, ase_atoms.cell[i])u"Å" for i = 0:2]
    periodicity  = pyconvert(Vector, ase_atoms.pbc)

    # Convention in ASE isolated cells have zero (undefined) cell vectors
    if all(iszero, cell_vectors)
        if any(periodicity)
            @warn "Ignoring ASE pbc settings when ASE cell vectors are zero (undefined)"
        end
        cϵll = IsolatedCell(3, typeof(1.0u"Å"))
    else
        cϵll = PeriodicCell(; cell_vectors, periodicity)
    end

    atnums     = pyconvert(Vector, ase_atoms.get_atomic_numbers())
    atmasses   = pyconvert(Vector, ase_atoms.get_masses())
    positions  = pyconvert(Matrix, ase_atoms.get_positions())
    velocities = pyconvert(Matrix, ase_atoms.get_velocities())
    magmoms    = pyconvert(Vector, ase_atoms.get_initial_magnetic_moments())
    charges    = pyconvert(Vector, ase_atoms.get_initial_charges())
    ase_info   = pyconvert(Dict{String,Any}, ase_atoms.info)

    particles = map(1:length(atnums)) do i
        Atom(atnums[i],
             positions[i, :]u"Å",
             velocities[i, :] * uVelocity;
             species=ChemicalSpecies(atnums[i]),
             mass=atmasses[i]u"u",
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

    PythonCall.pyconvert_return(FlexibleSystem(particles, cϵll; info...))
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

    ase_cell = zeros(3, 3)
    if cell(system) isa IsolatedCell
        pbc = [false, false, false]
    else
        pbc = periodicity(system)
        for (i, v) in enumerate(cell_vectors(system))
            ase_cell[i, 1:D] = ustrip.(u"Å", v)
        end
    end

    for at = 1:n_atoms
        spec = species(system, at)
        if spec.n_neutrons ≥ 0
            @warn("Atom $at has a non-default neutron count. This information " *
                  "is dropped when converting to ASE.")
        end
    end

    numbers = [atomic_number(species(system, at)) for at = 1:n_atoms]
    masses  = [ustrip(u"u", mass(system, at))     for at = 1:n_atoms]

    positions = zeros(n_atoms, 3)
    for at = 1:n_atoms
        positions[at, 1:D] = ustrip.(u"Å", position(system, at))
    end

    velocities = nothing
    if n_atoms > 0 && !ismissing(velocity(system, 1))
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
        if key in (:position, :velocity, :species, :mass)
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
        if k in (:cell_vectors, :periodicity)
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
              cell=ase_cell, pbc, velocities, info)
end

# TODO Could have a convert_ase(Vector{AbstractSystem}) to make an ASE trajectory

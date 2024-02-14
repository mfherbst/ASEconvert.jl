import AtomsCalculators

export ASEcalculator

"""
    ASEcalculator{T}

This structure wraps python ASE calculator to AtomsCalculators
interface compatible structure.

# Example
```julia
using AtomsCalculators
using ASEconvert
using PythonCall

fname = "path to eam potential file"
EAM = pyimport("ase.calculators.eam")
eam_cal = ASEconvert.ASEcalculator(EAM.EAM(potential=fname))

atoms_ase = ase.build.bulk("Ni") * pytuple((4, 3, 2))
atoms_ab = pyconvert(AbstractSystem, atoms_ase)

AtomsCalculators.potential_energy(atoms_ab, eam_cal)
AtomsCalculators.forces(atoms_ab, eam_cal)
```

"""
mutable struct ASEcalculator{T}
    ase_python_calculator::T
end


AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
    system,
    calc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(system)
    e = calc.ase_python_calculator.get_potential_energy(ase_system)
    return pyconvert(Float64, e) * u"eV"
end


AtomsCalculators.@generate_interface function AtomsCalculators.forces(
    system,
    calc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(system)
    f = calc.ase_python_calculator.get_forces(ase_system)
    # We need to convert from Python row major to Julia collumn major
    # and from there to Vector with SVector{3,Float64} element.
    # We'll do reallocation in the end to make sure that allignement is correct
    tmp =  pyconvert(Array, f)
    FT = AtomsCalculators.promote_force_type(system, calc)
    tmp2 = reinterpret(FT, tmp')
    return Vector(vec(tmp2))
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(
    system,
    calc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(system)
    tmp = calc.ase_python_calculator.get_stress(ase_system)
    cons = ase.constraints
    stress = cons.voigt_6_to_full_3x3_stress(tmp) * ( - ase_system.get_volume() )
    return pyconvert(Array, stress) * u"eV"
end
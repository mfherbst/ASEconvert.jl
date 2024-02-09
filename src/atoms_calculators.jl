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
    ab_system,
    ASEcalc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(ab_system)
    e = ASEcalc.ase_python_calculator.get_potential_energy(ase_system)
    return pyconvert(Float64, e) * u"eV"
end


AtomsCalculators.@generate_interface function AtomsCalculators.forces(
    ab_system,
    ASEcalc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(ab_system)
    f = ASEcalc.ase_python_calculator.get_forces(ase_system)
    tmp =  pyconvert(Array, f)
    tmp2 = reinterpret(AtomsCalculators.promote_force_type(ab_system,  ASEcalc), tmp') # |> Vector
    # This is a total hack, but better this than have complications later
    return vec(tmp2) |> Vector
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(
    ab_system,
    ASEcalc::ASEcalculator;
    kwargs...
)
    ase_system = convert_ase(ab_system)
    tmp = ASEcalc.ase_python_calculator.get_stress(ase_system)
    cons = ase.constraints
    stress = cons.voigt_6_to_full_3x3_stress(tmp) * ( - ase_system.get_volume() )
    return pyconvert(Array, stress) * u"eV"
end
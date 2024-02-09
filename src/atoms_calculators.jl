import AtomsCalculators

export ASEcalculator

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
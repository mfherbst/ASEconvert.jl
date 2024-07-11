import AtomsCalculators
import AtomsCalculators: @generate_interface, calculate, Energy, Forces, Virial


"""
    ASEcalculator

This structure wraps python ASE calculator to AtomsCalculators
interface compatible structure.

# Example
```julia
using AtomsCalculators
using ASEconvert
using PythonCall

potential = "path to eam potential file"
EAM = pyimport("ase.calculators.eam")
eam_cal = ASEconvert.ASEcalculator(EAM.EAM(potential))

atoms_ase = ase.build.bulk("Ni") * pytuple((4, 3, 2))
atoms_ab = pyconvert(AbstractSystem, atoms_ase)

AtomsCalculators.potential_energy(atoms_ab, eam_cal)
AtomsCalculators.forces(atoms_ab, eam_cal)
```

"""
mutable struct ASEcalculator
    calculator::PythonCall.Py
end

@generate_interface function calculate(::Energy, system, calc::ASEcalculator,
                                       ps=nothing, st=nothing; kwargs...)
    e = calc.calculator.get_potential_energy(convert_ase(system))
    pyconvert(Float64, e) * u"eV"
end

@generate_interface function calculate(::Forces, system, calc::ASEcalculator,
                                       ps=nothing, st=nothing; kwargs...)
    f = calc.calculator.get_forces(convert_ase(system))
    # We need to convert from Python row major to Julia collumn major
    # and from there to Vector with SVector{3,Float64} element.
    # We'll do reallocation in the end to make sure that allignement is correct
    tmp =  pyconvert(Array, f)
    FT = AtomsCalculators.promote_force_type(system, calc)
    tmp2 = reinterpret(FT, tmp')
    Vector(vec(tmp2))
end

@generate_interface function calculate(::Virial, system, calc::ASEcalculator,
                                       ps=nothing, st=nothing; kwargs...)
    ase_system = convert_ase(system)
    tmp = calc.calculator.get_stress(ase_system)
    Ω = ase_system.get_volume()
    stress = -Ω * ase.constraints.voigt_6_to_full_3x3_stress(tmp)
    pyconvert(Array, stress) * u"eV"
end

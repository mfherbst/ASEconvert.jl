import AtomsCalculators
import AtomsCalculators: @generate_interface, calculate, Energy, Forces, Virial


"""
    ASEcalculator

This structure wraps python ASE calculator to AtomsCalculators
interface compatible structure.

# Example
```julia
using ASEconvert
using AtomsBuilder
using AtomsCalculators
using PythonCall

ase_emt = pyimport("ase.calculators.emt")
calc_emt = ASEcalculator(ase_emt.EMT())

system = bulk(:Cu) * (4, 3, 2)
AtomsCalculators.potential_energy(system, calc_emt)
AtomsCalculators.forces(system, calc_emt)
AtomsCalculators.virial(system, calc_emt)
```
"""
mutable struct ASEcalculator
    calculator::PythonCall.Py
end

# TODO The ASE calculator could change in-place, so it therefore itself
#      is basically the state. The state handling should overall therefore
#      be done more thoroughly here.

AtomsCalculators.energy_unit(::ASEcalculator) = u"eV"
AtomsCalculators.length_unit(::ASEcalculator) = u"Å"

@generate_interface function AtomsCalculators.calculate(::AtomsCalculators.Energy,
        system, calc::ASEcalculator, parameters=nothing, state=nothing; kwargs...)
    edata = calc.calculator.get_potential_energy(convert_ase(system))
    (; energy=pyconvert(Float64, edata)u"eV", state)
end

@generate_interface function AtomsCalculators.calculate(::AtomsCalculators.Forces,
        system, calc::ASEcalculator, parameters=nothing, state=nothing; kwargs...)
    fdata = calc.calculator.get_forces(convert_ase(system))
    # We need to convert from Python row major to Julia column major
    # and from there to Vector with SVector{3,Float64} element.
    # We do reallocation in the end to ensure we have the correct memory allignement
    tmp  = pyconvert(Array, fdata)
    FT   = AtomsCalculators.promote_force_type(system, calc)
    tmp2 = reinterpret(FT, tmp')
    (; forces=Vector(vec(tmp2)), state)
end

@generate_interface function AtomsCalculators.calculate(::AtomsCalculators.Virial,
        system, calc::ASEcalculator, parameters=nothing, state=nothing; kwargs...)
    ase_system = convert_ase(system)
    sdata = calc.calculator.get_stress(ase_system)
    Ω = ase_system.get_volume()
    stress = -Ω * ase.constraints.voigt_6_to_full_3x3_stress(sdata)
    (; virial=pyconvert(Array, stress) * u"eV", state)
end

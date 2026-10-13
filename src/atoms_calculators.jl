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
struct ASEcalculator
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
    # ASE returns an (n_atoms, 3) array in eV/Å. Transpose to a contiguous 3×n_atoms
    # matrix, so that each column can be reinterpreted as one force vector of type FT
    # (an SVector{3} carrying the eV/Å unit), then copy into a plain Vector{FT}.
    FT = AtomsCalculators.promote_force_type(system, calc)
    forces_3xn = permutedims(pyconvert(Matrix{Float64}, fdata))
    (; forces=Vector(reinterpret(reshape, FT, forces_3xn)), state)
end

@generate_interface function AtomsCalculators.calculate(::AtomsCalculators.Virial,
        system, calc::ASEcalculator, parameters=nothing, state=nothing; kwargs...)
    ase_system = convert_ase(system)
    sdata = calc.calculator.get_stress(ase_system)
    Ω = ase_system.get_volume()
    stress = -Ω * ase.stress.voigt_6_to_full_3x3_stress(sdata)
    (; virial=pyconvert(Array, stress) * u"eV", state)
end

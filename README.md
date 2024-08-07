# ASEconvert
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mfherbst.github.io/ASEconvert.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mfherbst.github.io/ASEconvert.jl/dev/)
[![Build Status](https://github.com/mfherbst/ASEconvert.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mfherbst/ASEconvert.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/mfherbst/ASEconvert.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mfherbst/ASEconvert.jl)

Light-weight module to install ASE
and provide routines for converting between the Atoms datastructure
from [ASE](https://wiki.fysik.dtu.dk/ase/index.html)
and atomistic data provided in
the [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) ecosystem.
For both the package makes use of the
[PythonCall](https://github.com/cjdoris/PythonCall.jl/).


This can be used for example as follows
```julia
using ASEconvert

# Make a silicon supercell using ASE
atoms_ase = ase.build.bulk("Si") * pytuple((4, 1, 1))

# Convert to an AtomsBase-compatible structure
atoms_ab = pyconvert(AbstractSystem, atoms_ase)

# Convert back to ASE and create a vacancy
newatoms_ase = convert_ase(atoms_ab)
newatoms_ase.pop(4)
```

### AtomsCalculators interface

You can use ASE calculators in julia, by wrapping them to a `ASEcalculator` structure.
Here is a brief example:

```julia
using AtomsCalculators
using AtomsBuilder
using ASEconvert
using PythonCall

# Setup calculator in ASE
potential = "path to eam potential file"
ase_calc = pyimport("ase.calculators.eam").EAM(potential)

# Convert into AtomsCalculator-compatible calculator
calc = ASEcalculator(ase_calc)

# Use it to compute a Nickel supercell
system = bulk(:Ni) * (4, 3, 2)
AtomsCalculators.potential_energy(system, calc)
AtomsCalculators.forces(system, calc)
AtomsCalculators.virial(system, calc)
```

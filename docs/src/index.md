```@meta
CurrentModule = ASEconvert
```

# ASEconvert

Light-weight module to install ASE
and provide routines for converting between the Atoms datastructure
from [ASE](https://wiki.fysik.dtu.dk/ase/index.html)
and atomistic data provided in
the [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) ecosystem.

## Automatic ASE installation
Using the mechanism provided by [PythonCall](https://github.com/cjdoris/PythonCall.jl)
and [CondaPkg](https://github.com/cjdoris/CondaPkg.jl) ASEconvert will automatically
take care of installing ASE and exporting a useful subset of its modules under the
`ase` variable. For example one may easily create bulk systems

```@example
using ASEconvert
ase.build.bulk("Mg")
```

or surfaces

```@example
using ASEconvert
ase.build.surface(ase.build.bulk("Mg"), (1, 1, 0), 4, 0, periodic=true)
```

## Conversion from ASE to AtomsBase

```julia
using ASEconvert
using AtomsBase
using DFTK

# Construct bulk magnesium using ASE and convert to atomsbase
mg_ase = ase.build.bulk("Mg")
mg_atb = pyconvert(AbstractSystem, mg_ase)

# Attach pseudopotentials, construct LDA DFT model and solve for DFT ground state
system = attach_psp(mg_atb; family="hgh", functional="lda")
model  = model_LDA(system)
basis  = PlaneWaveBasis(model; Ecut=10, kgrid=(1, 1, 1))
scfres = self_consistent_field(basis)
```

## Conversion from AtomsBase to ASE

```@example
using ASEconvert
using ExtXYZ

# Read an extxyz file using ExtXYZ.jl.
frame_xyz = Atoms(read_frame("Mn3Si.extxyz"))

# The returned ExtXYZ.Atoms object is an AtomsBase.AbstractSystem,
# thus we can directly convert it to an ase.Atoms using `convert_ase`
# and write it again as an ASE json file
ase.io.write("out.json", convert_ase(frame_xyz))
```

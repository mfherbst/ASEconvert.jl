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

```@example dftk
using ASEconvert
using DFTK

# Construct bulk magnesium using ASE and convert to atomsbase
mg_ase = ase.build.bulk("Mg")
mg_atb = pyconvert(AbstractSystem, mg_ase)
```

```@example dftk
# Attach pseudopotentials, construct LDA DFT model and solve for DFT ground state
system = attach_psp(mg_atb; Mg="hgh/lda/mg-q2")
model  = model_LDA(system; temperature=1e-3, smearing=Smearing.MarzariVanderbilt())
basis  = PlaneWaveBasis(model; Ecut=20, kgrid=(4, 4, 4))
scfres = self_consistent_field(basis)

scfres.energies
```

## Conversion from AtomsBase to ASE

```@example extxyz
using ASEconvert
using AtomsIO

# Read an extxyz file using AtomsIO.jl.
system = load_system("Mn3Si.extxyz")
```

This example uses [AtomsIO](https://github.com/mfherbst/AtomsIO.jl)
to read the extended XYZ file file `Mn3Si.extxyz`. The data is returned
as a subtype of `AtomsBase.AbstractSystem`
(in this case an `ExtXYZ.Atoms` from [ExtXYZ](https://github.com/libAtoms/ExtXYZ.jl)).
We can thus directly convert this system to an `ase.Atoms` using [`convert_ase`](@ref)
and write it again as an ASE json file

```@example extxyz
ase.io.write("out.json", convert_ase(system));
```

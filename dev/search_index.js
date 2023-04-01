var documenterSearchIndex = {"docs":
[{"location":"apireference/","page":"API reference","title":"API reference","text":"CurrentModule = ASEconvert","category":"page"},{"location":"apireference/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"apireference/","page":"API reference","title":"API reference","text":"Modules = [ASEconvert]","category":"page"},{"location":"apireference/#ASEconvert.ase","page":"API reference","title":"ASEconvert.ase","text":"Global constant representing the ase python module available from Julia.\n\n\n\n\n\n","category":"constant"},{"location":"apireference/#ASEconvert.convert_ase-Union{Tuple{AbstractSystem{D}}, Tuple{D}} where D","page":"API reference","title":"ASEconvert.convert_ase","text":"convert_ase(system::AbstractSystem)\n\nConvert a passed system (which satisfies the AtomsBase.AbstractSystem interface) to an ase.Atoms datastructure. Conversions to other ASE objects from equivalent Julia objects may be added as additional methods in the future.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ASEconvert","category":"page"},{"location":"#ASEconvert","page":"Home","title":"ASEconvert","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Light-weight module to install ASE and provide routines for converting between the Atoms datastructure from ASE and atomistic data provided in the AtomsBase ecosystem.","category":"page"},{"location":"#Automatic-ASE-installation","page":"Home","title":"Automatic ASE installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Using the mechanism provided by PythonCall and CondaPkg ASEconvert will automatically take care of installing ASE and exporting a useful subset of its modules under the ase variable. For example one may easily create bulk systems","category":"page"},{"location":"","page":"Home","title":"Home","text":"using ASEconvert\nase.build.bulk(\"Mg\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"or surfaces","category":"page"},{"location":"","page":"Home","title":"Home","text":"using ASEconvert\nase.build.surface(ase.build.bulk(\"Mg\"), (1, 1, 0), 4, 0, periodic=true)","category":"page"},{"location":"#Conversion-from-ASE-to-AtomsBase","page":"Home","title":"Conversion from ASE to AtomsBase","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using ASEconvert\nusing DFTK\n\n# Construct bulk magnesium using ASE and convert to atomsbase\nmg_ase = ase.build.bulk(\"Mg\")\nmg_atb = pyconvert(AbstractSystem, mg_ase)","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Attach pseudopotentials, construct LDA DFT model and solve for DFT ground state\nsystem = attach_psp(mg_atb; Mg=\"hgh/lda/mg-q2\")\nmodel  = model_LDA(system; temperature=1e-3, smearing=Smearing.MarzariVanderbilt())\nbasis  = PlaneWaveBasis(model; Ecut=20, kgrid=(4, 4, 4))\nscfres = self_consistent_field(basis)\n\nscfres.energies","category":"page"},{"location":"#Conversion-from-AtomsBase-to-ASE","page":"Home","title":"Conversion from AtomsBase to ASE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using ASEconvert\nusing AtomsIO\n\n# Read an extxyz file using AtomsIO.jl.\nsystem = load_system(\"Mn3Si.extxyz\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"This example uses AtomsIO to read the extended XYZ file file Mn3Si.extxyz. The data is returned as a subtype of AtomsBase.AbstractSystem (in this case an ExtXYZ.Atoms from ExtXYZ). We can thus directly convert this system to an ase.Atoms using convert_ase and write it again as an ASE json file","category":"page"},{"location":"","page":"Home","title":"Home","text":"ase.io.write(\"out.json\", convert_ase(system));","category":"page"}]
}
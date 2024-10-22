module ASEconvert
using PythonCall
using AtomsBase
using Unitful
using UnitfulAtomic

export ase
export pyconvert       # Reexport from PythonCall
export AbstractSystem  # Reexport from AtomsBase

export ASEcalculator
export convert_ase

"""
Global constant representing the `ase` python module available from Julia.
"""
const ase = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(ase, pyimport("ase"))
    PythonCall.pyconvert_add_rule("ase.atoms:Atoms", AbstractSystem, ase_to_system)

    # Make a bunch of submodules available
    for sub in ("ase.io", "ase.build", "ase.lattice", "ase.visualize")
        pyimport(sub)
    end
end


include("ase_conversions.jl")
include("atoms_calculators.jl")


end

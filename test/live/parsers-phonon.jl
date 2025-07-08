include("./../../src/parsers.jl")
include("./../../src/phonon.jl")

qedata = Parsers.dfpt_qeOutputData("../../examples/espresso.ifc2")
ebdata = Parsers.ebInputData("../../examples/input.nml")

ifc2 = qedata.properties["ifc2"]

# Build some example data to test
super_mult = [4, 4, 4]
unit_mult = [2, 2, 2]
lattvecs = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
basis = [0.5 0.5 0.5; 0.0 0.0 0.0]
numbasisatoms = size(basis)[1]

basisconnectors = Phonon.build_basisconnectors(numbasisatoms, basis)

unit_points = Phonon.build_unitcell_points(unit_mult, super_mult, lattvecs)

ultracons = Phonon.get_origin_ultraconnectors(numbasisatoms, unit_points, basisconnectors)

for i in 1:size(ultracons)[1]
    println(ultracons[i, 1], " ", ultracons[i, 2], " ", ultracons[i, 3])
end

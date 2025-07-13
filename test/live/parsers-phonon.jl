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

super_points, super_sqmods = Phonon.build_unitcell_points(unit_mult, super_mult, lattvecs)
unit_points = Phonon.build_unitcell_points(unit_mult, super_mult, lattvecs)

shiftercons = Phonon.get_shiftercons(numbasisatoms, unit_points, basisconnectors)

# for i in 1:size(shiftercons)[1]
#     println(shiftercons[i, 1], " ", shiftercons[i, 2], " ", shiftercons[i, 3])
# end

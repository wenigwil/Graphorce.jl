include("./../../src/parsers.jl")
include("./../../src/phonon.jl")

qedata = Parsers.dfpt_qeOutputData("../../examples/espresso.ifc2")
ebdata = Parsers.ebInputData("../../examples/input.nml")

ifc2 = qedata.properties["ifc2"]

# Build some example data to test
super_mult = [4, 4, 4]
unit_mult = [2, 2, 2]
lattvecs = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
basis = [0.0 0.0 0.0; 0.6 0.6 0.6]
numbasisatoms = size(basis)[1]

basisconnectors = Phonon.build_basisconnectors(numbasisatoms, basis)

super_points, super_sqmods = Phonon.build_supercell_points(unit_mult, super_mult, lattvecs)
unit_points = Phonon.build_unitcell_points(unit_mult, super_mult, lattvecs)

shiftercons = Phonon.get_shiftercons(numbasisatoms, unit_points, basisconnectors)

weight_map = Phonon.get_weight_map(shiftercons, super_points, super_sqmods)
demux_unit_addr = Phonon.get_demux_unit_addr(super_mult, unit_mult)

# Printing non-zero 
for iat in 1:numbasisatoms
    for jat in 1:numbasisatoms
        for unit_addr in 1:size(unit_points)[1]
            if !iszero(weight_map[unit_addr, jat, iat])
                println(
                    "Non-Zero Element at index ",
                    (unit_addr, jat, iat),
                    ". With value ",
                    1 / weight_map[unit_addr, jat, iat],
                    ". Demuxed Index is ",
                    (demux_unit_addr[unit_addr, :], jat, iat),
                )
            end
        end
    end
end

scell = 6
for i in (-2 * scell):(2 * scell)
    x = mod(i + 1, scell)
    if x == 0
        x += scell
    end
    println(x, "\t", i)
end

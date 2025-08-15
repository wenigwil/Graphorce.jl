println("TEST: including master source")
include("../src/Graphorce.jl")

println("TEST: using .Graphorce")
using .Graphorce

println("TEST: Parsing example-files");
ebdata = ebInputData("examples/input.nml")
qedata = qeIfc2Output("examples/espresso.ifc2");

println("TEST: Constructing the deconvolution data")
deconv = DeconvData(ebdata);

# Quick and dirty function for reading a highsymmetry path
function read_highsympath(file::AbstractString)
    sympathfile = open(file)

    numlines = parse(Int64, readline(sympathfile))

    sympath = Matrix{Float64}(undef, (numlines, 3))
    for line in 1:numlines
        test = parse.(Float64, split(readline(sympathfile)))
        sympath[line, :] = test
    end

    close(sympathfile)
    return sympath
end

# q-points in crystal coordinates
qpoints_cryst = read_highsympath("examples/fcc_highsympath.txt");

# lattice vectors
lattvecs = ebdata.crystal_info["lattvecs"];

# interatomic force constants
ifc2 = qedata.properties["ifc2"];

# translate between the atomic index and species that sits there
basisatoms2species = ebdata.crystal_info["atomtypes"];

# masses of each species
species2masses = ebdata.crystal_info["masses"];

# println(qpoints_cryst[3, :])

dynmat = build_dynamical_matrix(
    deconv,
    qpoints_cryst,
    lattvecs,
    ifc2,
    basisatoms2species,
    species2masses,
)

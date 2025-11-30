# Include like this to test from root of project
include("../../src/Phunky.jl")
using .Phunky

phonons = Phonons("examples/input.nml", "examples/force.fc3")

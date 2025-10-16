# Include like this to test from root of project
include("../../src/Graphorce.jl")
using .Graphorce

phonons = Phonons("examples/input.nml", "examples/force.fc3")

import Plots

include("../src/Graphorce.jl")

using .Graphorce

harmonic_phonons = LatticeVibrations(
    "examples/input.nml",
    "examples/espresso.ifc2",
    "examples/fcc_highsympath.txt",
)

fullq_freqs = harmonic_phonons.fullq_freqs
sympath = harmonic_phonons.sympath

numbands = size(fullq_freqs, 2)

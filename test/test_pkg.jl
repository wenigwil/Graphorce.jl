import Plots

include("../src/Graphorce.jl")

using .Graphorce

highsympoints = begin
    [
        0.00 0.00 0.00 # Gamma
        0.00 0.50 0.50 # X
        0.25 0.75 0.50 # W
        0.375 0.25 0.375 # K 
        0.0 0.0 0.0 # Gamma
        0.5 0.5 0.5 # L
    ]
end

harmonic_phonons = LatticeVibrations(
    "examples/input.nml",
    "examples/espresso.ifc2",
    "examples/fcc_highsympath.txt",
)

fullq_freqs = harmonic_phonons.fullq_freqs
sympath = harmonic_phonons.sympath

numbands = size(fullq_freqs, 2)

import Plots
import LinearAlgebra as LinAlg
include("../src/Graphorce.jl")

using .Graphorce

# Seekpath cF1 high symmetry points
seek_path_basis = [
    -1.0 1.0 1.0
    1.0 -1.0 1.0
    1.0 1.0 -1.0
]

point_labels = ["Gamma", "K", "L", "U", "W", "X", "W2"]
seek_path_points = [
    0.0 0.0 0.0
    0.375 0.375 0.75
    0.5 0.5 0.5
    0.625 0.25 0.625
    0.5 0.25 0.75
    0.5 0.0 0.5
    0.75 0.25 0.5
]

seek_path_1 = seek_path_points[[1, 6, 4], :]
seek_path_2 = seek_path_points[[2, 1, 3, 5, 6, 7], :]

# println("SeeK Path Points")
# for i in axes(seek_path_points, 1)
#     println(point_labels[i], '\t', seek_path_points[i, :])
# end

harmonic_phonons = LatticeVibrations(
    "examples/input.nml",
    "examples/espresso.ifc2",
    seek_path_1,
    seek_path_2,
);

fullq_freqs = harmonic_phonons.fullq_freqs;
sympath = harmonic_phonons.sympath;

# Plots.plot(sympath, fullq_freqs; grid = false)
# Plots.savefig("phonon-disp.pdf")

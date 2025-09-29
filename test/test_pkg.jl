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

point_labels = ["Γ", "K", "L", "U", "W", "X", "W2"]
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
    numpoints_per_section = 50,
);

fullq_freqs = harmonic_phonons.fullq_freqs;
qpoints_cryst = harmonic_phonons.qpoints_cryst;
stitches = harmonic_phonons.stitches;

distances, xticks_labels, xticks_pos =
    Graphorce.path_to_xaxis(qpoints_cryst, stitches, seek_path_points, point_labels)

p = Plots.plot(
    distances,
    fullq_freqs;
    legend_position = false,
    grid = false,
    xlims = (minimum(distances), maximum(distances)),
    ylims = (minimum(fullq_freqs), 1.05 * maximum(fullq_freqs)),
    framestyle = :semi,
    tickdirection = :out,
    color = :blue,
    lw = 0.7,
    xticks = (xticks_pos, xticks_labels),
    ylabel = "Frequency ω [THz]",
)
Plots.vline!(xticks_pos[2:(end - 1)], lc = :black, lw = 0.9)

# Plots.savefig("phonon-disp.pdf")

import Plots;
Plots.pgfplotsx();
import LinearAlgebra as LinAlg
include("../src/Graphorce.jl")

using .Graphorce

# Seekpath cF1 high symmetry points
point_labels = ["Γ", "K", "L", "U", "W", "X", raw"$\mathrm{W}_2$"]
seek_path_points = [
    0.0 0.0 0.0
    0.375 0.375 0.75
    0.5 0.5 0.5
    0.625 0.25 0.625
    0.5 0.25 0.75
    0.5 0.0 0.5
    0.75 0.25 0.5
]
# Build a "walkable" path from the seekpath symmtry points
seek_path_1 = seek_path_points[[1, 6, 4], :]
seek_path_2 = seek_path_points[[2, 1, 3, 5, 6, 7], :]

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

extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

p = Plots.plot(
    distances,
    fullq_freqs;
    legend_position = false,
    grid = false,
    xlims = (minimum(distances), maximum(distances)),
    ylims = (minimum(fullq_freqs), 1.05 * maximum(fullq_freqs)),
    framestyle = :box,
    tickdirection = :out,
    color = :blue,
    lw = 0.7,
    xticks = (xticks_pos, xticks_labels),
    ylabel = "Frequency ω [THz]",
    tex_output_standalone = true,
    xtickfontsize = 12,
    ytickfontsize = 12,
    extra_kwargs = Dict(:subplot => extra_dict),
)

Plots.vline!(xticks_pos[2:(end - 1)], lc = :black, lw = 0.9)
Plots.savefig("phonon-disp.pdf")

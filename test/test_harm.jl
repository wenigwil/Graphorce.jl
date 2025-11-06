import Plots;
using BenchmarkTools;
Plots.pgfplotsx();
using Graphorce

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
# Build a "walkable" path from the seekpath symmetry points
seek_path_1 = seek_path_points[[1, 6, 4], :]
seek_path_2 = seek_path_points[[2, 1, 3, 5, 6, 7], :]

# Read-in the main description of the system
@info "Reading input.nml..."
ebinput = ebInputData("examples/input.nml")

# Construct a ifc2-deconvolution based on the system description
@info "Calculating deconvolution..."
deconvolution = DeconvData(ebinput)

# Read-in the second order force-constants
@info "Reading ifc2..."
qeinput = qeIfc2Output("examples/espresso.ifc2")

# Construct a qpoint list from a path and supply it with things to make it plottable
@info "Building path..."

sympath = Sympath(
    seek_path_points,
    point_labels,
    seek_path_1,
    seek_path_2;
    numpoints_per_section = 2000,
)

print("Generated a path of length", size(sympath.qpoints), "\n")

# Main computation of the harmonic phonon properties
@info "Calculating Lattice Vibrations..."
lattvibr = @btime LatticeVibrations(ebinput, qeinput, deconvolution, sympath.qpoints);

# extra_dict =
#     Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")
#
# # Making the names alittle shorter
# distances = sympath.distances
# freqs = lattvibr.fullq_freqs
#
# p = Plots.plot(
#     distances,
#     freqs;
#     legend_position = false,
#     grid = false,
#     xlims = (minimum(distances), maximum(distances)),
#     ylims = (minimum(freqs), 1.05 * maximum(freqs)),
#     framestyle = :box,
#     tickdirection = :out,
#     color = :blue,
#     lw = 0.7,
#     xticks = (sympath.xticks_pos, sympath.xticks_labels),
#     ylabel = "Frequency ω [THz]",
#     tex_output_standalone = true,
#     xtickfontsize = 12,
#     ytickfontsize = 12,
#     extra_kwargs = Dict(:subplot => extra_dict),
# )
#
# Plots.vline!(sympath.xticks_pos[2:(end - 1)], lc = :black, lw = 0.9)
# Plots.savefig("phonon-disp.pdf")

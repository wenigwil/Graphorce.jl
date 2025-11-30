using Phunky

# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")

# Build a symmetry path along which the phonon lifetime will be computed
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
# Walk the walk. Only small for now
route = seek_path_points[[2, 1], :]
sympath = Sympath(seek_path_points, point_labels, route; numpoints_per_section = 50)

q1_cryst = sympath.qpoints
brillouin_sampling = (5, 5, 5)

states =
    HarmonicStatesData(ebdata, sodata, deconvolution, q1_cryst; brillouin_sampling)

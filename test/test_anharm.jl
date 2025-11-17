using Graphorce

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
sympath = Sympath(seek_path_points, point_labels, route; numpoints_per_section = 3)

# Read a phonopy ifc3 file
todata = Ifc3Output("examples/force.fc3")

cont_freqs = collect(range(12.5, 15.0, 20))
kbT = 25e-3 # 300K * k_B in eV
smearing = 0.06
phonons = Phonons(
    ebdata,
    deconvolution,
    sodata,
    todata,
    sympath.qpoints,
    cont_freqs,
    kbT,
    smearing;
    brillouin_sampling = (3, 3, 3),
)

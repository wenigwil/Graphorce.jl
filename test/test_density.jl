using Graphorce
using Plots

# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")

dense = DensityOfStates(ebdata, sodata, deconvolution, 2000, (80, 80, 80), 0.7e-6)

density = dense.density
cont_energies = dense.cont_energies

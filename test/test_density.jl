using Graphorce

# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")

dense = DensityOfStates(ebdata, sodata, deconvolution, 300, (3, 3, 3), 0.01)

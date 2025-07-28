include("../src/Graphorce.jl")

using .Graphorce

ebdata = ebInputData("examples/input.nml")
qedata = qeIfc2Output("examples/espresso.ifc2")

deconv = Graphorce.DeconvUtils(ebdata)

println(typeof(deconv.ifc2_2nd_atom_unitaddr))

module Phunky

import Logging
import LinearAlgebra as LinAlg

include("constants.jl")
include("misc.jl")
export Sympath
include("crystal.jl")
include("parsers/elphbolt-nml.jl")
include("parsers/qe-ifc2-out.jl")
include("parsers/thirdorder-ifc3-out.jl")
export ebInputData
export qeIfc2Output
export Ifc3Output

include("deconvolute.jl")
export DeconvData

include("harmonic.jl")
export LatticeVibrations
export DensityOfStates

include("state.jl")
include("anharmonic.jl")
export Phonons

include("tests.jl")

end # module Phunky

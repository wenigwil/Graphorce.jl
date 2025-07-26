import Logging

include("../../src/parsers.jl")

data = Parsers.ebInputData("./../../examples/input.nml")

println("#########eb_input_test.jl#########")
println("#########eb_input_test.jl#########")
println("#########eb_input_test.jl#########")
println("#########eb_input_test.jl#########")
println("")
for (k, v) in data.crystal_info
    println(k, " ", v, " TYPE = ", typeof(v))
end
println("")
for (k, v) in data.numerics
    println(k, " ", v, " TYPE = ", typeof(v))
end

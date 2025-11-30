include("../../src/parsers.jl")

dfpt_data = Parsers.dfpt_qeOutputData("../../examples/espresso.ifc2")
ebinput_data = Parsers.ebInputData("../../examples/input.nml")

ci_comp = ["masses", "epsilon", "born", "polar", "basis"]
allocs_comp = ["numelements", "numatoms"]

for attr in allocs_comp
    println(attr)
    println("ebInputData: ", ebinput_data.allocations[attr])
    println("qeOutputData: ", dfpt_data.properties[attr], "\n")
end

for attr in ci_comp
    println(attr)
    println("ebInputData: ", ebinput_data.crystal_info[attr])
    println("qeOutputData: ", dfpt_data.properties[attr], "\n")
end

struct Phonons
    function Phonons(
        ebinputfilepath::AbstractString,
        ifc3outputfilepath::AbstractString,
    )
        ebdata = ebInputData(ebinputfilepath)
        type2mass = ebdata.crystal_info["masses"]
        # This will convert the index of an atom in the unitcell to its type. The 
        # type is represented as an integer and can be converted to the mass or the 
        # element/species.
        atindex2type = ebdata.crystal_info["atomtypes"]
        lattvecs = ebdata.crystal_info["lattvecs"]

        # Third Order Force Constants Data
        todata = Ifc3Output(ifc3outputfilepath)
        ifc3_tensor = todata.properties["ifc3_tensor"]
        numtriplets = todata.properties["numtriplets"]
        # Convert the triplet index to needed quantities
        trip2position_j = todata.properties["trip2position_j"]
        trip2position_k = todata.properties["trip2position_k"]
        trip2atomindeces = todata.properties["trip2atomindices"]

        # No force constants without the mass prefactor
        for itrip in 1:numtriplets
            ifc3_tensor[:, :, :, itrip] ./= begin
                sqrt(
                    type2mass[atindex2type[trip2atomindeces[1, itrip]]] *
                    type2mass[atindex2type[trip2atomindeces[2, itrip]]] *
                    type2mass[atindex2type[trip2atomindeces[3, itrip]]],
                )
            end
        end
    end
end

"""
Snapping positions such that their coordinates are the result of an integer linear
combination integer multiples of basisvectors from a given basis.
"""
function snap_to_basis!(basis::Matrix{Float64}, positions::Matrix{Float64})
    # Solve a system of linear equations to get the coefficients that make up the 
    # positions as linear combinations
    positions_frac = \(basis, positions)
end

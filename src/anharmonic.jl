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
        # Positions in Angstrom
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

        # Snapping the cell origin coordinates of the two displaced atoms to the 
        # exact lattice vectors. Converted to the units of lattvecs (nm)
        snap_to_lattvecs!(lattvecs, trip2position_k ./ 10)
        snap_to_lattvecs!(lattvecs, trip2position_j ./ 10)
    end
end

"""
Snapping positions such that their coordinates are the result of an linear
combination of integer multiples of basisvectors from given lattice vectors.

Lattice Vectors have to be in the form `lattvecs = [ a1 a2 a3 ]`
"""
function snap_to_lattvecs!(
    lattvecs::Matrix{Float64},
    positions::Matrix{Float64},
)::Matrix{Float64}
    # Solve a system of linear equations to get the coefficients that make up the 
    # positions as linear combinations

    # Getting the positions in fractional coordinates as integers
    # Rounding like fortrans `anint()`. Ties are rounded away from zero
    positions_frac = round.(\(lattvecs, positions), RoundNearestTiesAway)

    # Overriding the original positions
    positions = lattvecs * positions_frac

    return
end

"""
Calculate the q-vectors in crystal coordinates for the third phonon in both 3-phonon
processes using the conservation of momentum. `q1` is intended to be sampled along
a symmetry path and `q2` comes from the full Brillouin zone

For the absorption it must be that `q3 = (q1 + q2) mod G` and for the emission it
must be that `q3 = (q1 - q2) mod G` where `mod G` wraps beyond the Brillouin Zone
boundaries.
"""
function fill_q3!(
    q1_cryst::Matrix{Float64},
    q2_cryst::Matrix{Float64},
    q3_emission::Matrix{Float64},
    q3_absorption::Matrix{Float64},
)::Matrix{Float64}
    iq3 = 0
    for iq1 in axes(q1_cryst)
        for iq2 in axes(q2_cryst)
            iq3 += 1
            q3_emission[:, iq3] =
                mod.(
                    view(q1_cryst, :, iq1) + view(q2_cryst, :, iq2),
                    ones(Float64, 3),
                )

            q3_absorption[:, iq3] =
                mod.(
                    view(q1_cryst, :, iq1) - view(q2_cryst, :, iq2),
                    ones(Float64, 3),
                )
        end
    end

    return
end

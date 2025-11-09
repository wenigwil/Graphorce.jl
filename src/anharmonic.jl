struct Phonons
    function Phonons(
        ebdata::ebInputData,
        deconvolution::DeconvData,
        sodata::qeIfc2Output,
        todata::Ifc3Output,
        q1_cryst::Matrix{Float64};
        bz_sampling::Tuple{Int64,Int64,Int64} = (30, 30, 30),
    )
        # System description
        numatoms = ebdata.allocations["numatoms"]
        type2mass = ebdata.crystal_info["masses"]
        atindex2type = ebdata.crystal_info["atomtypes"]
        # The lattice vectors are in nm when they come out of the input.nml
        lattvecs = ebdata.crystal_info["lattvecs"]
        reclattvecs = calc_reciprocal_lattvecs(lattvecs)

        # Third Order Force Constants Data
        ifc3_tensor = todata.properties["ifc3_tensor"]
        numtriplets = todata.properties["numtriplets"]
        # Extract the data to convert the triplet index to
        # Positions in Angstrom
        trip2position_j = todata.properties["trip2position_j"]
        trip2position_k = todata.properties["trip2position_k"]
        # Atomic Indices of the whole triplet
        trip2atomindeces = todata.properties["trip2atomindices"]

        # We will include a mass normalization into the ifc3
        @info "Mass-normalizing the ifc3..."
        for itrip in 1:numtriplets
            ifc3_tensor[:, :, :, itrip] ./= begin
                sqrt(
                    type2mass[atindex2type[trip2atomindeces[1, itrip]]] *
                    type2mass[atindex2type[trip2atomindeces[2, itrip]]] *
                    type2mass[atindex2type[trip2atomindeces[3, itrip]]],
                )
            end
        end

        @info "Snapping read-in positions to direct lattice grid..."
        # The phonopy ifc3-file gives us the cell coordinates of the two displaced 
        # atoms. We make sure these are EXACTLY (numerics, huh) on our grid defined 
        # by the lattice vectors.
        snap_to_lattvecs!(lattvecs, trip2position_k ./ 10)
        snap_to_lattvecs!(lattvecs, trip2position_j ./ 10)

        @info "Converting the supplied q1 to cartesian coordinates..."
        # Get the q1 in cartesian coordinates
        q1_cart = q1_cryst * permutedims(reclattvecs)
        numq1 = size(q1_cart, 1)

        @info "Building uniformly sampled Brioullin Zone..."
        # Build the q2 by sampling the Brillouin Zone first in crystal coordinates
        # and then converting to cartesian
        q2_cryst = sample_cube(bz_sampling)
        numq2 = size(q2_cryst, 1)

        @info "Converting the calculated q2 to cartesian coordinates..."
        # q2_cryst[i,:] yields the i-th q2, so that we have to transpose
        q2_cart = q2_cryst * permutedims(reclattvecs)

        @info """
        Building q3s for emission and absorption process...
            """ numq1 numq2 (numq1 * numq2)
        # Now we will compute the both sets of q3 for absorption and emission
        # For the absorption we have q3 = ( q1 + q2 ) mod G
        # For the emission we have q3 = ( q1 - q2 ) mod G
        # We will do the computation in crystal coordinates so folding becomes easy
        q3_emission = Matrix{Float64}(undef, (numq1 * numq2, 3))
        q3_absorption = Matrix{Float64}(undef, (numq1 * numq2, 3))
        fill_q3!(q1_cryst, q2_cryst, q3_emission, q3_absorption)
        # Convert to cartesian coordinates
        q3_emission = q3_emission * permutedims(reclattvecs)
        q3_absorption = q3_absorption * permutedims(reclattvecs)
        numq3 = numq1 * numq2
        numallq = numq1 + numq2 + 2 * numq1 * numq2
        @info "Total number of qpoints is..." numallq

        # Put it all together for eigenvector calculation
        # Stacked vertically (dim = 1) in order q1 then q2 then q3_e then q3_a
        allq = vcat(q1_cart, q2_cart, q3_emission, q3_absorption)

        @info "Calculating all eigenvectors and frequencies..."
        harmonic = LatticeVibrations(ebdata, sodata, deconvolution, allq)

        # Split it apart again reshape
        q1_freqs = view(harmonic.fullq_freqs, 1:numq1, :)
        q1_freqs = reshape(q1_freqs, (numq1 * 3 * numatoms))
        q2_freqs = view(harmonic.fullq_freqs, (numq1 + 1):numq2, :)
        q2_freqs = reshape(q2_freqs, (numq2 * 3 * numatoms))
        q3_emission_freqs = view(harmonic.fullq_freqs, (numq2 + 1):numq3, :)
        q3_emission_freqs = reshape(q3_emission_freqs, (numq3 * 3 * numatoms))
        q3_absorption_freqs = view(harmonic.fullq_freqs, (numq3 + 1):numallq, :)
        q3_absorption_freqs = reshape(q3_absorption_freqs, (numq3 * 3 * numatoms))

        # TODO: Now the eigenvectors are in the form of eigvecs[iq,s,α,k] and we only 
        # have to reshape them into eigvecs[λ,α,k]. This 
        q1_eigvecs = view(harmonic.eigdisplacement, 1:numq1, :, :, :)
        q1_eigvecs = reshape(q1_eigvecs, (numq1 * 3 * numatoms, 3, numatoms))
        q2_eigvecs = view(harmonic.eigdisplacement, (numq1 + 1):numq2, :, :, :)
        q2_eigvecs = reshape(q1_eigvecs, (numq2 * 3 * numatoms, 3, numatoms))

        q3_emission_eigvecs =
            view(harmonic.eigdisplacement, (numq2 + 1):numq3, :, :, :)
        q3_emission_eigvecs =
            reshape(q3_emission_eigvecs, (numq3 * 3 * numatoms, 3, numatoms))

        q3_absorption_eigvecs =
            view(harmonic.eigdisplacement, (numq3 + 1):numallq, :, :, :)
        q3_absorption_eigvecs =
            reshape(q3_absorption_eigvecs, (numq3 * 3 * numatoms, 3, numatoms))
    end
end

"""
Snapping positions such that their coordinates are the result of an linear
combination of integer multiples of basisvectors from given lattice vectors.

Lattice Vectors have to be in the form `lattvecs = [ a1 a2 a3 ]`
"""
function snap_to_lattvecs!(lattvecs::Matrix{Float64}, positions::Matrix{Float64})
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
)
    iq3 = 0
    for iq1 in axes(q1_cryst, 1)
        for iq2 in axes(q2_cryst, 1)
            iq3 += 1
            q3_emission[iq3, :] =
                mod.(
                    view(q1_cryst, iq1, :) + view(q2_cryst, iq2, :),
                    ones(Float64, 3),
                )

            q3_absorption[iq3, :] =
                mod.(
                    view(q1_cryst, iq1, :) - view(q2_cryst, iq2, :),
                    ones(Float64, 3),
                )
        end
    end

    return
end

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
        for itrip in axes(ifc3_tensor, 1)
            ifc3_tensor[itrip, :, :, :] ./= begin
                sqrt(
                    type2mass[atindex2type[trip2atomindeces[itrip, 1]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 2]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 3]]],
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

        @info "Building uniformly sampled Brioullin Zone..."
        # Build the q2 by sampling the Brillouin Zone first in crystal coordinates
        # and then converting to cartesian
        q2_cryst = sample_cube(bz_sampling)

        numq1 = size(q1_cart, 1)
        numq2 = size(q2_cryst, 1)
        numq12 = numq1 + numq2
        numq3 = numq1 * numq2
        numq123 = numq12 + numq3
        numallq = numq12 + 2 * numq3

        @info "Converting the calculated q2 to cartesian coordinates..."
        # q2_cryst[i,:] yields the i-th q2, so that we have to transpose
        q2_cart = q2_cryst * permutedims(reclattvecs)

        @info """
        Building q3s for emission and absorption process...
            """ numq1 numq2 numq3
        @info "Total number of qpoints is..." numallq
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

        # Put it all together for eigenvector calculation
        # Stacked vertically (dim = 1) in order q1 then q2 then q3_e then q3_a
        allq = vcat(q1_cart, q2_cart, q3_emission, q3_absorption)

        @info "Calculating all eigenvectors and frequencies..."
        harmonic = LatticeVibrations(ebdata, sodata, deconvolution, allq)

        # Split it apart again and reshape
        # Frequency reshaping goes from ω[iq,branch] to ω[λ]
        q1_freqs = view(harmonic.fullq_freqs, 1:numq1, :)
        q1_freqs = reshape(q1_freqs, (numq1 * 3 * numatoms))

        q2_freqs = view(harmonic.fullq_freqs, (numq1 + 1):numq12, :)
        q2_freqs = reshape(q2_freqs, (numq2 * 3 * numatoms))

        q3_emission_freqs = view(harmonic.fullq_freqs, (numq12 + 1):numq123, :)
        q3_emission_freqs = reshape(q3_emission_freqs, (numq3 * 3 * numatoms))

        q3_absorption_freqs = view(harmonic.fullq_freqs, (numq123 + 1):numallq, :)
        q3_absorption_freqs = reshape(q3_absorption_freqs, (numq3 * 3 * numatoms))

        # Eigenvector reshaping goes from eigvecs[iq,branch,α,k] to eigvecs[λ,α,k]
        # where λ conforms with mux2to1(s,iq,numq) from misc.jl 
        q1_eigvecs = view(harmonic.eigdisplacement, 1:numq1, :, :, :)
        q1_eigvecs = reshape(q1_eigvecs, (numq1 * 3 * numatoms, 3, numatoms))

        q2_eigvecs = view(harmonic.eigdisplacement, (numq1 + 1):numq12, :, :, :)
        q2_eigvecs = reshape(q1_eigvecs, (numq2 * 3 * numatoms, 3, numatoms))

        q3_emission_eigvecs =
            view(harmonic.eigdisplacement, (numq12 + 1):numq123, :, :, :)
        q3_emission_eigvecs =
            reshape(q3_emission_eigvecs, (numq3 * 3 * numatoms, 3, numatoms))

        q3_absorption_eigvecs =
            view(harmonic.eigdisplacement, (numq123 + 1):numallq, :, :, :)
        q3_absorption_eigvecs =
            reshape(q3_absorption_eigvecs, (numq3 * 3 * numatoms, 3, numatoms))
    end
end

function calc_Λplus(
    λ::Int64,
    ω::Float64,
    kbT::Float64,
    numatoms::Int64,
    q2::Matrix{Float64},
    q3::Matrix{Float64},
    q2_freqs::Vector{Float64},
    q3_freqs::Vector{Float64},
    q1_eigvecs::Array{ComplexF64,3},
    q2_eigvecs::Array{ComplexF64,3},
    q3_eigvecs::Array{ComplexF64,3},
    ifc3_tensor::Array{Float64,3},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64};
    smearing::Float64,
)
    numq2 = size(q2, 1)
    numq3 = size(q3, 1)

    Λplus = 0.0
    for λ′ in 1:(numq2 * 3 * numatoms)
        for λ′′ in 1:(numq3 * 3 * numatoms)
            ω′ = view(q2_freqs[λ′])
            ω′′ = view(q3_freqs[λ′′])

            statistics = bose(ω′, kbT) - bose(ω′′, kbT)
            Λplus += @views begin
                statistics / (ω′ * ω′′) *
                calc_V2(
                    λ,
                    λ′,
                    λ′′,
                    q2,
                    q3,
                    q1_eigvecs,
                    q2_eigvecs,
                    q3_eigvecs,
                    ifc3_tensor,
                    trip2atomindices,
                    trip2position_j,
                    trip2position_k,
                ) *
                δ(ω, ω′′ - ω′; smearing)
            end
        end
    end
end

"""
Calculate the squared matrix element V for a three-phonon scattering (absorption or
emission). The interaction strengths are indexed by three phonon-state indices λ, λ'
and λ". Each index is the result of multiplexing the q-point index `iq` and the
branch index `s` using the `mux2to1(s,iq,numq)`-function.

# Arguments

  - `λ::Int64`: Collective state index of phonon-1 in the process. Multiplexing
    conforms with the `reshape()`-function and consists of the q-point index `iq` and
    the branch index `s` such that `λ = (s-1)*numq + iq` where `numq` is the maximum
    of `iq`. `iq` is considered to be a fast index and `s` the slow one.
  - `λ′::Int64`: Collective state index of phonon-2 in the process. Same as `λ`
  - `λ′′::Int64`: Collective state index of phonon-3 in the process. Same as `λ`
  - `q2::Matrix{Float64}`: q-points that are the result of sampling the whole
    Brillouin zone. `q2[iq,:]` will yield the `iq`-th q-point. Demuxing `λ′` into its
    slow and fast indices by using `demux1to2()` will yield `iq` to index this list
    of q-points.
  - `q3::Matrix{Float64}`: List of q-points as `q2`. This list must be generated from
    the states of phonon-1 and phonon-2 by utilizing the momentum conservation. Thus
    they (with their harmonic `q3_eigvecs`) determine whether the `calc_V2()`
    calculates the matrix element for the absorption or emission process.
  - `q1_eigvecs::Array{ComplexF64,3}`: Output from `LatticeVibrations()` upon input
    with the q-points of phonon-1.
  - `q2_eigvecs::Array{ComplexF64,3}`: Output from `LatticeVibrations()` upon input
    with the q-points of phonon-2.
  - `q3_eigvecs::Array{ComplexF64,3}`: Output from `LatticeVibrations()` upon input
    with the q-points of phonon-3.
  - `ifc3_tensor::Array{Float64,3}`: Ifc3-tensor generated by ShengBTEs thirdorder.py
    indexed by a triplet-index `itrip` and three cartesian indices `α`,`α′` and `α′′`.
  - `trip2atomindices::Matrix{Int64}`: Translation of `itrip` to the atomindex of
    both displaced atoms in the calculation of the `ifc3_tensor`.
  - `trip2position_j::Matrix{Float64}`: Unitcell positions of one of the two
    displaced atoms in the calculation of `ifc3_tensor`. Produced by thirdorder.py
    but snapped to the direct lattice grid for exact matching.
  - `trip2position_k::Matrix{Float64}`: Unitcell positions of one of the two
    displaced atoms in the calculation of `ifc3_tensor`. Produced by thirdorder.py
    but snapped to the direct lattice grid for exact matching.

# Output

  - `V2::ComplexF64`: Matrix element.
"""
function calc_V2(
    λ::Int64,
    λ′::Int64,
    λ′′::Int64,
    q2::Matrix{Float64},
    q3::Matrix{Float64},
    q1_eigvecs::Array{ComplexF64,3},
    q2_eigvecs::Array{ComplexF64,3},
    q3_eigvecs::Array{ComplexF64,3},
    ifc3_tensor::Array{Float64,3},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    # Get a view of the relevant eigenvectors as W_λ[cart_index, atom_index]
    W_λ = view(q1_eigvecs, λ, :, :)
    W_λ′ = view(q2_eigvecs, λ′, :, :)
    W_λ′′ = view(q3_eigvecs, λ′′, :, :)

    # Get the q-points corresponding to λ′ and λ′′
    q′′ = view(q3, demux1to2(λ′′, size(W_λ′′, 1))[2], :)
    q′ = view(q2, demux1to2(λ′, size(W_λ′, 1))[2], :)

    V = 0.0 + im * 0.0
    for α in axes(W_λ, 1)
        for α′ in axes(W_λ′, 1)
            for α′′ in axes(W_λ′′, 1)
                for itrip in axes(ifc3_tensor, 1)
                    V += @views begin
                        W_λ[α, trip2atomindices[itrip, 1]] *
                        W_λ′[α′, trip2atomindices[itrip, 2]] *
                        W_λ′′[α′′, trip2atomindices[itrip, 3]] *
                        ifc3_tensor[itrip, α, α′, α′′] *
                        exp(im * LinAlg.dot(q′, trip2position_j[itrip, :])) *
                        exp(im * LinAlg.dot(q′′, trip2position_k[itrip, :]))
                    end
                end
            end
        end
    end

    V2 = V * conj(V)
    return V2
end

"""
Snapping positions such that their coordinates are the result of an linear
combination of integer multiples of basisvectors from given lattice vectors. It is
assumed that `positions[i,:]` will yield the `i`-th position.

Lattice Vectors have to be in the form `lattvecs = [ a1 a2 a3 ]`
"""
function snap_to_lattvecs!(lattvecs::Matrix{Float64}, positions::Matrix{Float64})
    # Solve a system of linear equations to get the coefficients that make up the 
    # positions as linear combinations

    # Getting the positions in fractional coordinates as integers
    # Rounding like fortrans `anint()`. Ties are rounded away from zero
    positions_frac =
        round.(\(lattvecs, permutedims(positions)), RoundNearestTiesAway)

    # Overriding the original positions
    positions = permutedims(lattvecs * positions_frac)

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

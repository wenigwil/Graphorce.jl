struct Phonons
    function Phonons(
        ebdata::ebInputData,
        deconvolution::DeconvData,
        sodata::qeIfc2Output,
        todata::Ifc3Output,
        q1_cryst::Matrix{Float64},
        cont_freqs::Vector{Float64},
        kbT::Float64,
        smearing::Float64;
        brillouin_sampling::Tuple{Int64,Int64,Int64} = (30, 30, 30),
    )
        # System description
        numatoms = ebdata.allocations["numatoms"]
        type2mass = ebdata.crystal_info["masses"]
        atindex2type = ebdata.crystal_info["atomtypes"]
        # The lattice vectors are in nm when they come out of the input.nml but we 
        # will do everything in Angstrom around here!
        lattvecs = ebdata.crystal_info["lattvecs"] * 10
        reclattvecs = calc_reciprocal_lattvecs(lattvecs)

        # Third Order Force Constants Data
        ifc3_tensor = todata.properties["ifc3_tensor"]
        # Extract the data to convert the triplet index to
        # Positions in Angstrom
        trip2position_j = todata.properties["trip2position_j"]
        trip2position_k = todata.properties["trip2position_k"]
        # Atomic Indices of the whole triplet
        trip2atomindeces = todata.properties["trip2atomindices"]

        # We will include a mass normalization into the ifc3
        for itrip in axes(ifc3_tensor, 1)
            ifc3_tensor[itrip, :, :, :] ./= begin
                sqrt(
                    type2mass[atindex2type[trip2atomindeces[itrip, 1]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 2]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 3]]],
                )
            end
        end

        # The phonopy ifc3-file gives us the cell coordinates of the two displaced 
        # atoms. We make sure these are EXACTLY (numerics, huh) on our grid defined 
        # by the lattice vectors.
        snap_to_lattvecs!(lattvecs, trip2position_k)
        snap_to_lattvecs!(lattvecs, trip2position_j)

        @info "Building HarmonicStatesData..."
        # Calculate and reshape the frequencies and eigenvectors of 3 phonons by a 
        # given sampling of the brillouin zone.
        states = HarmonicStatesData(
            ebdata,
            sodata,
            deconvolution,
            q1_cryst;
            brillouin_sampling,
        )

        # Convert needed q-points into cartesian coordinates
        q2_cart = states.q2_cryst * permutedims(reclattvecs)
        q3_abso_cart = Array{Float64,3}(undef, size(states.q3_abso_cryst))
        q3_emit_cart = Array{Float64,3}(undef, size(states.q3_emit_cryst))
        for iq1 in axes(q1_cryst, 1)
            q3_abso_cart[:, :, iq1] .=
                states.q3_abso_cryst[:, :, iq1] * permutedims(reclattvecs)
            q3_emit_cart[:, :, iq1] .=
                states.q3_emit_cryst[:, :, iq1] * permutedims(reclattvecs)
        end

        numq1 = size(q1_cryst, 1)
        numq2 = size(states.q2_cryst, 1)
        numfreq = size(cont_freqs, 1)

        @info "Beginning lifetime calculation..."
        # Calculating the lifetime
        scattering_rate = Array{Float64,3}(undef, (numfreq, numq1, 3 * numatoms))
        for λ in axes(q1_cryst, 1)
            println("At λ=", λ)
            s1, iq1 = demux1to2(λ, numq1)
            ω = states.q1_freqs[λ]
            for ifreq in 1:numfreq
                if isapprox(ω, 0, atol = smearing)
                    scattering_rate[ifreq, iq1, s1]
                end
                println("\t At ifreq=", ifreq)
                scattering_rate[ifreq, iq1, s1] = begin
                    1 / (numq2 * ω) * (
                        calc_Λplus(
                            λ,
                            cont_freqs[ifreq],
                            kbT,
                            smearing,
                            states,
                            q2_cart,
                            q3_abso_cart,
                            ifc3_tensor,
                            trip2atomindeces,
                            trip2position_j,
                            trip2position_k,
                        ) +
                        0.5 * calc_Λminus(
                            λ,
                            cont_freqs[ifreq],
                            kbT,
                            smearing,
                            states,
                            q2_cart,
                            q3_emit_cart,
                            ifc3_tensor,
                            trip2atomindeces,
                            trip2position_j,
                            trip2position_k,
                        )
                    )
                end
                println("\tLifetime here is ", scattering_rate[ifreq, iq1, s1])
            end
        end
    end
end

function calc_Λminus(
    λ::Int64,
    ω::Float64,
    kbT::Float64,
    smearing::Float64,
    states::HarmonicStatesData,
    q2_cart::Matrix{Float64},
    q3_emit_cart::Array{Float64,3},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    numq1 = size(states.q1_cryst, 1)
    numq2 = size(states.q2_cryst, 1)
    _, iq1 = demux1to2(λ, numq1)

    q2_freqs = states.q2_freqs
    q3_emit_freqs = states.q3_emit_freqs

    Λminus = 0.0
    for λ′ in axes(states.q2_evec, 1)
        for λ′′ in axes(states.q2_evec, 1)
            _, iq′ = demux1to2(λ′, numq2)
            _, iq′′ = demux1to2(λ′′, numq2)

            q′ = q2_cart[iq′, :]
            q′′ = q3_emit_cart[iq′′, :, iq1]
            W_λ = states.q1_evec[λ, :, :]
            W_λ′ = states.q2_evec[λ′, :, :]
            W_λ′′ = states.q3_emit_evec[λ′′, :, :, iq1]
            ω′ = q2_freqs[λ′]
            ω′′ = q3_emit_freqs[λ′′, iq1]

            if isapprox(ω′, 0, atol = smearing) || isapprox(ω′, 0, atol = smearing)
                continue
            end
            statistics = begin
                bose(ω′ / RydtoTHz * rydberg_ev, kbT) +
                bose(ω′′ / RydtoTHz * rydberg_ev, kbT) +
                1
            end

            Λminus += begin
                statistics / (ω′′ * ω′) *
                calc_V2(
                    q′,
                    q′′,
                    W_λ,
                    W_λ′,
                    W_λ′′,
                    ifc3_tensor,
                    trip2atomindices,
                    trip2position_j,
                    trip2position_k,
                ) *
                δ(ω, ω′′ + ω′; smearing)
            end
        end
    end

    return Λminus
end

function calc_Λplus(
    λ::Int64,
    ω::Float64,
    kbT::Float64,
    smearing::Float64,
    states::HarmonicStatesData,
    q2_cart::Matrix{Float64},
    q3_abso_cart::Array{Float64,3},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    numq1 = size(states.q1_cryst, 1)
    numq2 = size(states.q2_cryst, 1)
    _, iq1 = demux1to2(λ, numq1)

    q2_freqs = states.q2_freqs
    q3_abso_freqs = states.q3_abso_freqs

    Λplus = 0.0
    for λ′ in axes(states.q2_evec, 1)
        for λ′′ in axes(states.q2_evec, 1)
            _, iq′ = demux1to2(λ′, numq2)
            _, iq′′ = demux1to2(λ′′, numq2)

            q′ = q2_cart[iq′, :]
            q′′ = q3_abso_cart[iq′′, :, iq1]
            W_λ = states.q1_evec[λ, :, :]
            W_λ′ = states.q2_evec[λ′, :, :]
            W_λ′′ = states.q3_abso_evec[λ′′, :, :, iq1]
            ω′ = q2_freqs[λ′]
            ω′′ = q3_abso_freqs[λ′′, iq1]

            if isapprox(ω′, 0, atol = smearing) || isapprox(ω′, 0, atol = smearing)
                continue
            end
            statistics = begin
                bose(ω′ / RydtoTHz * rydberg_ev, kbT) -
                bose(ω′′ / RydtoTHz * rydberg_ev, kbT)
            end

            Λplus += begin
                statistics / (ω′′ * ω′) *
                calc_V2(
                    q′,
                    q′′,
                    W_λ,
                    W_λ′,
                    W_λ′′,
                    ifc3_tensor,
                    trip2atomindices,
                    trip2position_j,
                    trip2position_k,
                ) *
                δ(ω, ω′′ - ω′; smearing)
            end
        end
    end

    return Λplus
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
    q′::Vector{Float64},
    q′′::Vector{Float64},
    W_λ::Matrix{ComplexF64},
    W_λ′::Matrix{ComplexF64},
    W_λ′′::Matrix{ComplexF64},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    V = 0.0 + im * 0.0
    for α in axes(W_λ, 1)
        for α′ in axes(W_λ′, 1)
            for α′′ in axes(W_λ′′, 1)
                for itrip in axes(ifc3_tensor, 1)
                    V += @views begin
                        W_λ[α, trip2atomindices[itrip, 1]] *
                        W_λ′[α′, trip2atomindices[itrip, 2]] *
                        conj(W_λ′′[α′′, trip2atomindices[itrip, 3]]) *
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

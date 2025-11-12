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

        # Calculating the lifetime
        numcontfreqs = size(cont_freqs, 1)
        Γ = Matrix{Float64}(undef, (numcontfreqs, 3 * numatoms * numq1))

        @info "Starting calculation of phonon lifetime..." numcontfreqs, numq1
        for λ in axes(Γ, 2)
            s, iq = demux1to2(λ, numq1)
            q1_cp = q1_cryst[iq, :]
            println("At λ=", λ, ". Corresponds to q_cart=", q1_cp)
            for ifreq in axes(Γ, 1)
                println("\tAt ifreq=", ifreq)
                Γ[ifreq, λ] = begin
                    1 / (numq1 * q1_freqs[λ]) * (calc_Λplus(
                        smearing,
                        λ,
                        cont_freqs[ifreq],
                        kbT,
                        numatoms,
                        q2_cart,
                        q3_absorption,
                        q2_freqs,
                        q3_absorption_freqs,
                        q1_eigvecs,
                        q2_eigvecs,
                        q3_absorption_eigvecs,
                        ifc3_tensor,
                        trip2atomindeces,
                        trip2position_j,
                        trip2position_k;
                    ))
                end
            end
        end
    end
end

function calc_Λplus(
    smearing::Float64,
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
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64};
)
    numq2 = size(q2, 1)
    numq3 = size(q3, 1)

    Λplus = 0.0
    Threads.@threads for λ′ in 1:(numq2 * 3 * numatoms)
        for λ′′ in 1:(numq3 * 3 * numatoms)
            ω′ = q2_freqs[λ′]
            ω′′ = q3_freqs[λ′′]

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
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    # Get a view of the relevant eigenvectors as W_λ[cart_index, atom_index]
    W_λ = view(q1_eigvecs, λ, :, :)
    W_λ′ = view(q2_eigvecs, λ′, :, :)
    W_λ′′ = view(q3_eigvecs, λ′′, :, :)

    # Get the q-points corresponding to λ′ and λ′′
    q′ = view(q2, demux1to2(λ′, size(W_λ′, 1))[2], :)
    q′′ = view(q3, demux1to2(λ′′, size(W_λ′′, 1))[2], :)

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

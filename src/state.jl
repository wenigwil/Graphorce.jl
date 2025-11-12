struct HarmonicStatesData
    q1_cryst::Matrix{Float64}
    q2_cryst::Matrix{Float64}
    q3_abso_cryst::Array{Float64,3}
    q3_emit_cryst::Array{Float64,3}
    q1_freqs::Vector{Float64}
    q2_freqs::Vector{Float64}
    q3_abso_freqs::Matrix{Float64}
    q3_emit_freqs::Matrix{Float64}
    q1_evec::Array{ComplexF64,3}
    q2_evec::Array{ComplexF64,3}
    q3_abso_evec::Array{ComplexF64,4}
    q3_emit_evec::Array{ComplexF64,4}

    function HarmonicStatesData(
        ebdata::ebInputData,
        sodata::qeIfc2Output,
        deconvolution::DeconvData,
        q1_cryst::Matrix{Float64};
        brillouin_sampling::Tuple{Int64,Int64,Int64} = (10, 10, 10),
    )
        numatoms = ebdata.allocations["numatoms"]
        @info "Sampling the Brillouin Zone..."
        q2_cryst = sample_cube(brillouin_sampling)

        numbranch = 3 * numatoms
        numq1 = size(q1_cryst, 1)
        numq2 = size(q2_cryst, 1)
        numq3 = numq1 * numq2
        numallq = numq1 + numq2 + 2 * numq3

        @info "Building q3 points for emission and absorption..." numq1 numq2 numq3
        q3_emit_cryst = Matrix{Float64}(undef, (numq3, 3))
        q3_abso_cryst = Matrix{Float64}(undef, (numq3, 3))
        fill_q3!(q1_cryst, q2_cryst, q3_abso_cryst, q3_emit_cryst)

        # Stack all qpoints to calculate everything in one go
        allq = vcat(q1_cryst, q2_cryst, q3_emit_cryst, q3_abso_cryst)

        @info "Calculating harmonic properties for all q-points..." numq3 numallq
        harmonic = LatticeVibrations(ebdata, sodata, deconvolution, allq)

        q3_emit_cryst = reshape(q3_emit_cryst, (numq2, numq1, 3))
        q3_emit_cryst = permutedims(q3_emit_cryst, (1, 3, 2))

        q3_abso_cryst = reshape(q3_abso_cryst, (numq2, numq1, 3))
        q3_abso_cryst = permutedims(q3_abso_cryst, (1, 3, 2))

        # Define the ends of all qpoint-set data
        # End for q2 data
        numq12 = numq1 + numq2
        # End for q3_abso data
        numq123 = numq12 + numq3

        # Split it apart again
        q1_freqs = harmonic.fullq_freqs[1:numq1, :]
        q2_freqs = harmonic.fullq_freqs[(numq1 + 1):numq12, :]
        q3_emit_freqs = harmonic.fullq_freqs[(numq12 + 1):numq123, :]
        q3_abso_freqs = harmonic.fullq_freqs[(numq123 + 1):numallq, :]

        # FREQUENCY RESHAPING goes from ω[iq,branch] to ω[λ]
        q1_freqs = reshape(q1_freqs, (numq1 * numbranch))
        q2_freqs = reshape(q2_freqs, (numq2 * numbranch))

        # EIGENVECTOR RESHAPING goes from eigvecs[iq,branch,α,k] to eigvecs[λ,α,k]
        # where λ conforms with mux2to1(s,iq,numq)
        q1_eigvecs = harmonic.eigdisplacement[1:numq1, :, :, :]
        q2_eigvecs = harmonic.eigdisplacement[(numq1 + 1):numq12, :, :, :]
        q3_emit_eigvecs = harmonic.eigdisplacement[(numq12 + 1):numq123, :, :, :]
        q3_abso_eigvecs = harmonic.eigdisplacement[(numq123 + 1):numallq, :, :, :]

        q1_eigvecs = reshape(q1_eigvecs, (numq1 * numbranch, 3, numatoms))
        q2_eigvecs = reshape(q2_eigvecs, (numq2 * numbranch, 3, numatoms))

        @info "Beginning reshaping of q3 points..."
        # We need to demux the fastest index in the q3_abso and q3_emit data which 
        # was the result of muxing iq1 (slow) and iq2 (fast). See the docs on 
        # fill_q3!()
        # q3_freqs[iq3, branch] to q3_freqs[λ′, iq1]
        q3_emit_freqs = reshape(q3_emit_freqs, (numq2, numq1, numbranch))
        q3_abso_freqs = reshape(q3_abso_freqs, (numq2, numq1, numbranch))
        q3_emit_freqs = permutedims(q3_emit_freqs, (1, 3, 2))
        q3_abso_freqs = permutedims(q3_abso_freqs, (1, 3, 2))
        q3_emit_freqs = reshape(q3_emit_freqs, (numq2 * numbranch, numq1))
        q3_abso_freqs = reshape(q3_abso_freqs, (numq2 * numbranch, numq1))

        # q3_eigvecs[iq3, branch, α, k] to q3_eigvecs[λ′, α, k, iq1]
        q3_emit_eigvecs =
            reshape(q3_emit_eigvecs, (numq2, numq1, numbranch, 3, numatoms))
        q3_abso_eigvecs =
            reshape(q3_abso_eigvecs, (numq2, numq1, numbranch, 3, numatoms))
        q3_emit_eigvecs = permutedims(q3_emit_eigvecs, (1, 3, 4, 5, 2))
        q3_abso_eigvecs = permutedims(q3_abso_eigvecs, (1, 3, 4, 5, 2))

        q3_emit_eigvecs =
            reshape(q3_emit_eigvecs, (numq2 * numbranch, 3, numatoms, numq1))
        q3_abso_eigvecs =
            reshape(q3_abso_eigvecs, (numq2 * numbranch, 3, numatoms, numq1))

        new(
            q1_cryst,
            q2_cryst,
            q3_abso_cryst,
            q3_emit_cryst,
            q1_freqs,
            q2_freqs,
            q3_abso_freqs,
            q3_emit_freqs,
            q1_eigvecs,
            q2_eigvecs,
            q3_abso_eigvecs,
            q3_emit_eigvecs,
        )
    end
end

"""
Calculate the q-vectors in crystal coordinates for the third phonon in both 3-phonon
processes using the conservation of momentum. `q1` is intended to be sampled along
a symmetry path and `q2` comes from the full Brillouin zone

For the absorption it must be that `q3 = (q1 + q2) mod G` and for the emission it
must be that `q3 = (q1 - q2) mod G` where `mod G` wraps beyond the Brillouin Zone
boundaries.

Both q3-matrices have a fast index `iq3` that is the result of multiplexing `iq1` and
`iq2` with each other. `iq2` is the fast index and `iq1` is the slow index.
"""
function fill_q3!(
    q1_cryst::Matrix{Float64},
    q2_cryst::Matrix{Float64},
    q3_abso::Matrix{Float64},
    q3_emit::Matrix{Float64},
)
    iq3 = 0
    for iq1 in axes(q1_cryst, 1)
        for iq2 in axes(q2_cryst, 1)
            iq3 += 1
            q3_abso[iq3, :] =
                mod.(
                    view(q1_cryst, iq1, :) + view(q2_cryst, iq2, :),
                    ones(Float64, 3),
                )
            q3_emit[iq3, :] =
                mod.(
                    view(q1_cryst, iq1, :) - view(q2_cryst, iq2, :),
                    ones(Float64, 3),
                )
        end
    end

    return
end

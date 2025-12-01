# This file consists of a collection of functions for testing Phunky.jl

"""
This function will test the reshaping that is implemented in the state.jl-file. It
reassures that the call vstack() with the slicing afterwards was done correctly.
"""
function test_anharm_reshape(bz_sampling::Tuple{Int64,Int64,Int64})
    ebdata = ebInputData("examples/input.nml")
    deconvolution = DeconvData(ebdata)
    sodata = qeIfc2Output("examples/espresso.ifc2")
    numatoms = ebdata.allocations["numatoms"]
    numbranches = 3 * numatoms

    point_labels = ["Γ", "K", "L", "U", "W", "X", raw"$\mathrm{W}_2$"]
    seek_path_points = [
        0.0 0.0 0.0
        0.375 0.375 0.75
        0.5 0.5 0.5
        0.625 0.25 0.625
        0.5 0.25 0.75
        0.5 0.0 0.5
        0.75 0.25 0.5
    ]
    route = seek_path_points[[2, 1, 3], :]
    sympath =
        Sympath(seek_path_points, point_labels, route; numpoints_per_section = 25)

    q1_cryst = sympath.qpoints
    numq1 = size(q1_cryst, 1)

    q2_cryst = sample_cube(bz_sampling)
    numq2 = size(q2_cryst, 1)

    numq3 = numq1 * numq2
    q3_abso_cryst = Matrix{Float64}(undef, (numq3, 3))
    q3_emit_cryst = Matrix{Float64}(undef, (numq3, 3))
    fill_q3!(q1_cryst, q2_cryst, q3_abso_cryst, q3_emit_cryst)

    # THIS IS WHERE THE REAL TEST ACTUALLY BEGINS
    # We will build the harmonic data by HarmonicStateData and for each qpoint-set 
    # alone and compare 

    states = HarmonicStatesData(
        ebdata,
        sodata,
        deconvolution,
        q1_cryst,
        brillouin_sampling = bz_sampling,
    )

    q1_data = LatticeVibrations(ebdata, sodata, deconvolution, q1_cryst)
    q2_data = LatticeVibrations(ebdata, sodata, deconvolution, q2_cryst)
    q3_abso_data = LatticeVibrations(ebdata, sodata, deconvolution, q3_abso_cryst)
    q3_emit_data = LatticeVibrations(ebdata, sodata, deconvolution, q3_emit_cryst)

    allq = vcat(q1_cryst, q2_cryst, q3_abso_cryst, q3_emit_cryst)
    allq_data = LatticeVibrations(ebdata, sodata, deconvolution, allq)
    copy_q1_evec = allq_data.eigdisplacement[(1:numq1), :, :, :]

    # copy_q1_data = LatticeVibrations(ebdata, sodata, deconvolution, q1_cryst)
    # copy_q1_evec = copy_q1_data.eigdisplacement

    # First we test the data corresponding to the q1-set coming from Sympath
    maxλ = numq1 * numbranches
    q1freq_match = Vector{Bool}(undef, maxλ)

    println("Checking the data related to the q1-set from Sympath...")
    println("\tChecking the frequencies...")
    for λ in 1:maxλ
        s, iq = demux1to2(λ, numq1)
        q1freq_match[λ] =
            isapprox(q1_data.fullq_freqs[iq, s], states.q1_freqs[λ], atol = 1e-9)
        # println("λ=", λ, "\t s=", s, "\t iq=", iq, "\t check=", q1freq_match[λ])
    end
    println("\tIs q1freq_match true everywhere? ", all(q1freq_match), "\n")

    println("\tChecking the eigenvectors...")
    q1eigvec_match = Array{ComplexF64}(undef, (maxλ, numatoms))
    for λ in 1:maxλ
        s, iq = demux1to2(λ, numq1)
        for k in 1:numatoms
            q1eigvec_match[λ, k] = LinAlg.norm(
                LinAlg.cross(
                    states.q1_evec[λ, :, k],
                    q1_data.eigdisplacement[iq, s, :, k],
                ),
            )
        end
    end
    # println("\tIs q1eigvec_match true everywhere? ", all(q1eigvec_match), "\n")

    println("\tChecking if the slicing is the problem...")
    q1copy_match = Array{Bool}(undef, (numq1, numbranches, numatoms))
    for λ in 1:maxλ
        s, iq = demux1to2(λ, numq1)
        for k in 1:numatoms
            parallel = LinAlg.cross(
                q1_data.eigdisplacement[iq, s, :, k],
                states.q1_evec[λ, :, k],
            )
            q1copy_match[iq, s, k] =
                all(isapprox.(parallel, 0.0 + im * 0.0, atol = 0.25))
            # if ~q1copy_match[iq, s, k]
            #     println("Not parallel at q=", q1_cryst[iq], " ", (iq, s, k))
            #     println(
            #         q1_data.eigdisplacement[iq, s, :, k],
            #         " ≈ ",
            #         copy_q1_evec[iq, s, :, k],
            #     )
            # end
        end
    end
    println(
        "\tIs q1copy_match true everywhere? ",
        sum(q1copy_match) / length(q1copy_match),
        "\n",
    )
    return q1copy_match
end

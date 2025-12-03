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
    # alone with LatticeVibrations() and compare 

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

    # allq = vcat(q1_cryst, q2_cryst, q3_abso_cryst, q3_emit_cryst)
    # allq_data = LatticeVibrations(ebdata, sodata, deconvolution, allq)
    # copy_q1_evec = allq_data.eigdisplacement[(1:numq1), :, :, :]

    # This is a specific identical copy which can be used to test a test.
    copy_q1_data = LatticeVibrations(ebdata, sodata, deconvolution, q1_cryst)

    # First we test the data corresponding to the q1-set coming from Sympath
    maxλ = numq1 * numbranches
    maxλ′ = numq2 * numbranches
    maxλ′′ = numq1 * numq2 * numbranches

    println("Checking the data related to the q1-set from Sympath...")
    println("\tChecking the frequencies...")
    q1freq_match = Vector{Bool}(undef, maxλ)
    for λ in 1:maxλ
        s, iq = demux1to2(λ, numq1)
        q1freq_match[λ] =
            isapprox(q1_data.fullq_freqs[iq, s], states.q1_freqs[λ], atol = 1e-6)
        # println("λ=", λ, "\t s=", s, "\t iq=", iq, "\t check=", q1freq_match[λ])
    end
    println("\tIs q1freq_match true everywhere? ", all(q1freq_match), "\n")

    println("\tChecking for the parallelity of eigenvectors...")
    q1evec_match = Array{Bool}(undef, (numq1, numbranches, numatoms))
    for λ in 1:maxλ
        s, iq = demux1to2(λ, numq1)
        for k in 1:numatoms
            parallel = LinAlg.cross(
                q1_data.eigdisplacement[iq, s, :, k],
                states.q1_evec[λ, :, k],
            )
            q1evec_match[iq, s, k] =
                all(isapprox.(parallel, 0.0 + im * 0.0, atol = 0.1))
        end
    end
    println(
        "\tIs q1evec_match true everywhere? ",
        round(100.0 * sum(q1evec_match) / length(q1evec_match), digits = 2),
        " %\n",
    )

    println("Checking the data related to the q2-set from Full BZ...")
    println("\tChecking the frequencies...")
    q2freq_match = Vector{Bool}(undef, maxλ′)
    for λ′ in 1:maxλ′
        s′, iq′ = demux1to2(λ′, numq2)
        q2freq_match[λ′] =
            isapprox(q2_data.fullq_freqs[iq′, s′], states.q2_freqs[λ′], atol = 1e-6)
    end
    println("\tIs q2freq_match true everywhere? ", all(q2freq_match), "\n")

    println("\tChecking for the parallelity of eigenvectors...")
    q2evec_match = Array{Bool}(undef, (numq2, numbranches, numatoms))
    for λ′ in 1:maxλ′
        s′, iq′ = demux1to2(λ′, numq2)
        for k′ in 1:numatoms
            parallel = LinAlg.cross(
                q2_data.eigdisplacement[iq′, s′, :, k′],
                states.q2_evec[λ′, :, k′],
            )
            q2evec_match[iq′, s′, k′] =
                all(isapprox.(parallel, 0.0 + im * 0.0, atol = 0.1))
        end
    end
    println(
        "\tIs q2evec_match true everywhere? ",
        round(100.0 * sum(q2evec_match) / length(q2evec_match), digits = 2),
        " %\n",
    )

    # Now what's left is checking if the permuting and reshaping worked for the data 
    # for both sets regarding the q3
    println("Checking the data related to the q3-absorption and -emission sets...")
    println("\tChecking the reshaping on the q3-vectors themselves...")

    q3abso_match = Array{Bool}(undef, (numq2, 3, numq1))
    q3emit_match = Array{Bool}(undef, (numq2, 3, numq1))
    for iq2 in 1:numq2
        for α in 1:3
            for iq1 in 1:numq1
                iq3 = mux2to1(iq1, iq2, numq2)
                q3abso_match[iq2, α, iq1] = isapprox(
                    q3_abso_cryst[iq3, α],
                    states.q3_abso_cryst[iq2, α, iq1],
                    atol = 1e-6,
                )

                q3emit_match[iq2, α, iq1] = isapprox(
                    q3_emit_cryst[iq3, α],
                    states.q3_emit_cryst[iq2, α, iq1],
                    atol = 1e-6,
                )
            end
        end
    end
    println("\tIs q3abso_match true everywhere? ", all(q3abso_match))
    println("\tIs q3emit_match true everywhere? ", all(q3emit_match), "\n")

    println("\tChecking the reshaping of the frequencies...")
    q3absofreq_match = Array{Bool}(undef, (maxλ′, numq1))
    q3emitfreq_match = Array{Bool}(undef, (maxλ′, numq1))

    for λ′ in 1:maxλ′
        s′, iq′ = demux1to2(λ′, numq2)
        for iq in 1:numq1
            iq′′ = mux2to1(iq, iq′, numq2)
            q3absofreq_match[λ′, iq] = isapprox(
                q3_abso_data.fullq_freqs[iq′′, s′],
                states.q3_abso_freqs[λ′, iq],
                atol = 1e-6,
            )

            q3emitfreq_match[λ′, iq] = isapprox(
                q3_emit_data.fullq_freqs[iq′′, s′],
                states.q3_emit_freqs[λ′, iq],
                atol = 1e-6,
            )
        end
    end
    println("\tIs q3absofreq_match true everywhere? ", all(q3absofreq_match))
    println("\tIs q3emitfreq_match true everywhere? ", all(q3emitfreq_match), "\n")

    println("\tChecking the reshaping of the q3 eigenvectors...")
    q3absoevec_match = Array{Bool}(undef, (maxλ′, numatoms, numq1))
    q3emitevec_match = Array{Bool}(undef, (maxλ′, numatoms, numq1))
    for λ′ in 1:maxλ′
        s′, iq′ = demux1to2(λ′, numq2)
        for k in 1:numatoms
            for iq in 1:numq1
                iq′′ = mux2to1(iq, iq′, numq2)
                parallel_emit = LinAlg.cross(
                    q3_abso_data.eigdisplacement[iq′′, s′, :, k],
                    states.q3_abso_evec[λ′, :, k, iq],
                )
                q3absoevec_match[λ′, k, iq] =
                    all(isapprox.(parallel_emit, 0.0 + im * 0.0, atol = 0.1))

                parallel_emit = LinAlg.cross(
                    q3_emit_data.eigdisplacement[iq′′, s′, :, k],
                    states.q3_emit_evec[λ′, :, k, iq],
                )
                q3emitevec_match[λ′, k, iq] =
                    all(isapprox.(parallel_emit, 0.0 + im * 0.0, atol = 0.1))
            end
        end
    end
    println(
        "\tIs q3absoevec_match true everywhere? ",
        round(100.0 * sum(q3absoevec_match) / length(q3absoevec_match), digits = 2),
        " %",
    )
    println(
        "\tIs q3emitevec_match true everywhere? ",
        round(100.0 * sum(q3emitevec_match) / length(q3emitevec_match), digits = 2),
        " %\n",
    )
end

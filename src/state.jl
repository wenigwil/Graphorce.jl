struct ThreePhononStates
    function ThreePhononStates(q1_cryst::Matrix{Float64}, q2_cryst::Matrix{Float64})
        # Define some numbers to make future things a little shorter
        numq1 = size(q1_cryst, 1)
        numq2 = size(q2_cryst, 1)
        numq12 = numq1 + numq2
        numq3 = numq1 * numq2
        numq123 = numq12 + numq3
        numallq = numq12 + 2 * numq3

        q3_emission = Matrix{Float64}(undef, (numq1 * numq2, 3))
        q3_absorption = Matrix{Float64}(undef, (numq1 * numq2, 3))

        fill_q3!(q1_cryst, q2_cryst, q3_emission, q3_absorption)
    end
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

"""
    calc_reciprocal_lattvecs(lattvecs::Matrix{Float64})

Calculate the reciprocal lattice vectors from a collection of direct lattice vectors.

# Arguments

  - `lattvecs::Matrix{Float64}`: Collection of direct lattice vectors. `lattvecs[:,i]`
    should yield the `i`-th lattice vector.

# Output

  - `reclattvecs::Matrix{Float64}`: Collection of reciprocal lattice vectors. The
    vectors are stored rowwise like `lattvecs`. The vectors already contain a
    prefactor of `2*pi`
"""
function calc_reciprocal_lattvecs(lattvecs::Matrix{Float64})
    reclattvecs = Matrix{Float64}(undef, size(lattvecs))
    for i in axes(reclattvecs, 1)
        j = mod(i, 3) + 1
        k = mod(j, 3) + 1
        reclattvecs[:, i] = LinAlg.cross(lattvecs[:, j], lattvecs[:, k])
    end
    volume = abs(LinAlg.dot(lattvecs[:, 1], reclattvecs[:, 1]))
    reclattvecs = (2 * pi / volume) * reclattvecs
    return reclattvecs
end

"""
Sample the Brillouin Zone in crystal coordinates that go from 0 to 1. Should be used
to sample the Brillouin Zone in combination with `calc_reciprocal_lattvecs`.

# Output

  - `points_cryst::Matrix{Float64}`: A list of uniformly sampled vectors from a
    1-cube.
"""
function sample_cube(sampling::Tuple{Int64,Int64,Int64})::Matrix{Float64}
    points_cryst = Matrix{Float64}(undef, (3, prod(sampling .+ 1)))

    # Sample a cube by dividing all axis in chunks between 0 and 1
    iq = 0
    for i in 0:sampling[1]
        for j in 0:sampling[2]
            for k in 0:sampling[3]
                iq += 1
                points_cryst[1, iq] = (i / sampling[1])
                points_cryst[2, iq] = (j / sampling[2])
                points_cryst[3, iq] = (k / sampling[3])
            end
        end
    end

    return points_cryst
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

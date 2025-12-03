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
to sample the Brillouin Zone in combination with `calc_reciprocal_lattvecs()` which
will convert the sampled points to cartesian coordinates

# Output

  - `points_cryst::Matrix{Float64}`: A list of uniformly sampled vectors from a
    1-cube. `points_cryst[i,:]` yields the `i`-th vector.
"""
function sample_cube(sampling::Tuple{Int64,Int64,Int64})::Matrix{Float64}
    # Sample a cube by dividing all axis in chunks between 0 and 1
    points_cryst = Matrix{Float64}(undef, (prod(sampling .+ 1), 3))

    iq = 0
    for i in 0:sampling[1]
        for j in 0:sampling[2]
            for k in 0:sampling[3]
                iq += 1
                points_cryst[iq, 1] = (i / sampling[1])
                points_cryst[iq, 2] = (j / sampling[2])
                points_cryst[iq, 3] = (k / sampling[3])
            end
        end
    end

    return points_cryst
end

function sample_cube_Γ_centered(sampling::Tuple{Int64,Int64,Int64})::Matrix{Float64}
    if ~all(iseven.(sampling))
        Logging.@error "sample_cube_Γ_centered: Sampling must be a collection of even integers!"
        return
    end

    # Sample a cube by dividing all axis in chunks between -1 and 1
    points_cryst = Matrix{Float64}(undef, (prod(sampling .+ 1), 3))

    sampling = div.(sampling, 2)
    iq = 0
    for i in (-sampling[1]):sampling[1]
        for j in (-sampling[2]):sampling[2]
            for k in (-sampling[3]):sampling[3]
                iq += 1
                points_cryst[iq, 1] = (i / sampling[1])
                points_cryst[iq, 2] = (j / sampling[2])
                points_cryst[iq, 3] = (k / sampling[3])
            end
        end
    end

    return points_cryst
end

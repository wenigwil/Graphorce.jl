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
    reclattvecs = Matrix{Float64}(undef, (3, 3))
    for i in 1:3
        j = mod(i, 3) + 1
        k = mod(j, 3) + 1
        reclattvecs[:, i] = LinAlg.cross(lattvecs[:, j], lattvecs[:, k])
    end
    volume = abs(LinAlg.dot(lattvecs[:, 1], reclattvecs[:, 1]))
    reclattvecs = (2 * pi / volume) * reclattvecs
    return reclattvecs
end

function read_highsympath(file::AbstractString)
    sympathfile = open(file)

    numlines = parse(Int64, readline(sympathfile))

    sympath = Matrix{Float64}(undef, (numlines, 3))
    for line in 1:numlines
        test = parse.(Float64, split(readline(sympathfile)))
        sympath[line, :] = test
    end

    close(sympathfile)
    return sympath
end

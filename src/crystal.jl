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

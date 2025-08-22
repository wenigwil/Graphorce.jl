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

function points_to_path(pointlist::Matrix{Float64}; numpoints_per_section = 50)
    numpoints = size(pointlist, 1)
    numsections = numpoints - 1
    numpoints_path = numsections * numpoints_per_section + 1

    path = zeros(Float64, (numpoints_path, 3))
    path[1, :] = pointlist[1, :]
    section = 0
    for is in 1:numsections
        start = 1 + numpoints_per_section * section
        stop = 1 + numpoints_per_section * (section + 1)

        direction = pointlist[is + 1, :] - pointlist[is, :]
        direction ./= numpoints_per_section
        i = 0
        for js in start:stop
            path[js, :] .= i * direction + path[start, :]
            i += 1
        end
        section += 1
    end

    return path
end

function path_to_distance(path::Matrix{Float64})
    numpoints = size(path, 1)

    distances = Vector{Float64}(undef, numpoints)

    distances[1] = 0.0
    for i in 1:(numpoints - 1)
        dist = LinAlg.norm(path[i, :] - path[i + 1, :])

        distances[i + 1] = distances[i] + dist
    end

    return distances
end
